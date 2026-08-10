"""
llm_common.py -- shared LLM plumbing for the reanalyzerGSE AI features
(multiqc_ai.py and llm_insight.py). One place for the OpenAI-compatible
/v1/chat/completions call, context-window capping, and secret redaction, so the
two callers cannot drift apart and the invariants below hold everywhere.

TEXT-ONLY IS A HARD INVARIANT.
    Every message we send must have plain-string content. chat_completion()
    rejects OpenAI multimodal "parts" arrays (image_url, input_audio, ...). If an
    insight seems to need a picture, send the DATA TABLE behind the picture as
    text instead -- the DE table, not the volcano PNG. No image/binary is ever
    sent to the model.

SEQUENTIAL IS A HARD INVARIANT.
    The local LLM is a single server that must not receive concurrent requests.
    chat_completion() serialises every call: an in-process lock (so no threaded
    fan-out) AND a cross-process file lock (flock on this host), so even separate
    wrapper / insight processes -- or a future accidental '&' -- cannot issue two
    queries to the same endpoint at once. Callers must still iterate
    sequentially; never wrap these calls in a thread or process pool.
    Set LLM_NO_LOCK=1 to disable the cross-process lock; LLM_LOCK_FILE=<path> to
    place it explicitly.

Nothing here is deployment-specific: the endpoint/model/key come from the caller
(CLI flag or LLM_* env var); none is hardcoded.
"""
import hashlib, json, os, re, socket, sys, tempfile, threading, time, urllib.error, urllib.request
try:
    import fcntl                       # POSIX advisory file locking; absent on non-Unix
except ImportError:                    # pragma: no cover
    fcntl = None


class LLMTimeout(Exception):
    """Raised by chat_completion() when the model does not respond within the
    per-request timeout. Callers treat this specially: the whole AI option is
    skipped and its boxes are filled with a 'timeout' placeholder rather than
    silently omitted (so the user knows it was tried and can retry)."""


def log(msg, file=sys.stderr):
    """Print log message prefixed with a timestamp [YYYY-MM-DD HH:MM:SS]."""
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", file=file, flush=True)


# An AI insight is meant to be a few sentences. A response far longer than that
# is almost always model degradation / hallucinated filler (long garbage prose),
# so we refuse to show it. Override with LLM_MAX_ANSWER_CHARS (0 disables).
def max_answer_chars():
    try:
        return int(os.environ.get("LLM_MAX_ANSWER_CHARS") or 4000)
    except ValueError:
        return 4000

def is_degraded(text):
    """True if `text` is implausibly long for a short insight (likely a runaway
    / hallucinated response that should not be surfaced)."""
    lim = max_answer_chars()
    return bool(lim) and len(text or "") > lim

def answer_tokens(usage):
    """Best 'tokens generated' figure for a tokens-per-second display: prefer the
    completion tokens, fall back to total, then 0."""
    for k in ("completion_tokens", "total_tokens"):
        v = (usage or {}).get(k)
        if isinstance(v, (int, float)) and v > 0:
            return int(v)
    return 0

# --------------------------------------------------------------- sequential I/O
_INPROC_LOCK = threading.Lock()

def _lock_path(endpoint):
    if os.environ.get("LLM_LOCK_FILE"):
        return os.environ["LLM_LOCK_FILE"]
    # One lock per endpoint host so calls to the SAME local LLM serialise while
    # unrelated endpoints don't block each other. Use a stable base dir (NOT
    # gettempdir(), which follows $TMPDIR and would fragment the lock between
    # processes that set TMPDIR differently).
    key = hashlib.sha1((endpoint or "default").encode()).hexdigest()[:12]
    base = "/tmp" if os.path.isdir("/tmp") and os.access("/tmp", os.W_OK) else tempfile.gettempdir()
    return os.path.join(base, f"reanalyzergse_llm_{key}.lock")

class _Serial:
    """Serialise LLM calls in-process (always) and across processes on this host
    (best-effort flock). flock is released automatically when the fd closes or
    the process dies, so a crash never leaves a stale lock and cannot deadlock."""
    def __init__(self, endpoint):
        self.endpoint = endpoint
        self.fd = None
    def __enter__(self):
        _INPROC_LOCK.acquire()
        if fcntl is not None and not os.environ.get("LLM_NO_LOCK"):
            try:
                self.fd = os.open(_lock_path(self.endpoint), os.O_CREAT | os.O_RDWR, 0o600)
                fcntl.flock(self.fd, fcntl.LOCK_EX)
            except OSError:
                if self.fd is not None:
                    os.close(self.fd)
                self.fd = None
        return self
    def __exit__(self, *exc):
        if self.fd is not None:
            try:
                fcntl.flock(self.fd, fcntl.LOCK_UN)
            except OSError:
                pass
            os.close(self.fd)
            self.fd = None
        _INPROC_LOCK.release()
        return False

# --------------------------------------------------------------- text-only chat
def _ensure_text_only(messages):
    """Enforce the text-in/text-out contract: content must be a plain string.
    Reject list/dict 'parts' payloads so no caller can smuggle multimodal input
    (images, audio) into what is meant to be a text-only pipeline."""
    for msg in messages:
        if not isinstance(msg.get("content"), str):
            raise ValueError(
                "llm_common: non-text (multimodal) message content is not allowed; "
                "send the underlying data table as text instead.")

def cap_text(text, context_window, floor=0):
    """Trim text to roughly the model's context budget (~3 chars/token), keeping
    at least `floor` chars. Appends a truncation marker when it cuts."""
    limit = max(floor, (context_window or 0) * 3)
    if limit and len(text) > limit:
        return text[:limit] + "\n[... truncated ...]"
    return text

_RETRY_BACKOFFS = [10, 30]  # seconds to wait before retry 2 and 3
_RETRYABLE_HTTP_CODES = {429, 500, 502, 503, 504}

def chat_completion(endpoint, model, api_key, messages, timeout=600):
    """One TEXT-ONLY, SERIALISED chat call with retry.  Returns (text, usage)
    where usage carries the server's token counts plus 'duration_s' (model time,
    excluding any time spent waiting for the sequential lock).

    Retries up to 3 times total on transient failures (timeouts, HTTP 429/5xx,
    connection errors) with exponential backoff (10 s, 30 s).  Permanent errors
    (HTTP 400/401/403) are raised immediately."""
    _ensure_text_only(messages)
    max_attempts = 1 + len(_RETRY_BACKOFFS)  # 3 total
    last_exc = None
    for attempt in range(max_attempts):
        req = urllib.request.Request(
            endpoint, data=json.dumps({"model": model, "messages": messages}).encode(),
            headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"})
        with _Serial(endpoint):                        # <= no two LLM queries ever overlap
            t0 = time.time()
            try:
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    body = json.loads(r.read())
                dur = time.time() - t0
                usage = dict(body.get("usage") or {})
                usage["duration_s"] = dur
                text = (body.get("choices") or [{}])[0].get("message", {}).get("content") or ""
                return text, usage
            except (socket.timeout, TimeoutError) as e:
                last_exc = LLMTimeout(f"no response within {timeout}s (attempt {attempt+1}/{max_attempts})")
                last_exc.__cause__ = e
            # HTTPError MUST be caught before URLError: it is a subclass of it, so the
            # broader clause would shadow this one and no HTTP status would ever be retried
            # (HTTPError.reason is a plain string, matching none of the isinstance checks below).
            except urllib.error.HTTPError as e:
                if e.code in _RETRYABLE_HTTP_CODES:
                    last_exc = e  # retryable server error
                else:
                    raise  # 400/401/403 etc. — permanent
            except urllib.error.URLError as e:
                reason = getattr(e, "reason", None)
                if isinstance(reason, (socket.timeout, TimeoutError)):
                    last_exc = LLMTimeout(f"no response within {timeout}s (attempt {attempt+1}/{max_attempts})")
                    last_exc.__cause__ = e
                elif isinstance(reason, ConnectionRefusedError):
                    last_exc = e  # retryable: server temporarily down
                else:
                    raise  # non-retryable URLError
        # Decide whether to retry
        if attempt < max_attempts - 1:
            wait = _RETRY_BACKOFFS[attempt]
            log(f"[LLM] Retryable error (attempt {attempt+1}/{max_attempts}): "
                f"{last_exc}. Waiting {wait}s before retry...")
            time.sleep(wait)
        else:
            break
    # All attempts exhausted
    if isinstance(last_exc, LLMTimeout):
        raise last_exc
    if last_exc is not None:
        raise last_exc
    raise RuntimeError("chat_completion: unexpected exit from retry loop")

# --------------------------------------------------------------- secret redaction
# The LLM endpoint (and API key) are deployment-sensitive and must not survive in
# any artifact left on disk. The MODEL NAME is intentionally kept -- useful
# provenance, not sensitive. redact_outputs() is called by the MultiQC wrapper
# immediately after MultiQC (before the slow per-section loop) + once more after.
REDACT = "[redacted]"

def secret_values(endpoint=None, api_key=None):
    """Literal strings to scrub, longest-first (so a longer secret is replaced
    before any substring of it). Endpoint + its bare host:port + API key; never
    the model. 'dummy' is the neutral api-key placeholder, not a secret."""
    vals = set()
    for v in (endpoint, api_key):
        if v and v != "dummy" and len(v) >= 4:
            vals.add(v)
    if endpoint:
        m = re.match(r"^\w+://([^/]+)", endpoint)   # bare host:port, in case only that leaks
        if m:
            vals.add(m.group(1))
    return sorted(vals, key=len, reverse=True)

def mask(text, secrets):
    """Replace each secret in `text` with REDACT; return (text, n_replacements)."""
    n = 0
    for s in secrets:
        if s and s in text:
            n += text.count(s)
            text = text.replace(s, REDACT)
    return text, n

def scrub_file(path, secrets):
    # surrogateescape round-trips any non-UTF8 bytes losslessly, so the big
    # self-contained HTML (base64 blobs etc.) is never corrupted -- only the
    # ASCII endpoint/key literals change.
    try:
        with open(path, encoding="utf-8", errors="surrogateescape", newline="") as fh:
            s = fh.read()
    except OSError:
        return 0
    new, n = mask(s, secrets)
    if n:
        try:
            with open(path, "w", encoding="utf-8", errors="surrogateescape", newline="") as fh:
                fh.write(new)
        except OSError:
            return 0
    return n

def redact_outputs(out_dir, endpoint=None, api_key=None):
    """Scrub the endpoint (its bare host:port too) and the API key from every
    text artifact MultiQC produced (the report HTML and everything under
    multiqc_data/). The model name is intentionally left in place."""
    secrets = secret_values(endpoint=endpoint, api_key=api_key)
    if not secrets:
        return
    exts = (".html", ".htm", ".json", ".txt", ".tsv", ".csv", ".log", ".yaml", ".yml", ".md")
    data_dir = os.path.join(out_dir, "multiqc_data")
    live_log = os.path.abspath(os.path.join(out_dir, "multiqc.log"))  # our stdout target: scrub LAST
    targets = []
    html = os.path.join(out_dir, "multiqc_report.html")
    if os.path.isfile(html):
        targets.append(html)
    if os.path.isdir(data_dir):
        for root, _dirs, files in os.walk(data_dir):
            for f in files:
                if f.lower().endswith(exts):
                    targets.append(os.path.join(root, f))
    targets = [t for t in targets if os.path.abspath(t) != live_log]
    files_hit = total = 0
    for path in targets:
        n = scrub_file(path, secrets)
        if n:
            files_hit += 1
            total += n
    log(f"[redact] AI endpoint scrubbed from report artifacts "
        f"({total} occurrence(s) in {files_hit} file(s)).", file=sys.stdout)
    # The live redirect target (out_dir/multiqc.log) last; flush our own output
    # first so it is on disk and itself gets scrubbed. Nothing is printed after.
    sys.stdout.flush()
    if os.path.isfile(live_log) and live_log.lower().endswith(exts):
        scrub_file(live_log, secrets)
