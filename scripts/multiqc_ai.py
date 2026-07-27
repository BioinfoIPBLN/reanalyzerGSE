#!/usr/bin/env python3
"""
multiqc_ai.py -- run MultiQC and add AI summaries, standalone (no sushi/ezRun).

Two independent modes, both work for ANY MultiQC module (not just FastQC):

  A) --builtin   : use MultiQC's own AI feature (one global report summary,
                   plus per-section summary buttons) via CLI flags.
  B) --per-section : custom pass. After MultiQC runs, discover every section
                   in multiqc_report.html, send each section's data file to an
                   OpenAI-compatible chat endpoint, and inject a pre-baked
                   summary inline under each section heading.

You can pass both.

The LLM is any OpenAI-compatible /v1/chat/completions endpoint (a local vLLM /
Ollama / llama.cpp server, or api.openai.com, etc.). No endpoint or model is
baked in -- you must supply them at runtime via CLI flags or env vars:
    LLM_ENDPOINT   OpenAI-compatible /v1/chat/completions URL   (or --llm-endpoint)
    LLM_MODEL      model name                                   (or --llm-model)
    LLM_API_KEY    Bearer token; "dummy" for servers that ignore it (or --llm-api-key)

Usage:
    export LLM_ENDPOINT=https://your-llm-host/v1/chat/completions
    export LLM_MODEL=your-model-name
    export LLM_API_KEY=dummy
    python multiqc_ai.py --analysis-dir ./results --out-dir ./multiqc --per-section

    # or entirely on the command line:
    python multiqc_ai.py --analysis-dir ./results --out-dir ./multiqc --per-section \
        --llm-endpoint https://your-llm-host/v1/chat/completions --llm-model your-model-name
    python multiqc_ai.py --analysis-dir ./results --out-dir ./multiqc --builtin
"""
import argparse, os, re, subprocess, sys, time

# Import the sibling shared module regardless of how this script is invoked
# (PATH, absolute path, symlink): put its own directory on sys.path first.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common

# ------------------------------------------------------- LLM config defaults
# Precedence for every parameter:  CLI flag  >  env var  >  fallback below.
# NOTE: endpoint and model have NO baked-in fallback on purpose -- no site-
# specific host or model name is hardcoded anywhere. They must be supplied at
# runtime (flag or env var) whenever an AI mode is used; main() enforces this.
# The remaining fallbacks are all generic/neutral, not tied to any deployment.
DEF_ENDPOINT       = os.environ.get("LLM_ENDPOINT")            # required at runtime
DEF_MODEL          = os.environ.get("LLM_MODEL")              # required at runtime
DEF_API_KEY        = os.environ.get("LLM_API_KEY", "dummy")
DEF_CONTEXT_WINDOW = int(os.environ.get("LLM_CONTEXT_WINDOW", "128000"))
DEF_PROVIDER       = os.environ.get("LLM_PROVIDER", "custom")   # MultiQC --ai-provider
DEF_TIMEOUT        = int(os.environ.get("LLM_TIMEOUT", "300"))   # 5 min/response cap

DEF_SYS_PROMPT = "\n".join([
    "You are an expert in bioinformatics, sequencing, and QC reports.",
    "You are given the data for a single section of a MultiQC report (a plot, table, or heatmap).",
    "Produce 1-2 concise bullet points summarising the key observations and any QC issues.",
    "If nothing is concerning, say so in one bullet.",
    "Use markdown. Highlight severity with directives like :span[39.2%]{.text-red}, .text-orange, .text-yellow, .text-green.",
    "Highlight sample names with :sample[name]{.text-red} etc.",
    "Use 4 spaces to indent nested lists. Do not add headers.",
])

# Secret redaction + the text-only, sequential LLM call now live in llm_common
# (imported above) and are shared with llm_insight.py.

# ------------------------------------------------------------- MultiQC step
def run_multiqc(cfg, analysis_dir, out_dir, builtin):
    tmp = os.path.join(os.getcwd(), ".multiqc_tmp")
    os.makedirs(tmp, exist_ok=True)
    env = dict(os.environ, TMPDIR=tmp)          # keep MultiQC's tempdir somewhere we own
    cmd = ["multiqc", "--outdir", out_dir, "-f", analysis_dir]
    if builtin:
        # MultiQC needs a key present even for a custom endpoint.
        env["OPENAI_API_KEY"] = cfg.api_key
        summary_flag = "--ai-summary-full" if cfg.ai_summary == "full" else "--ai-summary"
        cmd += [
            summary_flag,                       # bake summary into the HTML
            "--ai-provider", cfg.provider,
            "--ai-model", cfg.model,
        ]
        # endpoint/context-window only apply to the "custom" provider
        if cfg.provider == "custom":
            cmd += [
                "--ai-custom-endpoint", cfg.endpoint,
                "--ai-custom-context-window", str(cfg.context_window),
            ]
    # Mask the endpoint before logging the command line (it carries
    # --ai-custom-endpoint <url> in builtin mode); the model is left visible.
    print("[MultiQC]", llm_common.mask(" ".join(cmd),
          llm_common.secret_values(endpoint=cfg.endpoint, api_key=cfg.api_key))[0], flush=True)
    subprocess.run(cmd, check=True, env=env)

# --------------------------------------------------------- LLM call helper
def ask_llm(cfg, section_id, data):
    """One section summary. Sequential + text-only via llm_common. A timeout is
    re-raised (the caller then marks this and the remaining sections as timed
    out); any other error is swallowed so one bad section never aborts the pass."""
    data = llm_common.cap_text(data, cfg.context_window)
    messages = [
        {"role": "system", "content": cfg.sys_prompt},
        {"role": "user",   "content": f"Section: {section_id}\n\nData:\n{data}"},
    ]
    try:
        text, usage = llm_common.chat_completion(
            cfg.endpoint, cfg.model, cfg.api_key, messages, timeout=cfg.timeout)
    except llm_common.LLMTimeout:
        raise
    except Exception as e:
        print(f"[AI] error on {section_id}: {e}", flush=True)
        return "_(no AI response)_", {}
    return (text or "_(no AI response)_"), usage

# Shared AI-box stamp: model + when + tokens-per-second + seconds, plus the
# standard "verify against the data" wording (same text as the report insights).
AI_TIMEOUT_MSG = "ai_insights timeout. Please try again"

def ai_footer(cfg, usage, data_filename=None):
    secs = float(usage.get("duration_s", 0.0) or 0.0)
    toks = llm_common.answer_tokens(usage)
    tps = (toks / secs) if secs > 0 else 0.0
    stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    table_str = f" — Table parsed: {os.path.basename(data_filename)}" if data_filename else ""
    return (f'<div class="text-muted" style="font-size:.85em;margin-top:8px;'
            f'border-top:1px solid #eee;padding-top:4px;">'
            f'\U0001f916 AI summary — Automatic interpretation by an LLM; always verify '
            f'against the data — Generated by {cfg.model} on {stamp}, '
            f'{tps:.0f} TPS, {secs:.0f} seconds{table_str}</div>')

def timeout_box():
    return (f'<div class="text-muted" style="font-size:.9em;padding:6px 0;">'
            f'\U0001f916 {AI_TIMEOUT_MSG}</div>')

# ------------------------------------------------ markdown -> minimal HTML
def style_directives(txt):
    # Convert markdown bold (**text** or __text__) to <strong>
    txt = re.sub(r"\*\*(.*?)\*\*", r"<strong>\1</strong>", txt)
    txt = re.sub(r"__(.*?)__", r"<strong>\1</strong>", txt)
    # Convert markdown italic (*text*) to <em>
    txt = re.sub(r"(?<!\*)\*(?!\*)(.*?)(?<!\*)\*(?!\*)", r"<em>\1</em>", txt)
    # Convert inline code (`text`) to <code>
    txt = re.sub(r"`([^`]+)`", r"<code>\1</code>", txt)

    for sev in ("red", "orange", "yellow", "green"):
        txt = re.sub(rf":span\[([^\]]*)\]\{{\.text-{sev}\}}", rf'<span class="text-{sev}">\1</span>', txt)
        txt = re.sub(rf":sample\[([^\]]*)\]\{{\.text-{sev}\}}",
                     rf'<span class="text-{sev}" style="font-weight:600;font-style:italic;">\1</span>', txt)
    txt = re.sub(r":span\[([^\]]*)\]\{[^}]*\}", r"\1", txt)
    txt = re.sub(r":sample\[([^\]]*)\]\{[^}]*\}",
                 r'<span style="font-weight:600;font-style:italic;">\1</span>', txt)
    return txt

def md_to_html(txt):
    txt = style_directives(txt)
    out, depth = [], 0
    for line in txt.split("\n"):
        indent = len(line) - len(line.lstrip())
        t = line.strip()
        m = re.match(r"[-*] (.*)", t)
        if m:
            target = indent // 4 + 1
            while depth < target: out.append("<ul>"); depth += 1
            while depth > target: out.append("</ul>"); depth -= 1
            out.append(f"<li>{m.group(1)}</li>")
        elif t:
            while depth > 0: out.append("</ul>"); depth -= 1
            out.append(f"<p>{t}</p>")
    while depth > 0: out.append("</ul>"); depth -= 1
    return "\n".join(out)

# ------------------------------------------ discover sections FROM the HTML
def discover_sections(html):
    """Return the set of MultiQC section anchor ids actually present in the
    report. Generic: works for every module, not just FastQC."""
    ids = set(re.findall(r'id="mqc-section-wrapper-([^"]+)"', html))
    ids |= set(re.findall(r'id="([^"]+)_ai_summary_response"', html))
    return ids

def data_file_for_section(data_dir, sec_id):
    """Best-effort map a MultiQC HTML section id to its exported data file.
    The file name often differs from the section id in case (featureCounts vs
    featurecounts), punctuation (-/_) and plural/singular (assignment vs
    assignments); some sections also carry a generic '<module>-section-N' id
    that encodes no file name at all (e.g. Salmon). Match tolerantly."""
    if sec_id == "general_stats_table":
        p = os.path.join(data_dir, "multiqc_general_stats.txt")
        return p if os.path.exists(p) else None
    try:
        files = [f for f in os.listdir(data_dir) if f.lower().endswith(".txt")]
    except OSError:
        return None

    def first_match(patterns):
        for pat in patterns:
            rx = re.compile(pat + r"$", re.IGNORECASE)
            hits = sorted(f for f in files if rx.match(f))
            if hits:
                return os.path.join(data_dir, hits[0])
        return None

    # Case-insensitive, punctuation- (-/_) and plural-tolerant spellings of the id.
    variants = set()
    for base in (sec_id, sec_id.replace("-", "_"), sec_id.replace("_", "-")):
        variants.add(base)
        variants.add(base[:-1] if base.endswith("s") else base + "s")
    precise = []
    for v in sorted(variants):
        e = re.escape(v)
        precise += [rf"{e}[-_]plot.*\.txt", rf"{e}[-_]table\.txt",
                    rf"{e}.*-heatmap\.txt", rf"{e}\.txt"]
    hit = first_match(precise)
    if hit:
        return hit

    # Check key distinguishing descriptors in sec_id (e.g. length, duplication, adapter)
    # to avoid mis-mapping when the module has multiple distinct plot files.
    key_descriptors = ["length", "duplication", "adapter", "overrepresented", "gc"]
    sec_low = sec_id.lower()
    required_terms = [d for d in key_descriptors if d in sec_low]
    if "n_content" in sec_low or "per_base_n" in sec_low:
        required_terms.append("n_content")

    # Fuzzy token fallback for multi-token / hyphenated ids
    skip = ("multiqc_citations", "multiqc_sources", "multiqc_software_versions",
            "multiqc_data_sources", "section_prompts", "llms-full", "multiqc_general_stats")
    toks = [t for t in re.split(r"[-_]", sec_id.lower())
            if t and t != "section" and not t.isdigit()]
    module = toks[0] if toks else ""
    if not toks or module == "multiqc":
        return None

    need = 2 if len(toks) >= 2 else 1
    best = None
    for f in files:
        low = f.lower()
        if module not in low or low.rsplit(".", 1)[0] in skip:
            continue
        if any(req not in low for req in required_terms):
            continue
        n = sum(1 for t in toks if t in low)
        if n < need:
            continue
        if re.search(r"(plot|heatmap)", low):      kind = 3
        elif "table" in low:                        kind = 2
        elif low.startswith("multiqc_" + module):   kind = 1
        else:                                        kind = 0
        key = (n, kind, -len(low))
        if best is None or key > best[0]:
            best = (key, f)
    return os.path.join(data_dir, best[1]) if best else None

# --------------------------------------------------- inject one summary
def inject(html, sec_id, summary_html):
    # 1) MultiQC's own pre-created empty AI div (present when its AI feature is on)
    empty = rf'<div class="ai-summary-response" id="{re.escape(sec_id)}_ai_summary_response"[^>]*></div>'
    if re.search(empty, html):
        html = re.sub(empty,
            f'<div class="ai-summary-response" id="{sec_id}_ai_summary_response" '
            f'style="margin-bottom:-5px;">{summary_html}</div>', html, count=1)
        # its wrapper ships hidden; reveal it since we pre-baked the answer
        html = re.sub(rf'(id="{re.escape(sec_id)}_ai_summary_wrapper"[^>]*style=")display:\s*none;',
                      r"\1display: block;", html, count=1)
        return html, True
    # 2) fallback: inject right after the section wrapper opening tag
    wrap = rf'(<div [^>]*id="mqc-section-wrapper-{re.escape(sec_id)}"[^>]*>)'
    if re.search(wrap, html):
        block = (f'\\1\n<div class="ai-summary-response" style="margin:10px 0;padding:8px 12px;'
                 f'border-left:3px solid #4a90d9;background:#f5f9fd;">{summary_html}</div>')
        return re.sub(wrap, block, html, count=1), True
    return html, False

# --------------------------------------------------------------- per-section
def per_section(cfg, out_dir):
    html_path = os.path.join(out_dir, "multiqc_report.html")
    data_dir  = os.path.join(out_dir, "multiqc_data")
    html = open(html_path, encoding="utf-8").read()
    sec_ids = discover_sections(html)
    print(f"[AI] {len(sec_ids)} sections found in report", flush=True)

    mappings_log = ["Section ID -> Matched Data File Mappings:", "=" * 60]
    for sec_id in sorted(sec_ids):
        f = data_file_for_section(data_dir, sec_id)
        status = os.path.basename(f) if f else "SKIPPED (no matching data file)"
        mappings_log.append(f"  - {sec_id:<45} -> {status}")
    print("[AI] Section -> Data File Mappings:\n" + "\n".join(mappings_log[2:]), flush=True)
    open(os.path.join(data_dir, "section_mappings.txt"), "w").write("\n".join(mappings_log))

    prompts_log = [f"Per-section AI prompts -- model: {cfg.model}", ""]
    timed_out = False   # once the LLM times out we stop calling and mark the rest
    for sec_id in sorted(sec_ids):
        f = data_file_for_section(data_dir, sec_id)
        if not f:
            print(f"[AI] no data file for {sec_id}; skipping", flush=True)
            continue
        if timed_out:
            # This AI option (MultiQC per-section) already timed out: every
            # remaining box shows the placeholder rather than a stale summary.
            html, _ = inject(html, sec_id, timeout_box())
            continue
        data = open(f, encoding="utf-8", errors="replace").read()
        try:
            text, usage = ask_llm(cfg, sec_id, data)
        except llm_common.LLMTimeout:
            print(f"[AI] TIMEOUT on {sec_id}; marking this and remaining sections", flush=True)
            timed_out = True
            html, _ = inject(html, sec_id, timeout_box())
            continue
        # Refuse an implausibly long (degraded / hallucinated) answer: omit the box.
        if llm_common.is_degraded(text):
            print(f"[AI] {sec_id}: response too long ({len(text)} chars); omitting", flush=True)
            continue
        html, ok = inject(html, sec_id, md_to_html(text) + "\n" + ai_footer(cfg, usage, data_filename=f))
        print(f"[AI] {sec_id}: {'injected' if ok else 'NO ANCHOR'} "
              f"({usage.get('duration_s',0):.1f}s, {usage.get('total_tokens',0)} tok)", flush=True)
        prompts_log += ["=" * 76, f"Section: {sec_id}", "=" * 76, "", "[USER DATA FILE]", f, ""]
    open(os.path.join(data_dir, "section_prompts.txt"), "w").write("\n".join(prompts_log))
    open(html_path, "w", encoding="utf-8").write(html)
    print(f"[AI] wrote {html_path}", flush=True)

# --------------------------------------------------------------------- main
class Config:
    """Resolved LLM parameters threaded through the run. Populated in main()
    from CLI flags (which themselves default to env vars, then to constants)."""
    def __init__(self, a):
        self.endpoint       = a.llm_endpoint
        self.model          = a.llm_model
        self.api_key        = a.llm_api_key
        self.context_window = a.llm_context_window
        self.provider       = a.ai_provider
        self.timeout        = a.llm_timeout
        self.ai_summary     = a.ai_summary          # "full" | "short"
        # system prompt: --sys-prompt-file wins over --sys-prompt over default
        if a.sys_prompt_file:
            self.sys_prompt = open(a.sys_prompt_file, encoding="utf-8").read()
        else:
            self.sys_prompt = a.sys_prompt

def main():
    ap = argparse.ArgumentParser(
        description="Run MultiQC and add AI summaries against any OpenAI-compatible endpoint.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # what to do
    ap.add_argument("--analysis-dir", default=".", help="dir MultiQC scans for tool outputs")
    ap.add_argument("--out-dir", default="multiqc", help="MultiQC output dir")
    ap.add_argument("--builtin", action="store_true", help="use MultiQC's own AI summary")
    ap.add_argument("--per-section", action="store_true", help="custom per-section injection")
    ap.add_argument("--skip-multiqc", action="store_true", help="report already exists; only inject")

    # LLM parameters -- CLI flag > env var > fallback (endpoint/model have none)
    g = ap.add_argument_group("LLM parameters (supply endpoint/model at runtime)")
    g.add_argument("--llm-endpoint", default=DEF_ENDPOINT,
                   help="REQUIRED for AI: OpenAI-compatible /v1/chat/completions URL (env LLM_ENDPOINT)")
    g.add_argument("--llm-model", default=DEF_MODEL,
                   help="REQUIRED for AI: model name (env LLM_MODEL)")
    g.add_argument("--llm-api-key", default=DEF_API_KEY,
                   help="Bearer token; 'dummy' for local servers (env LLM_API_KEY)")
    g.add_argument("--llm-context-window", type=int, default=DEF_CONTEXT_WINDOW,
                   help="context window in tokens (env LLM_CONTEXT_WINDOW)")
    g.add_argument("--llm-timeout", type=int, default=DEF_TIMEOUT,
                   help="per-request timeout, seconds (env LLM_TIMEOUT)")
    g.add_argument("--ai-provider", default=DEF_PROVIDER,
                   help="MultiQC --ai-provider: custom|openai|anthropic|seqera (env LLM_PROVIDER)")
    g.add_argument("--ai-summary", choices=("full", "short"), default="full",
                   help="MultiQC built-in summary length (--builtin only)")
    g.add_argument("--sys-prompt", default=DEF_SYS_PROMPT, help="system prompt for per-section calls")
    g.add_argument("--sys-prompt-file", help="read system prompt from a file (wins over --sys-prompt)")

    a = ap.parse_args()

    # An AI mode needs a model, and (for the custom provider) an endpoint.
    # Nothing is hardcoded, so fail loudly rather than dialing an unset host.
    if a.builtin or a.per_section:
        missing = []
        if not a.llm_model:
            missing.append("--llm-model (or env LLM_MODEL)")
        needs_endpoint = a.per_section or (a.builtin and a.ai_provider == "custom")
        if needs_endpoint and not a.llm_endpoint:
            missing.append("--llm-endpoint (or env LLM_ENDPOINT)")
        if missing:
            ap.error("AI mode requested but not configured; set: " + ", ".join(missing))

    cfg = Config(a)
    ai_mode = a.builtin or a.per_section
    if not a.skip_multiqc:
        run_multiqc(cfg, a.analysis_dir, a.out_dir, a.builtin)
    # Scrub the sensitive endpoint that MultiQC's builtin AI baked into the
    # report NOW -- before the (potentially long, one-LLM-call-per-section)
    # per-section loop -- so the raw endpoint sits on disk for seconds, not for
    # the whole generation. This also makes redaction crash-safe: if per_section
    # dies mid-run, the endpoint is already gone.
    if ai_mode:
        llm_common.redact_outputs(a.out_dir, endpoint=cfg.endpoint, api_key=cfg.api_key)
    if a.per_section:
        per_section(cfg, a.out_dir)
        # per_section writes only the model (kept) + QC data, never the endpoint;
        # scrub once more so the guarantee holds regardless of future changes.
        llm_common.redact_outputs(a.out_dir, endpoint=cfg.endpoint, api_key=cfg.api_key)

if __name__ == "__main__":
    main()
