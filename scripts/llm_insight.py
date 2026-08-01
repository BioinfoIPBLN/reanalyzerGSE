#!/usr/bin/env python3
"""
llm_insight.py - generate a short, human-readable AI insight for one pipeline
artifact (a DE table, an enrichment result, sample metadata, ...) and write it
next to the artifact so the HTML reports can embed it.

It is the generic counterpart of multiqc_ai.py: same OpenAI-compatible endpoint,
same LLM_* environment variables, same opt-in behaviour. If no endpoint/model is
configured it exits 0 without writing anything, so callers can invoke it
unconditionally. Any runtime error (network, HTTP, bad response) is also
non-fatal: a warning is printed and nothing is written, so the reports simply
omit the insight box rather than the pipeline failing.

    LLM_ENDPOINT          OpenAI-compatible /v1/chat/completions URL  (or --endpoint)
    LLM_MODEL             model name                                  (or --model)
    LLM_API_KEY           Bearer token ("dummy" for servers that ignore it)
    LLM_CONTEXT_WINDOW    approx token budget, used to cap input size (default 32000)

Usage:
    llm_insight.py --input <file> [<file> ...] --task <dge|enrichment|design|counts|generic> \
                   [--title "Comparison 1"] [--out <file>] [--max-rows N | 0 = all rows]

With several --input files the script runs a sequential MAP-REDUCE: it summarises
each file in its own LLM call (one at a time), then makes a final call that
synthesises those summaries into a single narrative. Every call is text-only and
serialised (only one query hits the local LLM at a time).

Exit codes: 0 = wrote an insight, or nothing to do / non-fatal error (box omitted);
42 = the LLM timed out (>--timeout s). The caller uses 42 to fill that AI option's
boxes with a 'timeout, please try again' placeholder.
"""
import argparse, os, sys, time

# Import the sibling shared module regardless of how this script is invoked.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common

# One shared instruction block, plus a task-specific framing. All prompts insist
# on brevity and on not inventing anything that is not in the supplied data.
COMMON = ("You are a careful bioinformatics assistant. Base every statement only "
          "on the data provided; never invent genes, terms, numbers or conditions. "
          "Write plain prose (no markdown headings, tables or bullet points), a few "
          "short sentences at most. Your text is shown to researchers under an "
          "'AI summary (verify)' label, so hedge biological interpretation. "
          "Always respond in English.")
TASKS = {
    "design": ("Given sample metadata, describe the experiment in 2-4 sentences: the "
               "organism, the groups/conditions being compared, the number of samples "
               "and replicates per group, and the library layout (single/paired) if "
               "evident. State only what the metadata supports."),
    "dge": ("Given a differential-expression table (genes with a log fold-change and an "
            "FDR/adjusted p-value), summarise in 2-4 sentences: how many genes are "
            "significant at FDR < 0.05, the up- vs down-regulated split, and name a few "
            "of the strongest up- and down-regulated genes by fold-change. Do not "
            "speculate on biology beyond naming the genes."),
    "enrichment": ("Given the top over-representation / enrichment results (GO, KEGG, "
                   "Reactome, ... terms with adjusted p-values) for one comparison, write "
                   "3-5 sentences describing the main biological themes among the enriched "
                   "terms and what hypothesis they suggest. Make clear this is a hypothesis "
                   "to be verified, not a conclusion."),
    "counts": ("Given a per-sample summary of expression (log2(TPM+0.1)) for one experiment - "
               "each sample's expression-value range/quartiles and how many genes fall in the "
               "Low/Medium/High expression categories, plus which condition each sample belongs "
               "to - write 2-4 sentences on the overall expression landscape: whether samples "
               "look comparable, whether any sample stands out as an outlier (unusually shifted "
               "range or category split), and whether samples of the same condition look "
               "consistent. Use only the numbers given; do not name individual genes."),
    "generic": "Summarise the following bioinformatics output concisely and factually.",
}

# For map-reduce over MANY tables (used by --task enrichment with several inputs):
# summarise each table on its own (MAP), then synthesise all summaries (REDUCE).
# Both steps stay text-only and are issued as separate SEQUENTIAL queries.
MAP = {
    "enrichment": ("Summarise this SINGLE enrichment result table in 1-2 sentences: name the most "
                   "significant enriched terms and the biological theme they suggest. Use only the "
                   "rows given. If the table is empty or nothing is significant, say so briefly."),
}
REDUCE = {
    "enrichment": ("You are given short summaries of several enrichment result tables (GO, KEGG, "
                   "Reactome, ...) for ONE comparison. Synthesise them into 3-5 sentences describing "
                   "the main biological themes across all tables and the hypothesis they suggest. Make "
                   "clear this is a hypothesis to be verified, not a conclusion. Do not introduce terms "
                   "that are not present in the summaries."),
}

def _map_prompt(task):
    return MAP.get(task, TASKS[task])

def _reduce_prompt(task):
    return REDUCE.get(task,
                      "Synthesise the following per-item summaries into a single concise, factual paragraph.")


def resolve_rows(raw, task, max_rows):
    """Return a bounded slice of the artifact to send. For DGE tables, surface top
    UP and DOWN regulated genes sorted by FDR (balanced in both senses). For enrichment
    tables, sort strictly by FDR. Concatenated multi-file inputs are passed through."""
    lines = raw.splitlines()
    if not lines:
        return raw
    header = lines[0]
    multi = any(l.startswith("### ") for l in lines[:5])
    if task == "dge" and not multi and "\t" in header:
        try:
            cols = [c.lower() for c in header.split("\t")]
            sig = next((i for i, c in enumerate(cols)
                        if any(k in c for k in ("fdr", "padj", "p.adj", "adj.p", "qvalue"))), None)
            if sig is not None:
                body = [l for l in lines[1:] if l.strip()]
                def parse_fdr(l):
                    parts = l.split("\t")
                    try:
                        return (float(parts[sig]), l)
                    except (ValueError, IndexError):
                        return (float("inf"), l)
                parsed = [parse_fdr(l) for l in body]
                # Keep all genes with FDR < 0.05, sorted by FDR ascending
                degs = [p for p in parsed if p[0] < 0.05]
                degs.sort(key=lambda x: x[0])
                if degs:
                    sel = [p[1] for p in degs]
                    note = f"# Data note: {len(degs)} gene(s) met FDR < 0.05 significance."
                else:
                    # Fallback safeguard if no gene reaches FDR < 0.05: take top 100 by FDR
                    parsed.sort(key=lambda x: x[0])
                    sel = [p[1] for p in parsed[:100]]
                    note = "# Data note: No genes met FDR < 0.05 significance. Showing a subset of the top 100 genes by FDR for exploratory trend analysis."
                return f"{note}\n{header}\n" + "\n".join(sel)
        except Exception:
            pass
    elif task == "enrichment" and not multi and "\t" in header:
        try:
            cols = header.split("\t")
            sig = next((i for i, c in enumerate(cols)
                        if any(k in c.lower() for k in ("fdr", "padj", "p.adj", "adj.p", "qvalue"))), None)
            if sig is not None:
                body = [l for l in lines[1:] if l.strip()]
                def key(l):
                    parts = l.split("\t")
                    try:
                        return float(parts[sig])
                    except (ValueError, IndexError):
                        return float("inf")
                body.sort(key=key)
                sel = body if max_rows <= 0 else body[:max_rows]
                return "\n".join([header] + sel)
        except Exception:
            pass
    return "\n".join(lines if max_rows <= 0 else lines[: max_rows + 1])


def _read_capped(path, task, max_rows, context_window):
    """Read one artifact, keep the rows we care about (all rows when max_rows<=0,
    still sorted so the most significant survive if the context window trims)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        raw = fh.read()
    return llm_common.cap_text(resolve_rows(raw, task, max_rows), context_window, floor=2000)


def _call_tps(u):
    """Per-call tokens-per-second from one usage dict (0 if no timing)."""
    sec = float(u.get("duration_s", 0.0) or 0.0)
    return (llm_common.answer_tokens(u) / sec) if sec > 0 else 0.0, sec


def _write_pertable(path, summaries):
    """Persist the per-table MAP summaries as a TSV the report can read:
    table<TAB>tps<TAB>seconds<TAB>summary (summary flattened to one line). No-op
    if `path` is falsy or there is nothing to write."""
    if not path or not summaries:
        return
    def flat(s):
        return " ".join((s or "").split())
    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("table\ttps\tseconds\tsummary\n")
            for name, text, tps, sec in summaries:
                fh.write(f"{name}\t{tps:.0f}\t{sec:.0f}\t{flat(text)}\n")
        llm_common.log(f"[llm_insight] wrote {len(summaries)} per-table note(s) to {path}")
    except OSError as e:
        llm_common.log(f"[llm_insight] could not write per-table notes to {path}: {e}")


def _summarize_one(a, path):
    """Single artifact -> one query -> (summary text, tps, seconds). tps is this
    one call's real generation rate."""
    data = _read_capped(path, a.task, a.max_rows, a.context_window)
    user = (f"{TASKS[a.task]}\n\n"
            + (f"Context: {a.title}\n\n" if a.title else "")
            + f"Data:\n{data}")
    text, u = llm_common.chat_completion(
        a.endpoint, a.model, a.api_key,
        [{"role": "system", "content": COMMON + "\n\n" + TASKS[a.task]},
         {"role": "user", "content": user}], timeout=a.timeout)
    tps, sec = _call_tps(u)
    return (text or "").strip(), tps, sec


def _get_path_label(path, base_dir=None):
    if base_dir and os.path.isabs(path) and os.path.isabs(base_dir):
        try:
            return os.path.relpath(path, base_dir)
        except ValueError:
            pass
    parts = path.replace("\\", "/").split("/")
    if "DGE" in parts:
        idx = parts.index("DGE")
        return "/".join(parts[idx+1:])
    if len(parts) >= 2:
        return "/".join(parts[-2:])
    return os.path.basename(path)


def _map_reduce(a, paths):
    """MANY tables -> one sequential query PER table (all rows), then one query
    that synthesises the per-table summaries into a single narrative. Every call
    goes through llm_common.chat_completion (which serialises), and the loop is
    serial by construction: the next query starts only after the previous returns.
    So no query is ever sent before the previous one has finished. Returns
    (narrative, tps, seconds): tps is the AVERAGE of the individual calls' own
    per-call generation rates (a real per-call speed, not a blended
    total/total), and seconds is the TOTAL model time across every map call +
    the final reduce call."""
    summaries = []    # (table_name, summary_text, tps, seconds) per mapped table
    tps_calls = []     # each call's own tokens-per-second
    seconds = 0.0      # total model time across all calls
    base_dir = getattr(a, "rel_dir", None)
    for i, path in enumerate(paths, 1):
        data = _read_capped(path, a.task, a.max_rows, a.context_window)
        name = _get_path_label(path, base_dir)
        if len([l for l in data.splitlines() if l.strip()]) <= 1:   # header-only / empty table
            llm_common.log(f"[llm_insight] skip {i}/{len(paths)}: {name} (no rows)")
            continue
        user = f"{_map_prompt(a.task)}\n\nTable: {name}\n\nData:\n{data}"
        text, u = llm_common.chat_completion(
            a.endpoint, a.model, a.api_key,
            [{"role": "system", "content": COMMON + "\n\n" + _map_prompt(a.task)},
             {"role": "user", "content": user}], timeout=a.timeout)
        tps, sec = _call_tps(u); seconds += sec
        if tps > 0: tps_calls.append(tps)
        # A single per-table summary should be short; drop a runaway one so it
        # cannot poison the reduce step.
        text = "" if llm_common.is_degraded(text) else (text or "").strip()
        llm_common.log(f"[llm_insight] mapped {i}/{len(paths)}: {name} ({tps:.0f} TPS, {sec:.0f}s)"
                       f"{'' if text else ' (empty)'}")
        if text:
            summaries.append((name, text, tps, sec))
    avg_tps = (sum(tps_calls) / len(tps_calls)) if tps_calls else 0.0
    # Persist the per-table (MAP) summaries so the report can show each one next to
    # its table, not just the synthesised narrative. The calls were made anyway.
    _write_pertable(getattr(a, "pertable_out", None), summaries)
    if not summaries:
        return "", avg_tps, seconds
    joined = "\n".join(f"- {n}: {t}" for n, t, _tp, _s in summaries)
    user = (f"{_reduce_prompt(a.task)}\n\n"
            + (f"Context: {a.title}\n\n" if a.title else "")
            + f"Per-table summaries:\n{joined}")
    text, u = llm_common.chat_completion(
        a.endpoint, a.model, a.api_key,
        [{"role": "system", "content": COMMON + "\n\n" + _reduce_prompt(a.task)},
         {"role": "user", "content": user}], timeout=a.timeout)
    tps, sec = _call_tps(u); seconds += sec
    if tps > 0: tps_calls.append(tps)
    avg_tps = (sum(tps_calls) / len(tps_calls)) if tps_calls else 0.0
    return (text or "").strip(), avg_tps, seconds


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, nargs="+",
                   help="one or more artifact files; multiple => sequential map-reduce")
    p.add_argument("--task", default="generic", choices=list(TASKS))
    p.add_argument("--title", default="")
    p.add_argument("--out", default=None)
    p.add_argument("--pertable-out", default=None,
                   help="map-reduce only: TSV to receive the per-table summaries "
                        "(table, tps, seconds, summary) so the report can show each")
    p.add_argument("--rel-dir", default=None,
                   help="Base directory to format relative paths for per-table labels")
    p.add_argument("--max-rows", type=int, default=40)
    p.add_argument("--endpoint", default=os.environ.get("LLM_ENDPOINT"))
    p.add_argument("--model", default=os.environ.get("LLM_MODEL"))
    p.add_argument("--api-key", default=os.environ.get("LLM_API_KEY", "dummy"))
    p.add_argument("--context-window", type=int,
                   default=int(os.environ.get("LLM_CONTEXT_WINDOW") or 32000))
    p.add_argument("--timeout", type=int, default=300)   # 5 min per LLM response (item: hard cap)
    a = p.parse_args()

    # Opt-in: silently do nothing if no LLM is configured.
    if not a.endpoint or not a.model:
        llm_common.log("[llm_insight] no LLM endpoint/model configured; skipping.")
        return 0
    inputs = [f for f in a.input if os.path.isfile(f)]
    if not inputs:
        llm_common.log(f"[llm_insight] no input file found ({', '.join(a.input)}); skipping.")
        return 0

    out = a.out or (inputs[0] + ".ai_insight.md")
    try:
        t0 = time.time()
        # One artifact -> single query; many -> sequential map-reduce (all serialised).
        text, tps, seconds = (_summarize_one(a, inputs[0]) if len(inputs) == 1
                              else _map_reduce(a, inputs))
        if not text:
            llm_common.log(f"[llm_insight] empty response for {out}; nothing written.")
            return 0
        # Guard against LLM degradation / hallucination: a wildly long "answer"
        # (garbage filler) is refused rather than shown.
        if llm_common.is_degraded(text):
            llm_common.log(f"[llm_insight] response for {out} exceeds {llm_common.max_answer_chars()} chars "
                           f"(likely degraded); nothing written.")
            return 0
        stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
        # Model name is kept (useful provenance); the endpoint is never included.
        # TPS = per-call generation rate (for map-reduce, the average of the
        # individual calls' own rates); seconds = total model time.
        header = (f"**AI summary** — Automatic interpretation by an LLM; MUST always verify against "
                  f"the data — Generated by `{a.model}` on {stamp}, "
                  f"{tps:.0f} TPS, {seconds:.0f} seconds")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(header + "\n\n" + text.strip() + "\n")
        llm_common.log(f"[llm_insight] wrote {out} ({time.time()-t0:.1f}s, {len(inputs)} table(s), "
                       f"{tps:.0f} TPS avg, {seconds:.0f}s total)")
    except llm_common.LLMTimeout as e:
        # Signal the orchestrator (reanalyzerGSE.sh) that THIS ai_insights option
        # timed out; it fills every box for the option with the timeout placeholder.
        llm_common.log(f"[llm_insight] TIMEOUT for {out}: {e}")
        return 42
    except Exception as e:
        llm_common.log(f"[llm_insight] error for {out}: {e}; skipping.")
        return 0
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
