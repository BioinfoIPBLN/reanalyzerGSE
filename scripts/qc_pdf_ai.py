#!/usr/bin/env python3
"""
qc_pdf_ai.py -- parse exported QC tables and generate AI commentaries PDF.

Discovers exported TSV tables in QC_and_others/tables/, queries the OpenAI-compatible
LLM endpoint (via llm_common.py), and uses matplotlib.backends.backend_pdf.PdfPages
to generate an annotated QC PDF report (<Project>_<NormType>_QC_commented.pdf).
"""

import argparse
import os
import re
import sys
import textwrap
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Import shared LLM module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common

# Config defaults
DEF_ENDPOINT       = os.environ.get("LLM_ENDPOINT")
DEF_MODEL          = os.environ.get("LLM_MODEL")
DEF_API_KEY        = os.environ.get("LLM_API_KEY", "dummy")
DEF_CONTEXT_WINDOW = int(os.environ.get("LLM_CONTEXT_WINDOW", "128000"))
DEF_TIMEOUT        = int(os.environ.get("LLM_TIMEOUT", "300"))

SYS_PROMPT_GENERIC = "\n".join([
    "You are an expert bioinformatician evaluating RNA-Seq Quality Control (QC) metrics.",
    "Summarise the table data in 2-4 short, clear bullet points.",
    "State key findings, flag potential anomalies or outliers, and indicate whether sample quality is acceptable.",
    "Do not invent facts not present in the data.",
])

SECTION_PROMPTS = {
    "00_sample_targets.tsv": "Analyze sample experimental design: group balance, replicate numbers, and potential confounding metadata.",
    "01_density_metrics.tsv": "Evaluate gene expression count density distribution, recommended cutoffs, and noise filtering efficacy.",
    "02_boxplots_norm.tsv": "Assess expression value distributions across samples after normalization (check median alignment and IQR consistency).",
    "03_library_sizes.tsv": "Compare sequencing library sizes across samples. Highlight any sample with unusually low or high depth.",
    "04_alignment_categories.tsv": "Evaluate alignment percentages (uniquely mapped vs multi-mapped vs unmapped). Highlight low mapping efficiency.",
    "05_spearman_correlation.tsv": "Analyze sample-to-sample Spearman correlation matrix. Identify strongly correlated groups and any low-correlation outliers.",
    "06_mds_coordinates.tsv": "Assess sample clustering in 2D MDS space. State whether samples cluster primarily by biological condition.",
    "07_pca_coordinates.tsv": "Evaluate PCA sample separation along PC1 and PC2 relative to variance explained.",
    "08_top250_variable_genes.tsv": "Comment on top 250 variable genes and expression heterogeneity.",
    "09_top_overrepresented_genes.tsv": "Evaluate top over-represented genes per sample. Flag potential rRNA, globin, or mitochondrial dominance.",
    "10_top100_spearman_corr.tsv": "Assess sample clustering based on the top 100 most variable genes.",
    "11_condition_pairwise_means.tsv": "Evaluate mean expression correlation across experimental conditions.",
    "12_gene_body_coverage.tsv": "Assess 5'-to-3' transcript coverage uniformity across gene bodies. Flag 3' bias or RNA degradation.",
}

def ask_llm(cfg, task_prompt, table_text):
    data_text = llm_common.cap_text(table_text, cfg.llm_context_window)
    messages = [
        {"role": "system", "content": SYS_PROMPT_GENERIC},
        {"role": "user",   "content": f"Task: {task_prompt}\n\nData:\n{data_text}"},
    ]
    try:
        text, usage = llm_common.chat_completion(
            cfg.llm_endpoint, cfg.llm_model, cfg.llm_api_key, messages, timeout=cfg.llm_timeout)
    except llm_common.LLMTimeout:
        raise
    except Exception as e:
        print(f"[QC PDF AI] Error: {e}", flush=True)
        return "_(no AI response)_", {}
    return (text or "_(no AI response)_"), usage

def render_ai_slide(title, ai_text, usage_info, model_name, data_filename=None):
    """Render a clean A4 PDF slide containing the AI commentary using Matplotlib."""
    fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 dimensions in inches
    ax.axis("off")

    # Outer decorative card box
    rect = plt.Rectangle((0.05, 0.05), 0.9, 0.9, transform=ax.transAxes,
                         facecolor="#f8f9fa", edgecolor="#007bff", linewidth=2, zorder=1)
    ax.add_patch(rect)

    # Title header
    ax.text(0.1, 0.88, f"🤖 AI Summary: {title}", transform=ax.transAxes,
            fontsize=16, fontweight="bold", color="#007bff", va="top")

    # Body text wrapping
    wrapped_lines = []
    for paragraph in ai_text.split("\n"):
        if paragraph.strip():
            wrapped_lines.extend(textwrap.wrap(paragraph.strip(), width=75))
            wrapped_lines.append("")
        else:
            wrapped_lines.append("")

    body_content = "\n".join(wrapped_lines)
    ax.text(0.1, 0.82, body_content, transform=ax.transAxes,
            fontsize=11, color="#333333", va="top", linespacing=1.4)

    # Footer provenance info
    secs = float(usage_info.get("duration_s", 0.0) or 0.0)
    toks = llm_common.answer_tokens(usage_info)
    tps = (toks / secs) if secs > 0 else 0.0
    stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())

    table_str = f" — Table parsed: {os.path.basename(data_filename)}" if data_filename else ""
    footer_text = (f"AI summary — Automatic interpretation by an LLM; always verify against the data\n"
                   f"Generated by {model_name} on {stamp}, {tps:.0f} TPS, {secs:.1f} seconds{table_str}")
    ax.text(0.1, 0.08, footer_text, transform=ax.transAxes,
            fontsize=8, color="#6c757d", style="italic", va="bottom")

    plt.tight_layout()
    return fig

def merge_pdfs_if_possible(orig_pdf, ai_pdf_pages_file, output_pdf):
    """Attempt to interleave original plot PDF and AI PDF pages using pypdf/PyPDF2."""
    try:
        try:
            import pypdf
            reader_orig = pypdf.PdfReader(orig_pdf)
            reader_ai = pypdf.PdfReader(ai_pdf_pages_file)
            writer = pypdf.PdfWriter()

            n_orig = len(reader_orig.pages)
            n_ai = len(reader_ai.pages)

            for i in range(max(n_orig, n_ai)):
                if i < n_orig:
                    writer.add_page(reader_orig.pages[i])
                if i < n_ai:
                    writer.add_page(reader_ai.pages[i])

            with open(output_pdf, "wb") as fh:
                writer.write(fh)
            print(f"[QC PDF AI] Successfully interleaved plot & AI pages into: {output_pdf}", flush=True)
            return True
        except ImportError:
            import PyPDF2
            reader_orig = PyPDF2.PdfReader(orig_pdf)
            reader_ai = PyPDF2.PdfReader(ai_pdf_pages_file)
            writer = PyPDF2.PdfWriter()

            n_orig = len(reader_orig.pages)
            n_ai = len(reader_ai.pages)

            for i in range(max(n_orig, n_ai)):
                if i < n_orig:
                    writer.add_page(reader_orig.pages[i])
                if i < n_ai:
                    writer.add_page(reader_ai.pages[i])

            with open(output_pdf, "wb") as fh:
                writer.write(fh)
            print(f"[QC PDF AI] Successfully interleaved plot & AI pages into: {output_pdf}", flush=True)
            return True
    except Exception as e:
        print(f"[QC PDF AI] Note: PDF interleaving unavailable ({e}). AI commentary PDF saved as standalone.", flush=True)
        return False

def main():
    p = argparse.ArgumentParser(description="QC Figures AI PDF generator.")
    p.add_argument("--tables-dir", required=True, help="Directory containing exported TSV tables.")
    p.add_argument("--pdf", default="", help="Original QC PDF file to pair with (optional).")
    p.add_argument("--llm-endpoint", default=DEF_ENDPOINT, help="OpenAI-compatible URL")
    p.add_argument("--llm-model", default=DEF_MODEL, help="Model name")
    p.add_argument("--llm-api-key", default=DEF_API_KEY, help="API key")
    p.add_argument("--llm-context-window", type=int, default=DEF_CONTEXT_WINDOW)
    p.add_argument("--llm-timeout", type=int, default=DEF_TIMEOUT)
    args = p.parse_args()

    if not args.llm_endpoint or not args.llm_model:
        print("[QC PDF AI] LLM_ENDPOINT and LLM_MODEL must be supplied. Exiting without changes.", flush=True)
        sys.exit(0)

    cfg = args
    tables_dir = os.path.abspath(args.tables_dir)

    if not os.path.isdir(tables_dir):
        print(f"[QC PDF AI] Tables directory not found: {tables_dir}", flush=True)
        sys.exit(0)

    tsv_files = sorted([f for f in os.listdir(tables_dir) if f.endswith(".tsv") or f.endswith(".txt")])
    if not tsv_files:
        print(f"[QC PDF AI] No exported tables found in {tables_dir}", flush=True)
        sys.exit(0)

    parent_dir = os.path.dirname(tables_dir)
    out_pdf_standalone = os.path.join(parent_dir, "QC_AI_commentary_slides.pdf")

    secrets = llm_common.secret_values(endpoint=cfg.llm_endpoint, api_key=cfg.llm_api_key)

    with PdfPages(out_pdf_standalone) as pdf:
        for tsv_name in tsv_files:
            tsv_path = os.path.join(tables_dir, tsv_name)
            task_prompt = SECTION_PROMPTS.get(tsv_name, SYS_PROMPT_GENERIC)
            title = tsv_name.replace(".tsv", "").replace(".txt", "").replace("_", " ").title()

            try:
                with open(tsv_path, encoding="utf-8", errors="replace") as fh:
                    table_text = fh.read()
            except OSError:
                continue

            if not table_text.strip():
                continue

            print(f"[QC PDF AI] Generating commentary for section: {tsv_name}...", flush=True)
            try:
                ai_text, usage = ask_llm(cfg, task_prompt, table_text)
                if llm_common.is_degraded(ai_text):
                    print(f"[QC PDF AI] Degraded response dropped for {tsv_name}", flush=True)
                    continue

                fig = render_ai_slide(title, ai_text, usage, cfg.llm_model, data_filename=tsv_name)
                pdf.savefig(fig)
                plt.close(fig)
            except llm_common.LLMTimeout:
                print(f"[QC PDF AI] Timeout generating AI commentary for {tsv_name}", flush=True)

    llm_common.scrub_file(out_pdf_standalone, secrets)
    print(f"[QC PDF AI] Generated AI commentary slides PDF: {out_pdf_standalone}", flush=True)

    # Interleave with original PDF if provided
    if args.pdf and os.path.isfile(args.pdf):
        commented_pdf = args.pdf.replace(".pdf", "_commented.pdf")
        if merge_pdfs_if_possible(args.pdf, out_pdf_standalone, commented_pdf):
            llm_common.scrub_file(commented_pdf, secrets)

if __name__ == "__main__":
    main()
