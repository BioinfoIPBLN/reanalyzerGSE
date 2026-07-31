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

def clean_markdown_text(text):
    """Clean markdown bold/italic asterisks and formatting markers for Matplotlib rendering."""
    # Convert bold/italic asterisks (***text*** or **text** or *text*)
    text = re.sub(r"\*{1,3}(.*?)\*{1,3}", r"\1", text)
    # Convert __bold__ or _italic_
    text = re.sub(r"_{1,2}(.*?)_{1,2}", r"\1", text)
    # Convert inline backticks (`code`)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    # Remove emoji characters that render as unformatted boxes in standard PDF fonts
    text = re.sub(r"[^\x00-\x7F\u00A0-\u024F\u2022]", "", text)
    return text

def render_ai_slide(title, ai_text, usage_info, model_name, data_filename=None):
    """Render a clean horizontal (landscape) A4 PDF slide containing the AI commentary using Matplotlib."""
    fig, ax = plt.subplots(figsize=(11.69, 8.27))  # A4 landscape dimensions in inches
    ax.axis("off")

    # Outer decorative card box (landscape margins)
    rect = plt.Rectangle((0.04, 0.05), 0.92, 0.90, transform=ax.transAxes,
                         facecolor="#f8f9fa", edgecolor="#007bff", linewidth=2, zorder=1)
    ax.add_patch(rect)

    # Clean title and text
    clean_title = clean_markdown_text(title).strip()

    # Title header
    ax.text(0.08, 0.88, f"AI Summary: {clean_title}", transform=ax.transAxes,
            fontsize=16, fontweight="bold", color="#007bff", va="top")

    # Body text wrapping
    cleaned_ai_text = clean_markdown_text(ai_text)
    wrapped_lines = []
    for paragraph in cleaned_ai_text.split("\n"):
        p = paragraph.strip()
        if not p:
            wrapped_lines.append("")
            continue
        # Format bullet points nicely
        if p.startswith("- ") or p.startswith("* "):
            bullet_p = "• " + p[2:].strip()
            wrapped_lines.extend(textwrap.wrap(bullet_p, width=95, subsequent_indent="  "))
        else:
            wrapped_lines.extend(textwrap.wrap(p, width=95))
        wrapped_lines.append("")

    body_content = "\n".join(wrapped_lines)
    ax.text(0.08, 0.80, body_content, transform=ax.transAxes,
            fontsize=11, color="#333333", va="top", linespacing=1.4)

    # Footer provenance info
    secs = float(usage_info.get("duration_s", 0.0) or 0.0)
    toks = llm_common.answer_tokens(usage_info)
    tps = (toks / secs) if secs > 0 else 0.0
    stamp = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())

    table_str = f" — Table parsed: {os.path.basename(data_filename)}" if data_filename else ""
    footer_text = (f"AI summary — Automatic interpretation by an LLM; MUST always verify against the data\n"
                   f"Generated by {model_name} on {stamp}, {tps:.0f} TPS, {secs:.1f} seconds{table_str}")
    ax.text(0.08, 0.08, footer_text, transform=ax.transAxes,
            fontsize=8, color="#6c757d", style="italic", va="bottom")

    plt.tight_layout()
    return fig

def merge_pdfs_if_possible(orig_pdf, ai_pdf_pages_file, output_pdf, tsv_to_slide_idx=None):
    """Attempt to insert original plot PDF and AI PDF pages mapped correctly using pypdf, PyPDF2, or fitz."""
    # Try pypdf with section-aware mapping
    try:
        import pypdf
        reader_orig = pypdf.PdfReader(orig_pdf)
        reader_ai = pypdf.PdfReader(ai_pdf_pages_file)
        writer = pypdf.PdfWriter()

        n_orig = len(reader_orig.pages)
        n_ai = len(reader_ai.pages)

        # Detect section starting pages by searching PDF page text for section numbers/keywords
        section_start_pages = {}
        for idx, page in enumerate(reader_orig.pages):
            txt = page.extract_text() or ""
            # Match titles like "Sample table", "Density plots", "Boxplots", "Library size", "PCA", "Scatter plots"
            for sec_idx, key in enumerate(tsv_to_slide_idx.keys() if tsv_to_slide_idx else []):
                clean_key = key.replace(".tsv", "").replace(".txt", "").split("_", 1)[-1]
                if clean_key.lower() in txt.lower() and sec_idx not in section_start_pages:
                    section_start_pages[sec_idx] = idx

        # Sort sections by page order
        ai_insert_map = {}
        if section_start_pages and tsv_to_slide_idx:
            sorted_secs = sorted(section_start_pages.items(), key=lambda x: x[1])
            for i, (sec_idx, start_p) in enumerate(sorted_secs):
                # Put AI slide right before the NEXT section starts, or at end for last section
                end_p = sorted_secs[i+1][1] - 1 if i + 1 < len(sorted_secs) else n_orig - 1
                ai_insert_map[end_p] = tsv_to_slide_idx.get(list(tsv_to_slide_idx.keys())[sec_idx])

        # Write pages and insert mapped AI slides
        ai_written = set()
        for i in range(n_orig):
            writer.add_page(reader_orig.pages[i])
            if i in ai_insert_map and ai_insert_map[i] is not None:
                slide_idx = ai_insert_map[i]
                if slide_idx < n_ai and slide_idx not in ai_written:
                    writer.add_page(reader_ai.pages[slide_idx])
                    ai_written.add(slide_idx)

        # Append any remaining AI slides at end
        for s_idx in range(n_ai):
            if s_idx not in ai_written:
                writer.add_page(reader_ai.pages[s_idx])

        with open(output_pdf, "wb") as fh:
            writer.write(fh)
        llm_common.log(f"[QC PDF AI] Successfully interleaved plot & AI pages (pypdf section-mapped) into: {output_pdf}")
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] pypdf mapping error: {e}")
        pass

    # Fallback fitz section-mapped
    try:
        import fitz
        doc_orig = fitz.open(orig_pdf)
        doc_ai = fitz.open(ai_pdf_pages_file)
        doc_out = fitz.open()

        n_orig = len(doc_orig)
        n_ai = len(doc_ai)

        section_start_pages = {}
        for idx in range(n_orig):
            txt = doc_orig[idx].get_text() or ""
            for sec_idx, key in enumerate(tsv_to_slide_idx.keys() if tsv_to_slide_idx else []):
                clean_key = key.replace(".tsv", "").replace(".txt", "").split("_", 1)[-1]
                if clean_key.lower() in txt.lower() and sec_idx not in section_start_pages:
                    section_start_pages[sec_idx] = idx

        ai_insert_map = {}
        if section_start_pages and tsv_to_slide_idx:
            sorted_secs = sorted(section_start_pages.items(), key=lambda x: x[1])
            for i, (sec_idx, start_p) in enumerate(sorted_secs):
                end_p = sorted_secs[i+1][1] - 1 if i + 1 < len(sorted_secs) else n_orig - 1
                ai_insert_map[end_p] = tsv_to_slide_idx.get(list(tsv_to_slide_idx.keys())[sec_idx])

        ai_written = set()
        for i in range(n_orig):
            doc_out.insert_pdf(doc_orig, from_page=i, to_page=i)
            if i in ai_insert_map and ai_insert_map[i] is not None:
                slide_idx = ai_insert_map[i]
                if slide_idx < n_ai and slide_idx not in ai_written:
                    doc_out.insert_pdf(doc_ai, from_page=slide_idx, to_page=slide_idx)
                    ai_written.add(slide_idx)

        for s_idx in range(n_ai):
            if s_idx not in ai_written:
                doc_out.insert_pdf(doc_ai, from_page=s_idx, to_page=s_idx)

        doc_out.save(output_pdf)
        llm_common.log(f"[QC PDF AI] Successfully interleaved plot & AI pages (fitz section-mapped) into: {output_pdf}")
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] fitz mapping error: {e}")
        pass

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
        llm_common.log("[QC PDF AI] LLM_ENDPOINT and LLM_MODEL must be supplied. Exiting without changes.")
        sys.exit(0)

    cfg = args
    tables_dir = os.path.abspath(args.tables_dir)

    if not os.path.isdir(tables_dir):
        llm_common.log(f"[QC PDF AI] Tables directory not found: {tables_dir}")
        sys.exit(0)

    tsv_files = sorted([f for f in os.listdir(tables_dir) if f.endswith(".tsv") or f.endswith(".txt")])
    if not tsv_files:
        llm_common.log(f"[QC PDF AI] No exported tables found in {tables_dir}")
        sys.exit(0)

    parent_dir = os.path.dirname(tables_dir)
    out_pdf_standalone = os.path.join(parent_dir, "QC_AI_commentary_slides.pdf")

    secrets = llm_common.secret_values(endpoint=cfg.llm_endpoint, api_key=cfg.llm_api_key)

    tsv_to_slide_idx = {}
    slide_counter = 0

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

            llm_common.log(f"[QC PDF AI] Generating commentary for section: {tsv_name}...")
            try:
                ai_text, usage = ask_llm(cfg, task_prompt, table_text)
                if llm_common.is_degraded(ai_text):
                    llm_common.log(f"[QC PDF AI] Degraded response dropped for {tsv_name}")
                    continue

                fig = render_ai_slide(title, ai_text, usage, cfg.llm_model, data_filename=tsv_name)
                pdf.savefig(fig)
                plt.close(fig)
                tsv_to_slide_idx[tsv_name] = slide_counter
                slide_counter += 1
            except llm_common.LLMTimeout:
                llm_common.log(f"[QC PDF AI] Timeout generating AI commentary for {tsv_name}")

    llm_common.scrub_file(out_pdf_standalone, secrets)
    llm_common.log(f"[QC PDF AI] Generated AI commentary slides PDF: {out_pdf_standalone}")

    # Interleave with original PDF if provided
    if args.pdf and os.path.isfile(args.pdf):
        commented_pdf = args.pdf.replace(".pdf", "_commented.pdf")
        if merge_pdfs_if_possible(args.pdf, out_pdf_standalone, commented_pdf, tsv_to_slide_idx=tsv_to_slide_idx):
            llm_common.scrub_file(commented_pdf, secrets)

if __name__ == "__main__":
    main()
