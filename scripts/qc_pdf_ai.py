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
    "Always respond in English.",
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
    "09_dendrogram_distances.tsv": "Analyze sample-to-sample Euclidean distance matrix and hierarchical clustering structure. Identify closely clustered sample groups and early-branching outliers.",
    "10_top_overrepresented_genes.tsv": "Evaluate top over-represented genes per sample. Flag potential rRNA, globin, or mitochondrial dominance.",
    "11_top100_spearman_corr.tsv": "Assess sample clustering based on the top 100 most variable genes.",
    "12_condition_pairwise_means.tsv": "Evaluate mean expression correlation across experimental conditions.",
    "13_gene_body_coverage.tsv": "Assess 5'-to-3' transcript coverage uniformity across gene bodies. Flag 3' bias or RNA degradation.",
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
    """Clean markdown bold/italic asterisks, LLM directives, and formatting markers for Matplotlib rendering."""
    # Convert bold/italic asterisks (***text*** or **text** or *text*)
    text = re.sub(r"\*{1,3}(.*?)\*{1,3}", r"\1", text)
    # Convert __bold__ or _italic_
    text = re.sub(r"_{1,2}(.*?)_{1,2}", r"\1", text)
    # Convert inline backticks (`code`)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    # Strip LLM severity directives (matplotlib cannot render them)
    # Canonical :span[text]{.text-color} or :sample[text]{.text-color} → text
    text = re.sub(r":(?:span|sample)\[([^\]]*)\]\s*\{[^}]*\}", r"\1", text)
    # Parenthesized variant :span[text] (.text-color) → text
    text = re.sub(r":(?:span|sample)\[([^\]]*)\]\s*\(\.text-\w+\)", r"\1", text)
    # Dangling :span[text] without any directive → text
    text = re.sub(r":(?:span|sample)\[([^\]]*)\]", r"\1", text)
    # Bare {.text-color} at end of text
    text = re.sub(r"\s*\{\s*\.text-(?:red|orange|yellow|green)\s*\}", "", text)
    # Bold-wrapped {.text-color}
    text = re.sub(r"\s*\*{1,3}\{\s*\.text-(?:red|orange|yellow|green)\s*\}\*{1,3}", "", text)
    # Bare .text-color (no delimiters)
    text = re.sub(r"\s+\.text-(?:red|orange|yellow|green)\b", "", text)
    # (.text-color) parenthesized standalone
    text = re.sub(r"\s*\(\.text-(?:red|orange|yellow|green)\)", "", text)
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
    """Merge original QC PDF with AI commentary slides: all original pages first,
    then all AI slides appended in section order. Previous attempts to interleave
    at exact page boundaries failed because R-generated PDF pages rarely contain
    searchable text for section matching. Appending in order is reliable and each
    AI slide's title already identifies which section it belongs to."""
    # Try pypdf first
    try:
        import pypdf
        reader_orig = pypdf.PdfReader(orig_pdf)
        reader_ai = pypdf.PdfReader(ai_pdf_pages_file)
        writer = pypdf.PdfWriter()
        for page in reader_orig.pages:
            writer.add_page(page)
        for page in reader_ai.pages:
            writer.add_page(page)
        with open(output_pdf, "wb") as fh:
            writer.write(fh)
        llm_common.log(f"[QC PDF AI] Merged {len(reader_orig.pages)} QC pages + "
                       f"{len(reader_ai.pages)} AI slides into: {output_pdf}")
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] pypdf merge error: {e}")

    # Fallback: fitz (PyMuPDF)
    try:
        import fitz
        doc_orig = fitz.open(orig_pdf)
        doc_ai = fitz.open(ai_pdf_pages_file)
        doc_out = fitz.open()
        doc_out.insert_pdf(doc_orig)
        doc_out.insert_pdf(doc_ai)
        doc_out.save(output_pdf)
        llm_common.log(f"[QC PDF AI] Merged {len(doc_orig)} QC pages + "
                       f"{len(doc_ai)} AI slides (fitz) into: {output_pdf}")
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] fitz merge error: {e}")


    return False

def _match_section_to_tsv(section_pdf_name, tsv_files):
    """Match a section PDF filename (e.g. '03_library_size.pdf') to a TSV file
    (e.g. '03_library_sizes.tsv') using the two-digit numeric prefix."""
    prefix = section_pdf_name[:2]  # e.g. "03"
    for tsv in tsv_files:
        if tsv[:2] == prefix:
            return tsv
    return None

def _concat_pdfs(pdf_list, output_path):
    """Concatenate a list of PDF file paths into a single output PDF."""
    try:
        import pypdf
        writer = pypdf.PdfWriter()
        for pdf_path in pdf_list:
            reader = pypdf.PdfReader(pdf_path)
            for page in reader.pages:
                writer.add_page(page)
        with open(output_path, "wb") as fh:
            writer.write(fh)
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] pypdf concat error: {e}")
    try:
        import fitz
        doc_out = fitz.open()
        for pdf_path in pdf_list:
            doc_in = fitz.open(pdf_path)
            doc_out.insert_pdf(doc_in)
        doc_out.save(output_path)
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] fitz concat error: {e}")
    return False

def _append_ai_to_section(section_pdf, ai_slide_pdf, ai_slide_idx, output_pdf):
    """Append a single AI slide (page ai_slide_idx from ai_slide_pdf) to section_pdf."""
    try:
        import pypdf
        writer = pypdf.PdfWriter()
        reader_sec = pypdf.PdfReader(section_pdf)
        reader_ai = pypdf.PdfReader(ai_slide_pdf)
        for page in reader_sec.pages:
            writer.add_page(page)
        if ai_slide_idx < len(reader_ai.pages):
            writer.add_page(reader_ai.pages[ai_slide_idx])
        with open(output_pdf, "wb") as fh:
            writer.write(fh)
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] pypdf append error: {e}")
    try:
        import fitz
        doc_sec = fitz.open(section_pdf)
        doc_ai = fitz.open(ai_slide_pdf)
        doc_sec.insert_pdf(doc_ai, from_page=ai_slide_idx, to_page=ai_slide_idx)
        doc_sec.save(output_pdf)
        return True
    except Exception as e:
        llm_common.log(f"[QC PDF AI] fitz append error: {e}")
    return False

def main():
    p = argparse.ArgumentParser(description="QC Figures AI PDF generator.")
    p.add_argument("--tables-dir", required=True, help="Directory containing exported TSV tables.")
    p.add_argument("--pdf", default="", help="Original QC PDF file to pair with (optional).")
    p.add_argument("--sections-dir", default="", help="Directory with per-section PDFs from R (enables per-section interleaving).")
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

    # --- Generate AI commentary slides ---
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

    # --- Per-section interleaving mode ---
    sections_dir = args.sections_dir
    if sections_dir and os.path.isdir(sections_dir):
        section_pdfs = sorted([f for f in os.listdir(sections_dir) if f.endswith(".pdf")])
        if not section_pdfs:
            llm_common.log(f"[QC PDF AI] No section PDFs found in {sections_dir}")
            return

        llm_common.log(f"[QC PDF AI] Per-section interleaving: {len(section_pdfs)} section PDFs, "
                       f"{len(tsv_to_slide_idx)} AI slides")

        # Build commented section PDFs (original + AI slide appended)
        commented_section_paths = []
        original_section_paths = []
        for sec_pdf_name in section_pdfs:
            sec_pdf_path = os.path.join(sections_dir, sec_pdf_name)
            original_section_paths.append(sec_pdf_path)

            # Find matching TSV → AI slide index
            matched_tsv = _match_section_to_tsv(sec_pdf_name, tsv_files)
            slide_idx = tsv_to_slide_idx.get(matched_tsv) if matched_tsv else None

            if slide_idx is not None:
                # Append AI slide to this section
                commented_path = sec_pdf_path.replace(".pdf", "_commented.pdf")
                if _append_ai_to_section(sec_pdf_path, out_pdf_standalone, slide_idx, commented_path):
                    commented_section_paths.append(commented_path)
                    llm_common.log(f"[QC PDF AI] Appended AI slide for {matched_tsv} to {sec_pdf_name}")
                else:
                    commented_section_paths.append(sec_pdf_path)  # fallback: use original
            else:
                commented_section_paths.append(sec_pdf_path)  # no AI slide for this section

        # Assemble monolithic *_QC.pdf from original sections
        if args.pdf:
            qc_pdf_path = args.pdf
        else:
            qc_pdf_path = os.path.join(parent_dir, "assembled_QC.pdf")
        if _concat_pdfs(original_section_paths, qc_pdf_path):
            llm_common.log(f"[QC PDF AI] Assembled monolithic QC PDF: {qc_pdf_path}")

        # Assemble *_QC_commented.pdf from commented sections
        if args.pdf:
            commented_pdf = args.pdf.replace(".pdf", "_commented.pdf")
        else:
            commented_pdf = os.path.join(parent_dir, "assembled_QC_commented.pdf")
        if _concat_pdfs(commented_section_paths, commented_pdf):
            llm_common.scrub_file(commented_pdf, secrets)
            llm_common.log(f"[QC PDF AI] Assembled commented QC PDF: {commented_pdf}")

        # Clean up sections directory
        import shutil
        try:
            shutil.rmtree(sections_dir)
            llm_common.log(f"[QC PDF AI] Cleaned up sections directory: {sections_dir}")
        except OSError as e:
            llm_common.log(f"[QC PDF AI] Could not remove sections directory: {e}")

    elif args.pdf and os.path.isfile(args.pdf):
        # Fallback: append all AI slides at end of monolithic PDF
        commented_pdf = args.pdf.replace(".pdf", "_commented.pdf")
        if merge_pdfs_if_possible(args.pdf, out_pdf_standalone, commented_pdf, tsv_to_slide_idx=tsv_to_slide_idx):
            llm_common.scrub_file(commented_pdf, secrets)

if __name__ == "__main__":
    main()

