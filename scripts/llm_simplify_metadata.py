#!/usr/bin/env python3
"""
llm_simplify_metadata.py - AI-assisted sample and condition name streamlining for reanalyzerGSE.

Given sample metadata files (samples_info.txt, phenodata_extracted.txt, etc.) in a reads_study_info
directory, this script calls the OpenAI-compatible LLM endpoint to generate short, clean,
discriminative sample names and condition labels.

Rules:
- Condition names must be concise (e.g. "NTG" or "MCK_dOTC" instead of 100-char descriptions).
- Sample names must be clean, discriminative, and preserve GSM/replicate tags.
- Output contains only alphanumeric characters and underscores [A-Za-z0-9_].
- Writes streamlined samples_info.txt, sample_names.txt, and design files.
"""

import argparse
import os
import sys
import re

# Import shared llm_common helper
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common

PROMPT_TEMPLATE = """You are a bioinformatics metadata assistant. Always respond in English.
Your task is to take raw, overly verbose GEO sample metadata and simplify both:
1. The condition/group names (column 3 of samples_info).
2. The descriptive sample names (column 2 of samples_info).

RULES:
- Condition names must be short, clean, and discriminative (e.g., "NTG" vs "MCK_dOTC" or "WT" vs "KO").
- Sample names must be concise, combining the condition, replicate tag, and GSM accession if present.
- Replicate tags MUST always use the exact capitalized standard format "_Rep1", "_Rep2", "_Rep3", etc. (e.g. "NTG_Rep1_GSM5770284", "MCK_dOTC_Rep2_GSM5770287"). Never use "rep1", "rep_1", "replicate1", or "biological_rep1".
- Strictly use ONLY alphanumeric characters and underscores [A-Za-z0-9_]. No spaces, hyphens, parentheses, or special characters.
- Keep total sample name length under 60 characters.
- Output MUST have the exact same number of lines as input.
- Output format per line (tab-separated):
  <SRR_ID_OR_ACCESSION>\t<SIMPLIFIED_SAMPLE_NAME>\t<SIMPLIFIED_CONDITION>

Here is the raw samples_info file:
{samples_info}

Here is additional phenotype context (if helpful):
{phenodata}
"""

def parse_args():
    parser = argparse.ArgumentParser(description="LLM-based sample and condition name simplification")
    parser.add_argument("--info-dir", required=True, help="Path to reads_study_info directory")
    parser.add_argument("--endpoint", help="LLM endpoint URL")
    parser.add_argument("--model", help="LLM model name")
    parser.add_argument("--api-key", help="API key")
    parser.add_argument("--timeout", type=int, default=120, help="LLM timeout in seconds")
    return parser.parse_args()

def normalize_replicates(name):
    # Convert variations like rep1, rep_1, replicate1, bio_rep1, biological_rep1 to _Rep1
    name = re.sub(r'(?i)(?:biological_)?(?:bio_)?(?:replicate|rep)[_.-]?([0-9]+)', r'_Rep\1', name)
    return name

def main():
    args = parse_args()
    info_dir = args.info_dir
    samples_info_path = os.path.join(info_dir, "samples_info.txt")
    phenodata_path = os.path.join(info_dir, "phenodata_extracted.txt")

    if not os.path.exists(samples_info_path):
        print(f"llm_simplify_metadata.py: {samples_info_path} not found. Skipping.")
        sys.exit(0)

    # Read raw samples_info
    with open(samples_info_path, "r", encoding="utf-8", errors="replace") as f:
        samples_info_text = f.read().strip()

    if not samples_info_text:
        sys.exit(0)

    # Read phenodata if available
    phenodata_text = ""
    if os.path.exists(phenodata_path):
        with open(phenodata_path, "r", encoding="utf-8", errors="replace") as f:
            phenodata_text = f.read().strip()[:4000] # cap context

    prompt = PROMPT_TEMPLATE.format(
        samples_info=samples_info_text,
        phenodata=phenodata_text if phenodata_text else "None"
    )

    print("llm_simplify_metadata.py: Requesting AI metadata simplification...")
    
    endpoint = args.endpoint or os.environ.get("LLM_ENDPOINT")
    model = args.model or os.environ.get("LLM_MODEL")
    api_key = args.api_key or os.environ.get("LLM_API_KEY", "dummy")

    if not endpoint or not model:
        print("llm_simplify_metadata.py: No LLM endpoint or model set. Keeping original metadata.")
        sys.exit(0)

    try:
        messages = [
            {"role": "system", "content": "You are a careful bioinformatics assistant. Always respond in English."},
            {"role": "user", "content": prompt}
        ]
        response_text, usage = llm_common.chat_completion(
            endpoint, model, api_key, messages, timeout=args.timeout
        )
    except Exception as e:
        print(f"llm_simplify_metadata.py: AI query error ({e}). Keeping original metadata.")
        sys.exit(0)

    if not response_text:
        print("llm_simplify_metadata.py: AI simplification returned no response. Keeping original metadata.")
        sys.exit(0)

    # Clean lines from response
    lines = [l.strip() for l in response_text.strip().splitlines() if l.strip() and "\t" in l]

    raw_lines = [l for l in samples_info_text.splitlines() if l.strip()]

    if len(lines) != len(raw_lines):
        print(f"llm_simplify_metadata.py: Response line count ({len(lines)}) mismatch with raw input ({len(raw_lines)}). Keeping original metadata.")
        sys.exit(0)

    # Validate and sanitize response
    new_samples_info = []
    new_sample_names = []
    new_conditions = []

    for i, line in enumerate(lines):
        parts = line.split("\t")
        if len(parts) < 3:
            print(f"llm_simplify_metadata.py: Invalid line output: '{line}'. Keeping original metadata.")
            sys.exit(0)

        srr_id = raw_lines[i].split("\t")[0]
        s_name = parts[1].strip()
        s_cond = re.sub(r'[^A-Za-z0-9_]', '_', parts[2].strip())

        # Normalize replicate patterns to _RepN
        s_name = normalize_replicates(s_name)

        s_name = re.sub(r'[^A-Za-z0-9_]', '_', s_name)

        # Clean consecutive underscores
        s_name = re.sub(r'_+', '_', s_name).strip('_')
        s_cond = re.sub(r'_+', '_', s_cond).strip('_')

        new_samples_info.append(f"{srr_id}\t{s_name}\t{s_cond}")
        new_sample_names.append(s_name)
        new_conditions.append(s_cond)

    # Backup original files
    os.system(f"cp '{samples_info_path}' '{samples_info_path}.raw_bak'")

    # Write streamlined samples_info.txt
    with open(samples_info_path, "w", encoding="utf-8") as f:
        f.write("\n".join(new_samples_info) + "\n")

    # Write streamlined sample_names.txt
    sample_names_path = os.path.join(info_dir, "sample_names.txt")
    with open(sample_names_path, "w", encoding="utf-8") as f:
        f.write("\n".join(new_sample_names) + "\n")

    # Overwrite all design_possible_full_*.txt and design_possible_*.txt files in info_dir
    unique_conds = list(dict.fromkeys(new_conditions))
    design_files_found = 0
    for fname in os.listdir(info_dir):
        if fname.startswith("design_possible_full_"):
            with open(os.path.join(info_dir, fname), "w", encoding="utf-8") as fh:
                fh.write("\n".join(new_conditions) + "\n")
            design_files_found += 1
        elif fname.startswith("design_possible_") and not fname.startswith("design_possible_full_"):
            with open(os.path.join(info_dir, fname), "w", encoding="utf-8") as fh:
                fh.write("\n".join(unique_conds) + "\n")
            design_files_found += 1

    # Ensure fallback design_possible_full_1.txt and design_possible_1.txt are written
    design_full_path = os.path.join(info_dir, "design_possible_full_1.txt")
    design_uniq_path = os.path.join(info_dir, "design_possible_1.txt")
    with open(design_full_path, "w", encoding="utf-8") as f:
        f.write("\n".join(new_conditions) + "\n")
    with open(design_uniq_path, "w", encoding="utf-8") as f:
        f.write("\n".join(unique_conds) + "\n")

    print(f"\n==========================================================================")
    print(f"  AI Metadata Streamlining Complete ({len(lines)} samples)")
    print(f"==========================================================================")
    print(f"Simplified Conditions/Groups: {', '.join(unique_conds)}\n")
    print(f"Streamlined Sample Mapping:")
    print(f"{'SRR ID':<15} {'Simplified Sample Name':<45} {'Condition'}")
    print(f"-" * 75)
    for line in new_samples_info:
        parts = line.split('\t')
        print(f"{parts[0]:<15} {parts[1]:<45} {parts[2]}")
    print(f"==========================================================================\n")

if __name__ == "__main__":
    main()
