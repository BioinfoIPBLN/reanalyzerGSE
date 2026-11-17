#!/usr/bin/env python3
"""
llm_simplify_metadata.py - AI-assisted sample and condition name streamlining for reanalyzerGSE.

Given sample metadata files (samples_info.txt, phenodata_extracted.txt, design_possible_full*.txt)
in a reads_study_info directory, this script calls the OpenAI-compatible LLM endpoint to generate
short, clean, discriminative sample names and condition labels.

Rules:
- Condition names must be concise (e.g. "NTG" or "MCK_dOTC" instead of 100-char descriptions).
- Sample names must be clean, discriminative, and preserve GSM/replicate tags.
- Output contains only alphanumeric characters and underscores [A-Za-z0-9_].
- Each design is relabelled level by level, so the grouping of the samples never changes.
- Writes streamlined samples_info.txt, sample_names.txt, and design files.
"""

import argparse
import os
import shutil
import sys
import re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common

PROMPT_TEMPLATE = """You are a bioinformatics metadata assistant. Always respond in English.
Your task is to take raw, overly verbose GEO sample metadata and simplify:
1. The sample names.
2. The labels of the groups in each grouping of the samples.

SAMPLES (ID, current sample name, current condition):
{samples}

GROUPINGS of the samples. Each grouping is a separate design and is relabelled on its own.
Each line: level ID, current label, number of samples with that label.
{groupings}

Additional phenotype context (if helpful):
{phenodata}

RULES:
- Group labels must be short, clean, and discriminative (e.g., "NTG" vs "MCK_dOTC" or "WT" vs "KO").
- Within one grouping, different levels MUST get different labels. Never merge or split levels.
- Use the same word for the same concept across groupings (e.g., "KO" in a genotype grouping and "KO_OGD" in a combined grouping).
- Sample names must be concise, combining the condition, replicate tag, and GSM accession if present, and must all be different.
- Replicate tags MUST always use the exact capitalized standard format "_Rep1", "_Rep2", "_Rep3", etc. (e.g. "NTG_Rep1_GSM5770284"). Never use "rep1", "rep_1", "replicate1", or "biological_rep1".
- Strictly use ONLY alphanumeric characters and underscores [A-Za-z0-9_]. No spaces, hyphens, parentheses, or special characters.
- Keep total sample name length under 60 characters.

OUTPUT: one line per sample ID and one line per level ID, exactly two tab-separated fields, nothing else:
<ID>\\t<NEW_NAME_OR_LABEL>
For example:
S1\\tKO_Rep1_GSM6285636
D1.L1\\tKO
"""

DESIGN_RE = re.compile(r"^design_possible_full_?([0-9]+).*\.txt$")


def parse_args():
    parser = argparse.ArgumentParser(description="LLM-based sample and condition name simplification")
    parser.add_argument("--info-dir", required=True, help="Path to reads_study_info directory")
    parser.add_argument("--endpoint", help="LLM endpoint URL")
    parser.add_argument("--model", help="LLM model name")
    parser.add_argument("--api-key", help="API key")
    parser.add_argument("--timeout", type=int, default=120, help="LLM timeout in seconds")
    return parser.parse_args()


def normalize_replicates(name):
    name = re.sub(r'(?i)(?:biological_)?(?:bio_)?(?:replicate|rep)[_.-]?([0-9]+)', r'_Rep\1', name)
    return name


def clean(text):
    text = re.sub(r'[^A-Za-z0-9_]', '_', text.strip())
    return re.sub(r'_+', '_', text).strip('_')


def design_files(info_dir):
    found = []
    for fname in os.listdir(info_dir):
        m = DESIGN_RE.match(fname)
        if m:
            found.append((int(m.group(1)), fname))
    return [f for _, f in sorted(found)]


def read_lines(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read().splitlines()


def write_lines(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def levels_of(values):
    return list(dict.fromkeys(v for v in values if v.strip()))


def partition_key(values):
    first = {}
    return tuple(first.setdefault(v, len(first)) for v in values)


def parse_response(text):
    answers = {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        parts = [p.strip() for p in line.split("\t")]
        key = parts[0].strip("`*-• ")
        if key and parts[-1]:
            answers[key] = parts[-1]
    return answers


def main():
    args = parse_args()
    info_dir = args.info_dir
    samples_info_path = os.path.join(info_dir, "samples_info.txt")
    phenodata_path = os.path.join(info_dir, "phenodata_extracted.txt")

    if not os.path.exists(samples_info_path):
        print(f"llm_simplify_metadata.py: {samples_info_path} not found. Skipping.")
        sys.exit(0)

    raw_lines = [l for l in read_lines(samples_info_path) if l.strip()]
    if not raw_lines:
        sys.exit(0)
    rows = [l.split("\t") for l in raw_lines]
    has_cond = all(len(r) >= 3 for r in rows)

    phenodata_text = ""
    if os.path.exists(phenodata_path):
        with open(phenodata_path, "r", encoding="utf-8", errors="replace") as f:
            phenodata_text = f.read().strip()[:4000]

    groupings = [{"id": f"D{k}", "title": fname, "file": fname, "values": read_lines(os.path.join(info_dir, fname))}
                 for k, fname in enumerate(design_files(info_dir), start=1)]
    cond_values = [r[2] for r in rows] if has_cond else []
    cond_link = None
    cond_group = None
    if has_cond:
        cond_key = partition_key(cond_values)
        cond_link = next((g for g in groupings if len(g["values"]) == len(cond_values)
                          and partition_key(g["values"]) == cond_key), None)
        if cond_link is None:
            cond_group = {"id": "C", "title": "condition column of samples_info", "file": None, "values": cond_values}
            groupings.insert(0, cond_group)
    for g in groupings:
        g["levels"] = levels_of(g["values"])

    samples_block = "\n".join(
        f"S{i}\t{r[1] if len(r) > 1 else r[0]}\t{r[2] if has_cond else '-'}" for i, r in enumerate(rows, start=1))
    groupings_block = "\n".join(
        f"{g['id']} ({g['title']}):\n" + "\n".join(
            f"{g['id']}.L{j}\t{lev}\t{sum(v == lev for v in g['values'])}" for j, lev in enumerate(g["levels"], start=1))
        for g in groupings) or "None"

    prompt = PROMPT_TEMPLATE.format(
        samples=samples_block,
        groupings=groupings_block,
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

    answers = parse_response(response_text)

    new_names = [clean(normalize_replicates(answers.get(f"S{i}", ""))) for i in range(1, len(rows) + 1)]
    names_ok = all(new_names) and len(set(new_names)) == len(new_names)
    if not names_ok:
        print("llm_simplify_metadata.py: the AI sample names are missing, empty or not unique. Keeping the original sample names.")
        new_names = [r[1] if len(r) > 1 else r[0] for r in rows]

    mapped = 0
    for g in groupings:
        labels = [clean(answers.get(f"{g['id']}.L{j}", "")) for j in range(1, len(g["levels"]) + 1)]
        if g["levels"] and all(labels) and len(set(labels)) == len(labels):
            g["map"] = dict(zip(g["levels"], labels))
            mapped += 1
        else:
            g["map"] = None
            print(f"llm_simplify_metadata.py: the AI labels for {g['title']} are missing, empty or merge different groups. Keeping its original labels.")

    if not names_ok and mapped == 0:
        print("llm_simplify_metadata.py: nothing usable in the AI answer. Keeping original metadata.")
        sys.exit(0)

    def relabel(g, values):
        if g["map"] is None:
            return values
        return [g["map"].get(v, v) if v.strip() else v for v in values]

    if cond_link is not None:
        new_conditions = relabel(cond_link, cond_link["values"])
    elif cond_group is not None:
        new_conditions = relabel(cond_group, cond_values)
    else:
        new_conditions = []

    shutil.copy2(samples_info_path, samples_info_path + ".raw_bak")
    new_samples_info = []
    for i, r in enumerate(rows):
        fields = [r[0], new_names[i]] + ([new_conditions[i]] if has_cond else []) + r[3:]
        new_samples_info.append("\t".join(fields))
    write_lines(samples_info_path, new_samples_info)
    write_lines(os.path.join(info_dir, "sample_names.txt"), new_names)

    for g in groupings:
        if g["file"] is None or g["map"] is None:
            continue
        new_values = relabel(g, g["values"])
        write_lines(os.path.join(info_dir, g["file"]), new_values)
        companion = os.path.join(info_dir, "design_possible" + g["file"][len("design_possible_full"):])
        if os.path.exists(companion):
            write_lines(companion, levels_of(new_values))

    if has_cond and not any(g["file"] for g in groupings):
        write_lines(os.path.join(info_dir, "design_possible_full_1.txt"), new_conditions)
        write_lines(os.path.join(info_dir, "design_possible_1.txt"), levels_of(new_conditions))

    print(f"\n==========================================================================")
    print(f"  AI Metadata Streamlining Complete ({len(rows)} samples)")
    print(f"==========================================================================")
    for g in groupings:
        if g["map"] is not None:
            print(f"{g['title']}: " + ", ".join(f"{old} -> {new}" for old, new in g["map"].items()))
    print(f"\nStreamlined Sample Mapping:")
    print(f"{'SRR ID':<15} {'Simplified Sample Name':<45} {'Condition'}")
    print(f"-" * 75)
    for line in new_samples_info:
        parts = line.split('\t')
        print(f"{parts[0]:<15} {parts[1]:<45} {parts[2] if len(parts) > 2 else ''}")
    print(f"==========================================================================\n")

if __name__ == "__main__":
    main()
