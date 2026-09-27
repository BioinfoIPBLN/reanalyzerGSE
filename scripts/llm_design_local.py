#!/usr/bin/env python3
"""
llm_design_local.py - LLM-proposed experimental design for local reads in reanalyzerGSE.

Reads the ordered sample names (one per line) from --samples or stdin, asks the
configured OpenAI-compatible LLM (LLM_ENDPOINT / LLM_MODEL / LLM_API_KEY) for one
condition per sample, and prints them as a comma-separated list on stdout.
Diagnostics go to stderr.

Exit codes: 0 design printed, 2 no LLM configured, 3 LLM request failed,
4 no usable answer.
"""

import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import llm_common

PROMPT_TEMPLATE = """You are assigning the experimental condition of each sample in an RNA-seq study{study}.
Below are the {n} sample names, one per line, in the order the pipeline will use. They come from the raw FASTQ file names, with the read suffixes (_1/_2, _R1/_R2) already removed.

Infer the biological condition (e.g. treatment, genotype, dose, time point, tissue) that each sample belongs to. Replicates of the same condition must get exactly the same label.

RULES:
- Return exactly {n} labels, one per sample, in the same order as the input.
- Drop replicate identifiers (e.g. Rep1, rep_2, replicate3) and sequencing-run artefacts (e.g. lane L001, sample-sheet index S12); keep everything that distinguishes one condition from another.
- Reuse the wording already present in the names rather than inventing new terms.
- Use ONLY letters, digits and underscores, start each label with a letter, and keep each under 40 characters.
- If the names give no hint of a condition, use the same label for every sample.

Respond with JSON only, in exactly this form and nothing else:
{{"conditions": ["label_for_sample_1", "label_for_sample_2", ...]}}

Sample names:
{samples}
"""


def parse_args():
    parser = argparse.ArgumentParser(description="LLM-proposed experimental design for local reads")
    parser.add_argument("--samples", help="File with one sample name per line (default: stdin)")
    parser.add_argument("--study", default="", help="Study name, given to the LLM as context")
    parser.add_argument("--timeout", type=int, default=300, help="LLM timeout in seconds")
    return parser.parse_args()


def extract_conditions(text):
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.S).strip()
    candidates = [text]
    for opener, closer in (("{", "}"), ("[", "]")):
        start, end = text.find(opener), text.rfind(closer)
        if start != -1 and end > start:
            candidates.append(text[start:end + 1])
    for candidate in candidates:
        try:
            data = json.loads(candidate)
        except ValueError:
            continue
        if isinstance(data, dict):
            data = data.get("conditions")
        if isinstance(data, list) and all(isinstance(x, (str, int, float)) for x in data):
            return [str(x) for x in data]
    return None


def clean_label(label):
    label = re.sub(r"[^A-Za-z0-9_]", "_", label.strip())
    label = re.sub(r"_+", "_", label).strip("_")
    if label and not label[0].isalpha():
        label = "G_" + label
    return label


def main():
    args = parse_args()
    if args.samples:
        with open(args.samples, encoding="utf-8", errors="replace") as fh:
            samples = [line.strip() for line in fh if line.strip()]
    else:
        samples = [line.strip() for line in sys.stdin if line.strip()]
    if not samples:
        llm_common.log("llm_design_local.py: no sample names were given")
        sys.exit(4)

    endpoint = os.environ.get("LLM_ENDPOINT")
    model = os.environ.get("LLM_MODEL")
    api_key = os.environ.get("LLM_API_KEY", "dummy")
    if not endpoint or not model:
        llm_common.log("llm_design_local.py: no LLM endpoint or model is configured")
        sys.exit(2)
    secrets = llm_common.secret_values(endpoint=endpoint, api_key=api_key)

    n = len(samples)
    messages = [
        {"role": "system", "content": "You are a careful bioinformatics assistant. Always respond in English."},
        {"role": "user", "content": PROMPT_TEMPLATE.format(
            n=n, study=f" named '{args.study}'" if args.study else "", samples="\n".join(samples))},
    ]
    conditions = None
    for attempt in (1, 2):
        try:
            text, _usage = llm_common.chat_completion(endpoint, model, api_key, messages, timeout=args.timeout)
        except Exception as e:
            llm_common.log("llm_design_local.py: the LLM request failed: " + llm_common.mask(str(e), secrets)[0])
            sys.exit(3)
        text = text or ""
        conditions = extract_conditions(text)
        if conditions is None:
            problem = "the answer contained no JSON list of conditions"
        elif len(conditions) != n:
            problem = f"it gave {len(conditions)} conditions for {n} samples"
        else:
            conditions = [clean_label(c) for c in conditions]
            problem = None if all(conditions) else "some conditions were empty"
        if problem is None:
            break
        llm_common.log(f"llm_design_local.py: attempt {attempt} not usable: {problem}")
        if attempt == 2:
            sys.exit(4)
        messages += [
            {"role": "assistant", "content": text},
            {"role": "user", "content": f"That answer is not usable: {problem}. Reply again with JSON only, "
                                        f"with exactly {n} labels in the input order."},
        ]

    if n > 1 and len(set(conditions)) == n:
        llm_common.log("llm_design_local.py: WARNING: every sample got its own condition, so no condition has replicates")
    print(",".join(conditions))


if __name__ == "__main__":
    main()
