#!/usr/bin/env python3
"""Stop a small RNA run when fastp kept almost nothing: the adapter was not found, or the reads are longer than small_rna_length."""
import glob
import json
import os
import sys

fastp_dir, max_len = sys.argv[1], sys.argv[2]
bad = []
for f in sorted(glob.glob(os.path.join(fastp_dir, "*_report.json"))):
    try:
        j = json.load(open(f))
    except Exception:
        continue
    before = j["summary"]["before_filtering"]
    after = j["summary"]["after_filtering"]
    fr = j.get("filtering_result", {})
    adapter = j.get("adapter_cutting", {}).get("read1_adapter_sequence", "none detected")
    b, a = before.get("total_reads", 0), after.get("total_reads", 0)
    if b and (a == 0 or fr.get("too_long_reads", 0) > 0.5 * b):
        bad.append(f"  {os.path.basename(f).replace('_report.json', '')}: {b} reads, {a} kept, {fr.get('too_long_reads', 0)} longer than {max_len} nt after trimming, "
                   f"{fr.get('too_short_reads', 0)} too short; mean length before trimming {before.get('read1_mean_length')} nt; adapter: {adapter}")
if bad:
    sys.stderr.write("ERROR: fastp kept almost no reads in small RNA mode:\n" + "\n".join(bad) + "\n"
                     "Most reads stayed longer than the maximum of small_rna_length, so the 3' adapter was not trimmed: set small_rna_adapter to the adapter of the kit "
                     "(e.g. TGGAATTCTCGGGTGCCAAGG for Illumina TruSeq Small RNA, AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC for NEBNext, AACTGTAGGCACCATCAAT for QIAseq), "
                     "or raise the maximum of small_rna_length if the reads were uploaded already trimmed to a fixed length.\n")
    sys.exit(1)
