#!/usr/bin/env python3
"""
validate_alignment_counts.py

Universal validation gate for reanalyzerGSE across STAR, HISAT2, Kallisto and bowtie1 (small RNA mode).
Verifies that 100% of expected samples produced valid BAMs (or pseudoalignment tables),
valid count matrices, and zero failed jobs in GNU parallel.

If any sample failed, extracts relevant log excerpts and exits with code 1,
preventing the pipeline from cascading into Step 4 with missing or corrupted data.
"""

import argparse
import glob
import json
import os
import sys


def read_samples_info(path):
    """Return the tab-separated fields of every non-empty row of a samples_info.txt file."""
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.strip():
                rows.append(line.split("\t"))
    return rows


def read_stem(name):
    """Sample name of a FASTQ file name, or None when it is not a FASTQ file."""
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            s = name[:-len(suffix)]
            for end in ["_1", "_2", "_R1", "_R2"]:
                if s.endswith(end):
                    s = s[:-len(end)]
            return s or None
    return None


def sample_names_on_disk(analysis_dir):
    """Sample names the raw reads and the aligner outputs are named after."""
    names = set()
    raw_reads_dir = os.path.join(analysis_dir, "raw_reads")
    if os.path.isdir(raw_reads_dir):
        names.update(s for s in map(read_stem, os.listdir(raw_reads_dir)) if s)
    for subdir, suffix in (("star_results", "_STAR.bam"), ("hisat2_results", "_hisat2.bam"), ("bowtie1_results", "_bowtie1.bam")):
        for p in glob.glob(os.path.join(analysis_dir, "miARma_out*", subdir, "*" + suffix)):
            names.add(os.path.basename(p)[:-len(suffix)])
    for p in glob.glob(os.path.join(analysis_dir, "miARma_out*", "kallisto_results", "*")):
        b = os.path.basename(p)
        names.add(b[:-len("_kallisto")] if b.endswith("_kallisto") else b)
    return names


def sample_column(rows, names):
    """Index of the samples_info.txt column the files are named after: the first column for
    local reads, the second (sample name, after the SRR accession) for GEO/SRA downloads."""
    best, best_hits = 0, -1
    for col in (0, 1):
        hits = sum(1 for r in rows if len(r) > col and r[col].strip() in names)
        if hits > best_hits:
            best, best_hits = col, hits
    return best


def parse_expected_samples(analysis_dir, samples_info_path=None):
    """Extract expected sample names from samples_info.txt or raw_reads directory.

    Returns (samples, rows_without_sample)."""
    default_info = os.path.join(analysis_dir, "reads_study_info", "samples_info.txt")
    for path in (samples_info_path, default_info):
        if path and os.path.exists(path):
            rows = read_samples_info(path)
            if rows:
                col = sample_column(rows, sample_names_on_disk(analysis_dir))
                samples = []
                rows_without_sample = 0
                for r in rows:
                    sample = r[col].strip() if len(r) > col else ""
                    if not sample:
                        rows_without_sample += 1
                    elif sample not in samples:
                        samples.append(sample)
                return samples, rows_without_sample

    expected = []
    raw_reads_dir = os.path.join(analysis_dir, "raw_reads")
    if os.path.exists(raw_reads_dir):
        for f in os.listdir(raw_reads_dir):
            s = read_stem(f)
            if s and s not in expected:
                expected.append(s)

    return sorted(expected), 0


def find_star_error(log_file):
    """Extract the primary fatal error line from a STAR log file."""
    if not os.path.exists(log_file):
        return ""
    try:
        with open(log_file, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if any(term in line for term in ["FATAL ERROR", "wrong read ID line format", "EXITING because of FATAL ERROR", "Segmentation fault"]):
                excerpt = line.strip()
                if i + 1 < len(lines) and lines[i + 1].strip().startswith("Offending line"):
                    excerpt += " | " + lines[i + 1].strip()
                return excerpt
    except Exception:
        pass
    return ""


def find_hisat2_error(log_file):
    """Extract error details from a HISAT2 log file."""
    if not os.path.exists(log_file):
        return ""
    try:
        with open(log_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if any(term in line.lower() for term in ["error", "segmentation fault", "failed", "cannot open"]):
                    return line.strip()
    except Exception:
        pass
    return ""


def find_bowtie1_error(log_file):
    """Extract error details from a bowtie1 stderr log (its normal summary lines mention reads that 'failed to align')."""
    if not os.path.exists(log_file):
        return ""
    try:
        with open(log_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                low = line.lower()
                if low.startswith("# reads") or low.startswith("reported "):
                    continue
                if any(term in low for term in ["error", "segmentation fault", "cannot open", "no such file", "could not locate"]):
                    return line.strip()
    except Exception:
        pass
    return ""


def find_kallisto_error(log_file):
    """Extract error details from a Kallisto log file."""
    if not os.path.exists(log_file):
        return ""
    try:
        with open(log_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if any(term in line.lower() for term in ["error", "zero targets", "failed", "exception"]):
                    return line.strip()
    except Exception:
        pass
    return ""


def check_parallel_logs(analysis_dir):
    """Scan all *_log_parallel.txt files for failed jobs."""
    failed_parallel = {}
    parallel_logs = glob.glob(os.path.join(analysis_dir, "miARma_out*", "*_log_parallel.txt"))
    for plog in parallel_logs:
        log_type = os.path.basename(plog)
        try:
            with open(plog, "r", encoding="utf-8", errors="replace") as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 7:
                        exitval = parts[6].strip()
                        cmd = parts[8] if len(parts) > 8 else ""
                        if exitval.isdigit() and int(exitval) != 0:
                            # Extract sample name if possible
                            sample_match = None
                            for word in cmd.split():
                                if "_STAR" in word or "_hisat2" in word or "_readcount" in word or "abundance" in word:
                                    s = os.path.basename(word).split("_")[0]
                                    sample_match = s
                                    break
                            failed_parallel[f"{log_type} (Exitval {exitval})"] = cmd[:100]
        except Exception:
            pass
    return failed_parallel


def has_alignment_output(analysis_dir, aligner):
    """True when the aligner produced something. Used to tell a run that has
    nothing to validate yet from one whose sample list went missing."""
    patterns = {
        "star":     [("star_results", "*_STAR.bam")],
        "hisat2":   [("hisat2_results", "*_hisat2.bam")],
        "kallisto": [("kallisto_results", "*")],
        "bowtie1":  [("bowtie1_results", "*_bowtie1.bam")],
    }.get(aligner, [("star_results", "*_STAR.bam")])
    for subdir, leaf in patterns:
        if glob.glob(os.path.join(analysis_dir, "miARma_out*", subdir, leaf)):
            return True
    return False


def validate_alignment_and_counts(analysis_dir, aligner="star", samples_info_path=None):
    """
    Validate alignment BAMs, read counts, and execution logs for all expected samples.
    Returns (is_valid: bool, report_str: str).
    """
    aligner = aligner.lower() if aligner else "star"
    expected_samples, rows_without_sample = parse_expected_samples(analysis_dir, samples_info_path)

    if rows_without_sample:
        return False, (
            f"[Validation ERROR] samples_info.txt has {rows_without_sample} row(s) with no sample name, "
            "which happens when the design lists more conditions than there are samples (for "
            "paired-end reads, one condition per sample, not one per FASTQ file). The sample "
            "conditions in that file are misaligned too, so please fix the design and rerun from "
            "the start rather than resuming.")

    if not expected_samples:
        if has_alignment_output(analysis_dir, aligner):
            return False, (
                "[Validation ERROR] Alignment output is present but the expected sample list "
                "could not be determined from samples_info.txt or raw_reads, so no sample could "
                "be checked. Refusing to continue into Step 4 with an unverifiable run.")
        return True, "[Validation Notice] No expected samples found in samples_info.txt or raw_reads. Skipping validation."

    failed_samples = {}
    successful_samples = []

    for sample in expected_samples:
        sample_errors = []
        log_excerpt = ""

        if aligner == "star":
            # 1. Check BAM file
            bam_pattern = os.path.join(analysis_dir, "miARma_out*", "star_results", f"{sample}_STAR.bam")
            bams = glob.glob(bam_pattern)
            if not bams:
                sample_errors.append(f"Missing BAM file (star_results/{sample}_STAR.bam)")
            else:
                bam_size = os.path.getsize(bams[0])
                if bam_size < 50000:  # < 50 KB means aborted/empty header-only BAM
                    sample_errors.append(f"Truncated/Empty BAM ({bam_size / 1024:.1f} KB)")

            # 2. Check readcount table
            rc_pattern = os.path.join(analysis_dir, "miARma_out*", "star_readcount_results", f"{sample}_star_readcount.tab")
            rcs = glob.glob(rc_pattern)
            if not rcs:
                sample_errors.append(f"Missing count table (star_readcount_results/{sample}_star_readcount.tab)")
            else:
                try:
                    with open(rcs[0], "r", encoding="utf-8", errors="replace") as f:
                        lines_count = sum(1 for _ in f)
                    if lines_count < 10:
                        sample_errors.append(f"Empty count table ({lines_count} lines)")
                except Exception as e:
                    sample_errors.append(f"Unreadable count table: {e}")

            # 3. Check STAR log for error messages
            star_log = os.path.join(analysis_dir, "miARma_out0", "star_results", f"{sample}_STAR_Log.out")
            log_excerpt = find_star_error(star_log)

        elif aligner == "hisat2":
            # 1. Check BAM file
            bam_pattern = os.path.join(analysis_dir, "miARma_out*", "hisat2_results", f"{sample}_hisat2.bam")
            bams = glob.glob(bam_pattern)
            if not bams:
                sample_errors.append(f"Missing BAM file (hisat2_results/{sample}_hisat2.bam)")
            else:
                bam_size = os.path.getsize(bams[0])
                if bam_size < 50000:
                    sample_errors.append(f"Truncated/Empty BAM ({bam_size / 1024:.1f} KB)")

            # 2. Check readcount table
            rc_pattern = os.path.join(analysis_dir, "miARma_out*", "hisat2_readcount_results", f"{sample}_hisat2_readcount.tab")
            rcs = glob.glob(rc_pattern)
            if not rcs:
                sample_errors.append(f"Missing count table (hisat2_readcount_results/{sample}_hisat2_readcount.tab)")
            else:
                try:
                    with open(rcs[0], "r", encoding="utf-8", errors="replace") as f:
                        lines_count = sum(1 for _ in f)
                    if lines_count < 10:
                        sample_errors.append(f"Empty count table ({lines_count} lines)")
                except Exception as e:
                    sample_errors.append(f"Unreadable count table: {e}")

            # 3. Check HISAT2 log
            hisat2_log = os.path.join(analysis_dir, "miARma_out0", "hisat2_results", f"{sample}_hisat2.log")
            log_excerpt = find_hisat2_error(hisat2_log)

        elif aligner == "bowtie1":
            # 1. Check BAM file
            bam_pattern = os.path.join(analysis_dir, "miARma_out*", "bowtie1_results", f"{sample}_bowtie1.bam")
            bams = glob.glob(bam_pattern)
            if not bams:
                sample_errors.append(f"Missing BAM file (bowtie1_results/{sample}_bowtie1.bam)")
            else:
                bam_size = os.path.getsize(bams[0])
                if bam_size < 5000:
                    sample_errors.append(f"Truncated/Empty BAM ({bam_size / 1024:.1f} KB)")

            # 2. Check readcount table
            rc_pattern = os.path.join(analysis_dir, "miARma_out*", "bowtie1_readcount_results", f"{sample}_bowtie1_readcount.tab")
            rcs = glob.glob(rc_pattern)
            if not rcs:
                sample_errors.append(f"Missing count table (bowtie1_readcount_results/{sample}_bowtie1_readcount.tab)")
            else:
                try:
                    with open(rcs[0], "r", encoding="utf-8", errors="replace") as f:
                        lines_count = sum(1 for _ in f)
                    if lines_count < 10:
                        sample_errors.append(f"Empty count table ({lines_count} lines)")
                except Exception as e:
                    sample_errors.append(f"Unreadable count table: {e}")

            # 3. Check bowtie1 log
            bowtie1_log = os.path.join(analysis_dir, "miARma_out0", "bowtie1_results", f"{sample}_bowtie1.log")
            log_excerpt = find_bowtie1_error(bowtie1_log)

        elif aligner == "kallisto":
            # 1. Check abundance.tsv
            ab_pattern = os.path.join(analysis_dir, "miARma_out*", "kallisto_results", sample, "abundance.tsv")
            abs_files = glob.glob(ab_pattern)
            if not abs_files:
                # Check alternate path
                ab_pattern2 = os.path.join(analysis_dir, "miARma_out*", "kallisto_results", f"{sample}_kallisto", "abundance.tsv")
                abs_files = glob.glob(ab_pattern2)

            if not abs_files:
                sample_errors.append(f"Missing abundance.tsv (kallisto_results/{sample}/abundance.tsv)")
            else:
                try:
                    with open(abs_files[0], "r", encoding="utf-8", errors="replace") as f:
                        lines_count = sum(1 for _ in f)
                    if lines_count < 10:
                        sample_errors.append(f"Empty abundance.tsv ({lines_count} lines)")
                except Exception as e:
                    sample_errors.append(f"Unreadable abundance.tsv: {e}")

            # 2. Check run_info.json
            json_pattern = os.path.join(analysis_dir, "miARma_out*", "kallisto_results", sample, "run_info.json")
            jsons = glob.glob(json_pattern)
            if jsons:
                try:
                    with open(jsons[0], "r", encoding="utf-8") as f:
                        data = json.load(f)
                    if data.get("n_pseudoaligned", 0) == 0 and data.get("n_processed", 0) > 0:
                        sample_errors.append("0 reads pseudoaligned (0.0%)")
                except Exception:
                    pass

            kallisto_log = os.path.join(analysis_dir, "miARma_out0", "kallisto_results", f"{sample}_kallisto.log")
            log_excerpt = find_kallisto_error(kallisto_log)

        if sample_errors:
            failed_samples[sample] = {
                "errors": sample_errors,
                "log_excerpt": log_excerpt
            }
        else:
            successful_samples.append(sample)

    # Check parallel job logs
    failed_parallel = check_parallel_logs(analysis_dir)

    is_valid = (len(failed_samples) == 0)

    # Build report string
    report_lines = []
    if is_valid:
        report_lines.append(f"\n[Validation OK] All {len(successful_samples)} expected sample(s) successfully completed {aligner.upper()} alignment and read counting.")
    else:
        report_lines.append("\n" + "=" * 80)
        report_lines.append(f"  ALIGNMENT / READCOUNT VALIDATION FAILED ({aligner.upper()})")
        report_lines.append("=" * 80)
        report_lines.append(f"Expected Samples   : {len(expected_samples)}")
        report_lines.append(f"Successful Samples : {len(successful_samples)}")
        report_lines.append(f"Failed Samples     : {len(failed_samples)} ({', '.join(sorted(failed_samples.keys()))})")
        report_lines.append("-" * 80)
        report_lines.append("Detailed Sample Failures:")
        for s, d in sorted(failed_samples.items()):
            report_lines.append(f"  * {s}:")
            for err in d["errors"]:
                report_lines.append(f"      - {err}")
            if d["log_excerpt"]:
                report_lines.append(f"      - Log Excerpt: {d['log_excerpt']}")

        if failed_parallel:
            report_lines.append("-" * 80)
            report_lines.append("Failed GNU Parallel Jobs:")
            for ptitle, pcmd in failed_parallel.items():
                report_lines.append(f"  * {ptitle}: {pcmd}")

        report_lines.append("=" * 80)
        report_lines.append(f"ERROR: Step 3b failed for {len(failed_samples)} sample(s).")
        report_lines.append("Pipeline execution halted before Step 4 to prevent cascaded errors.\n")

    return is_valid, "\n".join(report_lines)


def main():
    parser = argparse.ArgumentParser(description="Validate alignment and count tables across samples.")
    parser.add_argument("-d", "--dir", required=True, help="Analysis output directory ($output_folder/$name).")
    parser.add_argument("-a", "--aligner", default="star", help="Aligner used: star, hisat2, kallisto.")
    parser.add_argument("-s", "--samples", default=None, help="Path to samples_info.txt.")
    args = parser.parse_args()

    is_valid, report = validate_alignment_and_counts(args.dir, args.aligner, args.samples)
    print(report)

    if not is_valid:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
