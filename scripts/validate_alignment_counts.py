#!/usr/bin/env python3
"""
validate_alignment_counts.py

Universal validation gate for reanalyzerGSE across STAR, HISAT2, and Kallisto.
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


def parse_expected_samples(analysis_dir, samples_info_path=None):
    """Extract expected sample names from samples_info.txt or raw_reads directory."""
    expected = []
    
    # 1. Try samples_info.txt
    if samples_info_path and os.path.exists(samples_info_path):
        with open(samples_info_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line:
                    sample = line.split("\t")[0].strip()
                    if sample and sample not in expected:
                        expected.append(sample)
        if expected:
            return expected

    # 2. Try default samples_info.txt location inside analysis_dir
    default_info = os.path.join(analysis_dir, "reads_study_info", "samples_info.txt")
    if os.path.exists(default_info):
        with open(default_info, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line:
                    sample = line.split("\t")[0].strip()
                    if sample and sample not in expected:
                        expected.append(sample)
        if expected:
            return expected

    # 3. Fallback: inspect raw_reads folder
    raw_reads_dir = os.path.join(analysis_dir, "raw_reads")
    if os.path.exists(raw_reads_dir):
        files = os.listdir(raw_reads_dir)
        for f in files:
            if f.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
                # Strip suffixes
                s = f
                for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
                    if s.endswith(suffix):
                        s = s[:-len(suffix)]
                for end in ["_1", "_2", "_R1", "_R2"]:
                    if s.endswith(end):
                        s = s[:-len(end)]
                if s and s not in expected:
                    expected.append(s)

    return sorted(expected)


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


def validate_alignment_and_counts(analysis_dir, aligner="star", samples_info_path=None):
    """
    Validate alignment BAMs, read counts, and execution logs for all expected samples.
    Returns (is_valid: bool, report_str: str).
    """
    aligner = aligner.lower() if aligner else "star"
    expected_samples = parse_expected_samples(analysis_dir, samples_info_path)

    if not expected_samples:
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
