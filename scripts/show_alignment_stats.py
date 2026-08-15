#!/usr/bin/env python3
"""Display and save a formatted alignment summary table for STAR, HISAT2, or Kallisto.

Parses aligner logs and displays:
  - Sample Name
  - Total Input Reads
  - Uniquely Mapped Reads & Percentage
  - Multi-Mapped Reads & Percentage
  - Overall Alignment Rate (with ANSI color coding)
"""

import argparse
import glob
import json
import os
import re
import sys


def format_number(num):
    """Format numbers cleanly (e.g. 1.25 M, 450.2 K, or exact count)."""
    if num >= 1e6:
        return f"{num / 1e6:.2f} M"
    elif num >= 1e3:
        return f"{num / 1e3:.1f} K"
    else:
        return f"{int(num)}"


def parse_alignment_stats(analysis_dir, aligner=""):
    """Parse log files from STAR, HISAT2, or Kallisto and return structured stats."""
    sample_stats = {}  # sample_name -> dict

    # 1. Search for HISAT2 logs
    hisat2_logs = glob.glob(f"{analysis_dir}/**/hisat2_results/*_hisat2.log", recursive=True)
    if not hisat2_logs:
        hisat2_logs = [f for f in glob.glob(f"{analysis_dir}/**/*_hisat2.log", recursive=True) if "sphinx" not in f]

    # 2. Search for STAR logs
    star_logs = glob.glob(f"{analysis_dir}/**/star_results/*Log.final.out", recursive=True)
    if not star_logs:
        star_logs = [f for f in glob.glob(f"{analysis_dir}/**/*Log.final.out", recursive=True) if "sphinx" not in f]

    # 3. Search for Kallisto logs
    kallisto_logs = glob.glob(f"{analysis_dir}/**/kallisto_results/*/run_info.json", recursive=True)
    if not kallisto_logs:
        kallisto_logs = [f for f in glob.glob(f"{analysis_dir}/**/run_info.json", recursive=True) if "sphinx" not in f]

    detected_aligner = aligner.upper() if aligner else "ALIGNMENT"

    if hisat2_logs and (not aligner or aligner.lower() == "hisat2"):
        detected_aligner = "HISAT2"
        for log_path in hisat2_logs:
            sample = os.path.basename(log_path).replace("_hisat2.log", "").replace(".log", "")
            if sample in sample_stats:
                continue

            total_reads = 0
            unique_reads = 0
            unique_pct = 0.0
            multi_reads = 0
            multi_pct = 0.0
            overall_rate = 0.0

            with open(log_path, "r", errors="ignore") as f:
                content = f.read()

            # Total pairs / reads
            m = re.search(r'Total (?:pairs|reads):\s*(\d+)', content)
            if m:
                total_reads = int(m.group(1))

            # Overall alignment rate
            m = re.search(r'Overall alignment rate:\s*([\d\.]+)%', content)
            if not m:
                m = re.search(r'([\d\.]+)%\s*overall alignment rate', content)
            if m:
                overall_rate = float(m.group(1))

            # Unique alignment
            m = re.search(r'Aligned (?:concordantly )?1 time:\s*(\d+)\s*\(([\d\.]+)%\)', content)
            if m:
                unique_reads = int(m.group(1))
                unique_pct = float(m.group(2))

            # Multi-mapping
            m = re.search(r'Aligned (?:concordantly )?>1 times:\s*(\d+)\s*\(([\d\.]+)%\)', content)
            if m:
                multi_reads = int(m.group(1))
                multi_pct = float(m.group(2))

            # Legacy HISAT2 summary format fallback
            if total_reads == 0:
                m = re.search(r'(\d+)\s*reads; of these:', content)
                if m:
                    total_reads = int(m.group(1))
                m = re.search(r'(\d+)\s*\(([\d\.]+)%\)\s*aligned exactly 1 time', content)
                if m:
                    unique_reads = int(m.group(1))
                    unique_pct = float(m.group(2))
                m = re.search(r'(\d+)\s*\(([\d\.]+)%\)\s*aligned >1 times', content)
                if m:
                    multi_reads = int(m.group(1))
                    multi_pct = float(m.group(2))

            sample_stats[sample] = {
                "sample": sample,
                "total": total_reads,
                "unique": unique_reads,
                "unique_pct": unique_pct,
                "multi": multi_reads,
                "multi_pct": multi_pct,
                "overall_pct": overall_rate,
            }

    elif star_logs and (not aligner or aligner.lower() == "star"):
        detected_aligner = "STAR"
        for log_path in star_logs:
            sample = os.path.basename(log_path).replace("_STAR_Log.final.out", "").replace("_Log.final.out", "").replace("Log.final.out", "")
            if sample in sample_stats:
                continue

            total_reads = 0
            unique_reads = 0
            unique_pct = 0.0
            multi_reads = 0
            multi_pct = 0.0

            with open(log_path, "r", errors="ignore") as f:
                for line in f:
                    if "Number of input reads" in line:
                        total_reads = int(line.split("|")[1].strip())
                    elif "Uniquely mapped reads number" in line:
                        unique_reads = int(line.split("|")[1].strip())
                    elif "Uniquely mapped reads %" in line:
                        unique_pct = float(line.split("|")[1].replace("%", "").strip())
                    elif "Number of reads mapped to multiple loci" in line:
                        multi_reads = int(line.split("|")[1].strip())
                    elif "% of reads mapped to multiple loci" in line:
                        multi_pct = float(line.split("|")[1].replace("%", "").strip())

            overall_rate = unique_pct + multi_pct
            sample_stats[sample] = {
                "sample": sample,
                "total": total_reads,
                "unique": unique_reads,
                "unique_pct": unique_pct,
                "multi": multi_reads,
                "multi_pct": multi_pct,
                "overall_pct": overall_rate,
            }

    elif kallisto_logs and (not aligner or aligner.lower() == "kallisto"):
        detected_aligner = "Kallisto"
        for log_path in kallisto_logs:
            sample = os.path.basename(os.path.dirname(log_path))
            if sample in sample_stats:
                continue
            try:
                with open(log_path, "r") as f:
                    info = json.load(f)
                total_reads = int(info.get("n_processed", 0))
                unique_reads = int(info.get("n_unique", 0))
                unique_pct = float(info.get("p_unique", 0.0))
                pseudo_reads = int(info.get("n_pseudoaligned", 0))
                pseudo_pct = float(info.get("p_pseudoaligned", 0.0))
                multi_reads = max(0, pseudo_reads - unique_reads)
                multi_pct = max(0.0, pseudo_pct - unique_pct)

                sample_stats[sample] = {
                    "sample": sample,
                    "total": total_reads,
                    "unique": unique_reads,
                    "unique_pct": unique_pct,
                    "multi": multi_reads,
                    "multi_pct": multi_pct,
                    "overall_pct": pseudo_pct,
                }
            except Exception:
                pass

    return detected_aligner, list(sample_stats.values())


def colorize(text, ansi_code, use_color=True):
    """Wrap text in ANSI color escape codes only if use_color is True."""
    if use_color:
        return f"\033[{ansi_code}m{text}\033[0m"
    return text


def display_and_save_stats(analysis_dir, aligner="", save_path=None, force_no_color=False):
    """Print a clean terminal table and optionally save to TSV."""
    detected_aligner, stats = parse_alignment_stats(analysis_dir, aligner)

    if not stats:
        print(f"\n[Alignment Summary] No alignment log files found in {analysis_dir}.")
        return

    # Check if stdout is an interactive terminal and no color-suppression is requested
    use_color = sys.stdout.isatty() and not force_no_color and "NO_COLOR" not in os.environ

    # Sort samples alphabetically
    stats = sorted(stats, key=lambda x: x["sample"])

    # Header
    print("\n" + "=" * 80)
    print(f"  ALIGNMENT SUMMARY TABLE ({detected_aligner})")
    print("=" * 80)
    print(f"{'Sample':<25} {'Total Reads':<14} {'Uniquely Mapped':<20} {'Multi-Mapped':<18} {'Overall Rate':<12}")
    print("-" * 80)

    total_input = 0
    total_uniq = 0
    total_multi = 0

    for s in stats:
        s_name = s["sample"][:23]
        tot_str = format_number(s["total"])
        uniq_str = f"{format_number(s['unique'])} ({s['unique_pct']:.1f}%)"
        multi_str = f"{format_number(s['multi'])} ({s['multi_pct']:.1f}%)"

        rate = s["overall_pct"]
        rate_clean = f"{rate:6.2f}%"

        # Color coding: Green >= 70%, Yellow 40-70%, Red < 40%
        if rate >= 70.0:
            rate_fmt = colorize(rate_clean, "1;32", use_color)
        elif rate >= 40.0:
            rate_fmt = colorize(rate_clean, "1;33", use_color)
        else:
            rate_fmt = colorize(rate_clean, "1;31", use_color)

        print(f"{s_name:<25} {tot_str:<14} {uniq_str:<20} {multi_str:<18} {rate_fmt}")
        total_input += s["total"]
        total_uniq += s["unique"]
        total_multi += s["multi"]

    print("-" * 80)
    avg_rate = ((total_uniq + total_multi) / total_input * 100) if total_input > 0 else 0.0
    tot_in_str = format_number(total_input)
    tot_un_str = f"{format_number(total_uniq)} ({total_uniq / total_input * 100:.1f}%)" if total_input > 0 else "0"
    tot_mu_str = f"{format_number(total_multi)} ({total_multi / total_input * 100:.1f}%)" if total_input > 0 else "0"
    avg_clean = f"{avg_rate:6.2f}%"
    avg_col = colorize(avg_clean, "1;32" if avg_rate >= 70 else "1;33" if avg_rate >= 40 else "1;31", use_color)
    print(f"{'TOTAL / MEAN':<25} {tot_in_str:<14} {tot_un_str:<20} {tot_mu_str:<18} {avg_col}")
    print("=" * 80 + "\n")

    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write("Sample\tTotal_Reads\tUniquely_Mapped_Reads\tUniquely_Mapped_Pct\tMulti_Mapped_Reads\tMulti_Mapped_Pct\tOverall_Rate_Pct\n")
            for s in stats:
                f.write(f"{s['sample']}\t{s['total']}\t{s['unique']}\t{s['unique_pct']:.2f}\t{s['multi']}\t{s['multi_pct']:.2f}\t{s['overall_pct']:.2f}\n")
            f.write(f"TOTAL_MEAN\t{total_input}\t{total_uniq}\t{(total_uniq/total_input*100) if total_input>0 else 0:.2f}\t{total_multi}\t{(total_multi/total_input*100) if total_input>0 else 0:.2f}\t{avg_rate:.2f}\n")


def main():
    parser = argparse.ArgumentParser(description="Display alignment percentage table per sample.")
    parser.add_argument("-d", "--dir", required=True, help="Analysis output directory.")
    parser.add_argument("-a", "--aligner", default="", help="Aligner name: star, hisat2, kallisto.")
    parser.add_argument("-s", "--save", default=None, help="Optional output TSV file path.")
    parser.add_argument("--no-color", action="store_true", help="Disable ANSI color output.")
    args = parser.parse_args()

    display_and_save_stats(args.dir, args.aligner, args.save, args.no_color)


if __name__ == "__main__":
    main()
