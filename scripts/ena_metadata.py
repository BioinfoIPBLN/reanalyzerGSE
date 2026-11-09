#!/usr/bin/env python3
"""ENA metadata for the accessions given to reanalyzerGSE in -i (PRJ*, [EDS]RP, SAM*, [EDS]RS, [EDS]RX, [EDS]RR).

--out-dir writes the files that the pipeline reads in reads_study_info; --print only shows the table;
--design-auto reads sample names from stdin and prints the design derived from them.
"""
import argparse
import csv
import io
import os
import re
import sys
import time
import urllib.parse
import urllib.request

FIELDS = ["run_accession", "sample_accession", "secondary_sample_accession", "experiment_accession",
          "study_accession", "secondary_study_accession", "sample_title", "sample_alias", "experiment_title",
          "experiment_alias", "run_alias", "sample_description", "description", "scientific_name", "tax_id",
          "library_layout", "library_strategy", "library_source", "library_selection", "instrument_platform",
          "instrument_model", "read_count", "base_count", "nominal_length", "fastq_ftp", "study_title", "first_public"]
RNA_STRATEGIES = {"RNA-Seq", "miRNA-Seq", "ncRNA-Seq", "ssRNA-seq", "FL-cDNA", "EST"}
ACCESSION_RE = re.compile(r"^(PRJ(EB|NA|DB)|[EDS]RP|SAM[DEN]A?|[EDS]R[SXR])[0-9]+$")
DB_TOKEN_RE = re.compile(r"^(GSM|GSE|SRR|ERR|DRR|SRS|ERS|DRS|SRX|ERX|DRX|SRP|ERP|DRP|SAMN|SAMEA|SAMD|PRJNA|PRJEB|PRJDB)[0-9]+$", re.I)
REPLICATE_RE = re.compile(r"(?:[ _.\-]+(?:replicate|rep|r|n)?|(?:replicate|rep)|)[ _.\-]*[0-9]+$", re.I)


def fetch_runs(accession, retries=4):
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
        {"accession": accession, "result": "read_run", "fields": ",".join(FIELDS), "format": "tsv", "limit": 0})
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "reanalyzerGSE ena_metadata.py"})
            with urllib.request.urlopen(req, timeout=180) as r:
                text = r.read().decode("utf-8", "replace")
            return list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt < retries - 1:
                time.sleep(15 * (attempt + 1))
    raise RuntimeError(f"ENA did not answer for {accession}: {last}")


def clean_name(text):
    s = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_")
    s = re.sub(r"_+", "_", s)
    s = re.sub(r"_([0-9]+)", r"-\1", s)
    if re.match(r"^[0-9]", s):
        s = "S" + s
    return s[:80] if s else ""


def sample_base_name(rows, organism):
    title = (rows[0].get("sample_title") or "").strip()
    gsm = ""
    m = re.match(r"^(GSM[0-9]+):\s*(.*)$", title)
    if m:
        gsm, title = m.group(1), m.group(2)
    strategy = (rows[0].get("library_strategy") or "").strip()
    for tail in (organism, strategy):
        if tail:
            title = re.sub(r";?\s*" + re.escape(tail) + r"\s*;?\s*$", "", title, flags=re.I).strip()
    if not title:
        title = (rows[0].get("sample_alias") or "").strip()
    if not title or re.match(r"^(SAM[DEN]A?|[EDS]RS)[0-9]+$", title):
        title = (rows[0].get("experiment_title") or "").strip()
    if not title:
        title = rows[0]["run_accession"]
    name = clean_name(title)
    if gsm:
        name = f"{name}_{gsm}" if name else gsm
    return name or rows[0]["run_accession"]


def group_samples(runs, organism):
    by_sample = {}
    for r in runs:
        key = r.get("sample_accession") or r.get("secondary_sample_accession") or r["run_accession"]
        by_sample.setdefault(key, []).append(r)
    samples = []
    for key, rows in by_sample.items():
        rows.sort(key=lambda x: x["run_accession"])
        samples.append({"sample_accession": key, "rows": rows, "base": sample_base_name(rows, organism)})
    samples.sort(key=lambda s: (s["base"], s["sample_accession"]))
    seen = {}
    for s in samples:
        seen[s["base"]] = seen.get(s["base"], 0) + 1
        s["name"] = s["base"] if seen[s["base"]] == 1 else f"{s['base']}-{seen[s['base']]}"
    samples.sort(key=lambda s: s["name"])
    return samples


def design_auto(names):
    conds = []
    for n in names:
        tokens = [t for t in re.split(r"[_\-]", n) if t and not DB_TOKEN_RE.match(t)]
        core = "_".join(tokens) if tokens else n
        cond = REPLICATE_RE.sub("", core)
        conds.append(cond if cond else core)
    groups = {}
    for c in conds:
        groups[c] = groups.get(c, 0) + 1
    if len(groups) < 2 or max(groups.values()) < 2:
        sys.stderr.write("design 'auto': dropping the replicate numbers of the sample names gives " +
                         ", ".join(f"{k} ({v})" for k, v in groups.items()) +
                         ", which is not a usable design (at least two conditions, one of them with replicates, are needed)\n")
        return None
    return conds


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", help="comma-separated ENA/SRA accessions")
    ap.add_argument("--out-dir", help="reads_study_info folder to write the pipeline files into")
    ap.add_argument("--organism", default="", help="scientific name that overrides the one reported by ENA")
    ap.add_argument("--print", action="store_true", help="only print the runs and samples")
    ap.add_argument("--design-auto", action="store_true", help="read sample names from stdin and print the auto design")
    a = ap.parse_args()

    if a.design_auto:
        names = [l.strip() for l in sys.stdin if l.strip()]
        conds = design_auto(names)
        if conds is None:
            sys.exit(2)
        print(",".join(conds))
        return

    if not a.accessions:
        ap.error("--accessions is required")
    accessions = [x.strip() for x in a.accessions.split(",") if x.strip()]
    bad = [x for x in accessions if not ACCESSION_RE.match(x)]
    if bad:
        sys.stderr.write("ERROR: not ENA/SRA accessions: " + ", ".join(bad) + "\n")
        sys.exit(1)
    runs = []
    seen_runs = set()
    for acc in accessions:
        rows = fetch_runs(acc)
        if not rows:
            sys.stderr.write(f"WARNING: ENA has no runs for {acc}\n")
        for r in rows:
            if r["run_accession"] not in seen_runs:
                seen_runs.add(r["run_accession"])
                runs.append(r)
    if not runs:
        sys.stderr.write("ERROR: no runs found in ENA for " + ", ".join(accessions) + "\n")
        sys.exit(1)

    strategies = sorted({(r.get("library_strategy") or "").strip() for r in runs})
    rna = [r for r in runs if (r.get("library_strategy") or "").strip() in RNA_STRATEGIES]
    if rna and len(rna) < len(runs):
        dropped = [r for r in runs if r not in rna]
        sys.stderr.write(f"NOTE: {len(dropped)} run(s) with a non-RNA library strategy are left out: " +
                         ", ".join(f"{r['run_accession']} ({r.get('library_strategy')})" for r in dropped) + "\n")
        runs = rna
    elif not rna:
        sys.stderr.write("WARNING: none of the runs is RNA-seq according to ENA (library_strategy: " +
                         ", ".join(strategies) + "); all of them are kept\n")

    organisms = sorted({(r.get("scientific_name") or "").strip() for r in runs if (r.get("scientific_name") or "").strip()})
    organism = a.organism.strip()
    if not organism:
        if len(organisms) == 1:
            organism = organisms[0]
        elif len(organisms) > 1:
            sys.stderr.write("ERROR: the runs come from several organisms (" + ", ".join(organisms) +
                             "). Give the accessions of one organism, or set 'organism' to the one to analyse.\n")
            sys.exit(1)
        else:
            sys.stderr.write("ERROR: ENA reports no organism for these runs; please set 'organism'.\n")
            sys.exit(1)
    layouts = sorted({(r.get("library_layout") or "").strip().upper() for r in runs})
    if len(layouts) != 1 or layouts[0] not in ("SINGLE", "PAIRED"):
        sys.stderr.write("ERROR: the runs mix library layouts (" + ", ".join(layouts) +
                         "). reanalyzerGSE needs one layout per run; give the run accessions of one layout.\n")
        sys.exit(1)
    layout = layouts[0]
    samples = group_samples(runs, organism)

    header = ["sample_name", "runs", "layout", "scientific_name", "library_strategy", "sample_accession", "sample_title"]
    table = []
    for s in samples:
        r0 = s["rows"][0]
        table.append([s["name"], ";".join(r["run_accession"] for r in s["rows"]), layout,
                      (r0.get("scientific_name") or "").strip(), (r0.get("library_strategy") or "").strip(),
                      s["sample_accession"], (r0.get("sample_title") or "").strip()])
    multi = [t for t in table if ";" in t[1]]
    out = sys.stdout
    if not a.out_dir:
        out.write("\t".join(header) + "\n")
        for t in table:
            out.write("\t".join(t) + "\n")
        out.write("\n")
    out.write(f"{len(runs)} run(s) in {len(samples)} sample(s); organism: {organism}; layout: {layout}\n")
    if multi:
        out.write("The runs of a sample are merged into one FASTQ per mate: " +
                  ", ".join(f"{t[0]} ({t[1].count(';') + 1} runs)" for t in multi) + "\n")
    if a.out_dir:
        os.makedirs(a.out_dir, exist_ok=True)
        with open(os.path.join(a.out_dir, "ena_run_metadata.tsv"), "w") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in runs:
                w.writerow(r)
        with open(os.path.join(a.out_dir, "ena_samples.tsv"), "w") as fh:
            fh.write("\t".join(header) + "\n")
            for t in table:
                fh.write("\t".join(t) + "\n")
        with open(os.path.join(a.out_dir, "srr_ids.txt"), "w") as fh:
            for s in samples:
                for r in s["rows"]:
                    fh.write(r["run_accession"] + "\n")
        with open(os.path.join(a.out_dir, "sample_names.txt"), "w") as fh:
            for s in samples:
                fh.write(s["name"] + "\n")
        with open(os.path.join(a.out_dir, "library_layout_info.txt"), "w") as fh:
            fh.write(layout + "\n")
        with open(os.path.join(a.out_dir, "organism.txt"), "w") as fh:
            fh.write(organism + "\n")
        titles = sorted({(r.get("study_title") or "").strip() for r in runs if (r.get("study_title") or "").strip()})
        with open(os.path.join(a.out_dir, "study_title.txt"), "w") as fh:
            fh.write((titles[0] if titles else ",".join(accessions)) + "\n")


if __name__ == "__main__":
    main()
