#!/usr/bin/env python3
"""Detect human orthologs of the analysed organism's proteins.

Wraps orthologr's ortho_detection methods (DIAMOND_RBH, DIAMOND_BH, RBH, BH) via
R_orthologs_human.R. Resolves the human reference proteome, maps the query proteins back to
gene identifiers using the pipeline annotation, and writes a gene-level query->human symbol
table for downstream use.

Outputs in --out-dir:
    orthologs_human.tsv   one row per ortholog pair, protein and gene level
    gene_to_human.tsv     best human symbol per query gene identifier or alias
    gene_to_protein.tsv   query gene -> protein identifiers used
    status.json           counts and the parameters used

Exit codes:
    0   ortholog pairs written
    3   ran cleanly but nothing usable (no pairs, or no gene could be mapped)
    4   a required external resource was unavailable (proteome download, tool missing)
    5   the query proteome or the annotation could not be parsed
"""

import argparse, collections, gzip, json, os, re, shutil, subprocess, sys, urllib.error, urllib.request

ORTHOLOGR_METHODS = ("DIAMOND_RBH", "DIAMOND_BH", "RBH", "BH")
UNIPROT_HUMAN = ("https://rest.uniprot.org/uniprotkb/stream"
                 "?format=fasta&compressed=true&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29")
ID_PREFIX_RE = re.compile(r"^(?:lcl\||gnl\|[^|]*\||ref\||sp\||tr\|)")
SYMBOL_KEYS = ("GN", "gene_symbol", "gene_name", "gene", "locus_tag")
GENE_ATTR_KEYS = ("gene_name", "gene", "gene_id", "Name", "locus_tag")
FEATURES_WITH_PROTEIN = ("CDS", "mRNA", "transcript", "exon")
VERSION_RE = re.compile(r"\.\d+$")
DUP_RE = re.compile(r"__dup\d+$")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def opener(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") \
        else open(path, "rt", errors="replace")


def read_fasta(path):
    header, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].rstrip("\n"), []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def base_id(pid):
    return DUP_RE.sub("", pid)


def strip_version(pid):
    return VERSION_RE.sub("", pid)


def clean_id(raw):
    tok = raw.split()[0] if raw.strip() else ""
    tok = ID_PREFIX_RE.sub("", tok)
    if "|" in tok:
        parts = [p for p in tok.split("|") if p]
        if parts:
            tok = parts[0]
    return tok


def attr_value(text, key, sep_hint=None):
    """Value of `key` in a GTF (key "value") / GFF (key=value) / FASTA (key=value or key:value) blob."""
    seps = sep_hint or "[ =:]"
    m = re.search(r'(?:^|[;\s])' + re.escape(key) + seps + r'+"?([^";\n]+?)"?\s*(?:;|$|\s)', text)
    return m.group(1).strip() if m else ""


def header_symbols(header):
    """Candidate gene identifiers in a FASTA header, most authoritative first."""
    out = []
    for key in SYMBOL_KEYS:
        val = attr_value(header, key)
        if val and val not in out:
            out.append(val)
    m = re.search(r"\[gene=([^\]]+)\]", header)
    if m and m.group(1).strip() not in out:
        out.insert(0, m.group(1).strip())
    return out


def with_aliases(names):
    out = []
    for n in names:
        for cand in (n, strip_version(n)):
            if cand and cand not in out:
                out.append(cand)
    return out


def sanitize_proteome(src, dest):
    """Write a FASTA with one bare identifier per header and BLAST/DIAMOND-safe residues."""
    seen, n_in, n_out = {}, 0, 0
    headers, genes = {}, {}
    with open(dest, "w", encoding="utf-8") as out:
        for header, seq in read_fasta(src):
            n_in += 1
            pid = clean_id(header)
            if not pid:
                continue
            seq = seq.upper().rstrip("*").replace("*", "X").replace(".", "X").replace("-", "")
            if not seq:
                continue
            if pid in seen:
                seen[pid] += 1
                pid = f"{pid}__dup{seen[pid]}"
            else:
                seen[pid] = 0
            headers[pid] = header
            cands = header_symbols(header)
            if cands:
                genes[pid] = cands
            out.write(f">{pid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
            n_out += 1
    log(f"  {os.path.basename(src)}: {n_in} record(s) read, {n_out} written to {os.path.basename(dest)}"
        f", {len(genes)} carrying a gene identifier in the header")
    return headers, genes


def download_human_proteome(cache_dir):
    dest = os.path.join(cache_dir, "human_reference_proteome_uniprot.faa.gz")
    if os.path.isfile(dest) and os.path.getsize(dest) > 1_000_000:
        log(f"  reusing the cached human proteome: {dest}")
        return dest
    os.makedirs(cache_dir, exist_ok=True)
    tmp = dest + ".part"
    log(f"  downloading the human reference proteome (UniProt UP000005640, reviewed) to {dest}")
    try:
        req = urllib.request.Request(UNIPROT_HUMAN, headers={"User-Agent": "reanalyzerGSE"})
        with urllib.request.urlopen(req, timeout=900) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
    except (OSError, urllib.error.URLError) as exc:
        if os.path.exists(tmp):
            os.remove(tmp)
        log(f"  [error] the human proteome could not be downloaded: {exc}")
        return None
    if os.path.getsize(tmp) < 1_000_000:
        os.remove(tmp)
        log("  [error] the downloaded human proteome is implausibly small; discarded.")
        return None
    os.replace(tmp, dest)
    return dest


def parse_annotation(path, preferred_attr):
    """protein_id -> ordered candidate gene identifiers, from a GTF/GFF/GFF3 (plain or gzipped)."""
    p2g = {}
    if not path or not os.path.isfile(path):
        log("  annotation: not provided or unreadable; falling back to the FASTA headers only")
        return p2g
    keys = [preferred_attr] + [k for k in GENE_ATTR_KEYS if k != preferred_attr]
    rows = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in FEATURES_WITH_PROTEIN:
                continue
            attrs = cols[8]
            prot = attr_value(attrs, "protein_id")
            if not prot:
                continue
            cands = []
            for key in keys:
                val = attr_value(attrs, key)
                if val and val not in cands:
                    cands.append(val)
            if not cands:
                continue
            rows += 1
            prev = p2g.setdefault(prot, [])
            for c in cands:
                if c not in prev:
                    prev.append(c)
    log(f"  annotation: {len(p2g)} protein_id(s) linked to a gene identifier over {rows} row(s)")
    return p2g


def build_protein_to_gene(query_headers, header_genes, annotation_path, preferred_attr):
    p2g_annot = parse_annotation(annotation_path, preferred_attr)
    lookup = {}
    for k, v in p2g_annot.items():
        lookup.setdefault(k, v)
        lookup.setdefault(strip_version(k), v)
    out, from_annot, from_header = {}, 0, 0
    for pid in query_headers:
        cands = lookup.get(base_id(pid)) or lookup.get(strip_version(pid))
        if cands:
            out[pid] = with_aliases(cands)
            from_annot += 1
            continue
        cands = header_genes.get(pid)
        if cands:
            out[pid] = with_aliases(cands)
            from_header += 1
    log(f"  query proteins linked to a gene: {from_annot} via the annotation, "
        f"{from_header} via the FASTA header, {len(query_headers) - len(out)} unlinked")
    return out


def run_orthologr(script_dir, query, subject, method, cores, out_tsv, work_dir, eval_thr, sens):
    rscript = os.path.join(script_dir, "R_orthologs_human.R")
    cmd = ["Rscript", rscript, query, subject, method, str(cores), out_tsv, work_dir, eval_thr, sens]
    log("  " + " ".join(cmd))
    return subprocess.call(cmd)


def read_pairs_tsv(path):
    pairs = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            if not ln.strip():
                continue
            row = dict(zip(header, ln.rstrip("\n").split("\t")))
            if row.get("query_id") and row.get("subject_id"):
                pairs.append((row["query_id"], row["subject_id"], row.get("perc_identity", ""),
                              row.get("evalue", ""), row.get("bit_score", "")))
    return pairs


def pair_score(perc_identity, bit_score):
    for v in (perc_identity, bit_score):
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return 0.0


def main():
    p = argparse.ArgumentParser(description="Detect human orthologs of the analysed organism's proteins.")
    p.add_argument("--faa", required=True, help="protein FASTA of the analysed organism (.faa or .faa.gz)")
    p.add_argument("--method", required=True, choices=list(ORTHOLOGR_METHODS))
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cache-dir", required=True,
                   help="reused across runs for the human proteome and the search databases")
    p.add_argument("--organism", default="")
    p.add_argument("--annotation", default="", help="pipeline GTF/GFF, used to map proteins back to genes")
    p.add_argument("--gene-attribute", default="gene_name",
                   help="annotation attribute holding the identifiers used in the count tables "
                        "(matches the pipeline's optionsFeatureCounts_seq; default gene_name)")
    p.add_argument("--human-faa", default="", help="human reference proteome; downloaded from UniProt when empty")
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--eval", dest="eval_thr", default="1E-5")
    p.add_argument("--sensitivity-mode", default="fast",
                   help="DIAMOND sensitivity: fast, mid-sensitive, sensitive, more-sensitive, "
                        "very-sensitive, ultra-sensitive")
    p.add_argument("--min-identity", type=float, default=0.0, help="drop pairs below this percentage identity")
    p.add_argument("--keep-work", action="store_true", help="do not delete the search databases and intermediates")
    p.add_argument("--force", action="store_true",
                   help="re-run the search even when cached ortholog pairs for this method exist")
    a = p.parse_args()

    script_dir = os.path.dirname(os.path.realpath(__file__))
    a.out_dir = os.path.abspath(a.out_dir)
    a.cache_dir = os.path.abspath(a.cache_dir)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.cache_dir, exist_ok=True)
    work_dir = os.path.join(a.cache_dir, f"work_{a.method}")
    os.makedirs(work_dir, exist_ok=True)

    if not os.path.isfile(a.faa):
        log(f"ERROR: the protein FASTA '{a.faa}' does not exist.")
        return 5

    log("Preparing the proteomes ...")
    query_clean = os.path.join(work_dir, "query.faa")
    query_headers, query_header_genes = sanitize_proteome(a.faa, query_clean)
    if not query_headers:
        log("ERROR: no protein record could be read from the provided FASTA.")
        return 5

    human_src = a.human_faa or download_human_proteome(a.cache_dir)
    if not human_src or not os.path.isfile(human_src):
        log("ERROR: no human reference proteome is available.")
        return 4
    human_clean = os.path.join(work_dir, "human.faa")
    human_headers, human_genes = sanitize_proteome(human_src, human_clean)
    if not human_headers:
        log("ERROR: no protein record could be read from the human proteome.")
        return 4
    if not human_genes:
        log("WARNING: no gene symbol could be parsed from the human proteome headers. The ENCODE "
            "lookup needs HGNC symbols, so please provide a proteome carrying 'GN=', 'gene_symbol:' "
            "or '[gene=]' in its headers.")

    log(f"Linking the query proteins to gene identifiers (organism: {a.organism or 'unspecified'}) ...")
    p2g = build_protein_to_gene(query_headers, query_header_genes, a.annotation, a.gene_attribute)

    pairs_tsv = os.path.join(a.out_dir, f"orthologr_pairs_{a.method}.tsv")
    cached_pairs = os.path.join(a.cache_dir, f"orthologr_pairs_{a.method}.tsv")
    pairs = None
    if not a.force and os.path.isfile(cached_pairs) and os.path.getsize(cached_pairs) > 0:
        pairs = read_pairs_tsv(cached_pairs)
        log(f"Reusing {len(pairs)} cached ortholog pair(s) from {cached_pairs} "
            "(pass --force to search again).")
        if not pairs:
            pairs = None

    if pairs is None:
        log(f"Running ortholog detection with {a.method} on {a.cores} core(s) ...")
        rc = run_orthologr(script_dir, query_clean, human_clean, a.method, a.cores,
                           pairs_tsv, work_dir, a.eval_thr, a.sensitivity_mode)
        if rc == 3:
            log("ERROR: the orthologr R package is missing from this installation.")
            return 4
        if rc != 0 or not os.path.isfile(pairs_tsv):
            log(f"ERROR: ortholog detection failed (exit {rc}).")
            return 4 if rc == 4 else 3
        pairs = read_pairs_tsv(pairs_tsv)
        try:
            with open(cached_pairs, "w", encoding="utf-8") as fh:
                fh.write("query_id\tsubject_id\tperc_identity\tevalue\tbit_score\n")
                for row in pairs:
                    fh.write("\t".join(str(x) for x in row) + "\n")
        except OSError as exc:
            log(f"  [warn] the ortholog pairs could not be cached: {exc}")

    if a.min_identity > 0:
        kept = []
        for q, s, pid, ev, bs in pairs:
            try:
                if float(pid) < a.min_identity:
                    continue
            except (TypeError, ValueError):
                pass
            kept.append((q, s, pid, ev, bs))
        log(f"  identity filter (>= {a.min_identity}%): {len(kept)}/{len(pairs)} pair(s) kept")
        pairs = kept

    status_path = os.path.join(a.out_dir, "status.json")
    if not pairs:
        log("No ortholog pair was found.")
        with open(status_path, "w", encoding="utf-8") as fh:
            json.dump({"organism": a.organism, "method": a.method, "pairs": 0,
                       "genes_with_ortholog": 0, "usable": False}, fh, indent=1)
        return 3

    rows = []
    for q, s, pid, ev, bs in pairs:
        q_genes = p2g.get(q) or p2g.get(base_id(q)) or []
        h_genes = human_genes.get(s) or human_genes.get(base_id(s)) or []
        rows.append({
            "query_protein": base_id(q),
            "query_gene": q_genes[0] if q_genes else "",
            "query_gene_aliases": ",".join(q_genes[1:]),
            "human_protein": base_id(s),
            "human_symbol": h_genes[0] if h_genes else "",
            "perc_identity": pid, "evalue": ev, "bit_score": bs, "method": a.method,
            "_q_genes": q_genes,
        })

    cols = ["query_protein", "query_gene", "query_gene_aliases", "human_protein", "human_symbol",
            "perc_identity", "evalue", "bit_score", "method"]
    with open(os.path.join(a.out_dir, "orthologs_human.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    best = {}
    for r in rows:
        if not r["human_symbol"]:
            continue
        sc = pair_score(r["perc_identity"], r["bit_score"])
        for gene in r["_q_genes"]:
            cur = best.get(gene)
            if cur is None or sc > cur[0]:
                best[gene] = (sc, r)
    with open(os.path.join(a.out_dir, "gene_to_human.tsv"), "w", encoding="utf-8") as fh:
        fh.write("query_gene\thuman_symbol\thuman_protein\tquery_protein\tperc_identity\tmethod\n")
        for gene in sorted(best):
            r = best[gene][1]
            fh.write("\t".join([gene, r["human_symbol"], r["human_protein"], r["query_protein"],
                                str(r["perc_identity"]), r["method"]]) + "\n")

    g2p = collections.defaultdict(set)
    for pid, cands in p2g.items():
        if cands:
            g2p[cands[0]].add(base_id(pid))
    with open(os.path.join(a.out_dir, "gene_to_protein.tsv"), "w", encoding="utf-8") as fh:
        fh.write("query_gene\tquery_proteins\n")
        for gene in sorted(g2p):
            fh.write(f"{gene}\t{','.join(sorted(g2p[gene]))}\n")

    n_genes = len({r["query_gene"] for r in rows if r["query_gene"]})
    with open(status_path, "w", encoding="utf-8") as fh:
        json.dump({"organism": a.organism, "method": a.method, "eval": a.eval_thr,
                   "sensitivity_mode": a.sensitivity_mode, "min_identity": a.min_identity,
                   "gene_attribute": a.gene_attribute,
                   "query_proteome": os.path.abspath(a.faa),
                   "human_proteome": os.path.abspath(human_src),
                   "query_proteins": len(query_headers), "human_proteins": len(human_headers),
                   "pairs": len(rows), "query_genes": n_genes,
                   "genes_with_ortholog": len(best), "usable": bool(best)}, fh, indent=1)

    if not a.keep_work:
        for pat in ("_blast_db", "_diamond_db", "_calculation"):
            shutil.rmtree(os.path.join(work_dir, pat), ignore_errors=True)

    log(f"\n{len(rows)} ortholog pair(s); {len(best)} query gene identifier(s) resolved to a human symbol.")
    if not best:
        log("No query gene could be linked to a human symbol: check that the annotation carries "
            "protein_id attributes matching the FASTA identifiers, and that the human proteome "
            "headers carry gene symbols.")
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
