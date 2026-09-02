#!/usr/bin/env python3
"""Look the DEGs of each comparison up in ENCODE, and report their regulatory landscape.

For every comparison in --dge-dir the strongest DEGs are taken, translated to a human gene
symbol (directly when the organism is human, otherwise through the human orthologs detected
by ortho_human.py), and queried against ENCODE:

  * ENCODE-rE2G element-gene links   -> the candidate regulatory elements of each gene
                                        (promoter/enhancer class, distance to TSS, score)
  * ENCODE TF ChIP-seq               -> the transcription factors with a peak in those elements
  * other ENCODE assays (optional)   -> DNase/ATAC/histone-mark experiments over the same regions

Each transcription factor is then checked against the differential expression table of the
same comparison, so that a TF whose own encoding gene (or its ortholog) is itself a DEG is
flagged in the TF_is_DEG column.

Outputs in --out-dir:
    status.json                 parameters, per-comparison counts, usable comparisons
    biosamples.tsv              the ENCODE rE2G annotations used
    unmapped.tsv                DEGs with no usable human symbol
    <comparison>/elements.tsv   regulatory elements per DEG
    <comparison>/tfs.tsv        TF ChIP-seq peaks per element
    <comparison>/tf_summary.tsv one row per TF, with its DEG status
    <comparison>/assays.tsv     other ENCODE assays over the same regions (when enabled)

Exit codes:
    0   at least one comparison produced usable results
    3   ran cleanly but nothing usable (no DEGs, none mapped, or no ENCODE elements)
    4   ENCODE unreachable, or no rE2G annotation could be resolved
"""

import argparse, collections, csv, gzip, http.client, json, os, re, sys, time
import urllib.error, urllib.parse, urllib.request

ENCODE = "https://www.encodeproject.org"
UA = {"User-Agent": "reanalyzerGSE-encode/1.0", "Accept": "application/json"}
RE2G_TYPE = "element gene regulatory interaction predictions"
LINK_OUTPUT_TYPE = "thresholded element gene links"
ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d{6,}", re.IGNORECASE)
TRANSIENT = (OSError, http.client.HTTPException, json.JSONDecodeError)
RETRY_PAUSES = [3.0, 8.0, 20.0]

ELEMENT_COLS = ["gene_id", "human_symbol", "sense", "logFC", "FDR", "element", "chrom", "start",
                "end", "class", "isSelfPromoter", "distanceToTSS", "rE2G_score", "cell_type",
                "re2g_accession"]
TF_COLS = ["gene_id", "human_symbol", "element", "class", "rE2G_score", "distanceToTSS", "TF",
           "n_biosamples", "biosamples", "n_elements_with_this_TF"]
TF_SUMMARY_COLS = ["TF", "n_target_genes", "n_at_promoter", "n_elements", "target_genes",
                   "TF_is_DEG", "TF_gene_id", "TF_query_gene", "TF_logFC", "TF_FDR",
                   "TF_in_DE_table", "TF_ortholog_identity"]
ASSAY_COLS = ["element", "chrom", "start", "end", "assay_title", "target", "n_experiments",
              "biosamples"]


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def get_json(url, retries=len(RETRY_PAUSES) + 1):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as fh:
                return json.loads(fh.read().decode("utf-8", "replace"))
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt == retries - 1:
                raise
            log(f"  retry {attempt + 1}/{retries - 1} after HTTP {exc.code}")
            time.sleep(RETRY_PAUSES[min(attempt, len(RETRY_PAUSES) - 1)])
        except TRANSIENT as exc:
            if attempt == retries - 1:
                raise
            log(f"  retry {attempt + 1}/{retries - 1} after {type(exc).__name__}: {exc}")
            time.sleep(RETRY_PAUSES[min(attempt, len(RETRY_PAUSES) - 1)])
    return None


def embedded(obj, key):
    v = obj.get(key)
    return v if isinstance(v, dict) else {}


def as_list(v):
    if isinstance(v, str):
        return [v]
    return list(v) if isinstance(v, (list, tuple)) else []


def write_tsv(path, cols, rows):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log(f"  wrote {path}  ({len(rows)} rows)")


######### DE tables

def read_table(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]
    return header, rows


def symbol_map(annot_path, gene_ids):
    """Gene_ID -> SYMBOL from the per-comparison annotation table, for ENSEMBL-style ids."""
    out = {}
    if not annot_path or not os.path.isfile(annot_path):
        return out
    if not any(ENSEMBL_RE.match(g or "") for g in list(gene_ids)[:20]):
        return out
    header, rows = read_table(annot_path)
    if "SYMBOL" not in header:
        return out
    for r in rows:
        sym = (r.get("SYMBOL") or "").strip()
        if sym and sym.upper() != "NA":
            out[r.get("Gene_ID", "")] = sym.split(",")[0].strip()
    return out


def parse_de(path, fdr_cut):
    """(all rows keyed by gene id, contrast label, significant rows) from one DE table."""
    header, rows = read_table(path)
    lfc_col = next((c for c in header if c.startswith("logFC")), None)
    fdr_col = next((c for c in header if c.lower() in ("fdr", "padj", "adj.p.val", "qvalue")), None)
    if not lfc_col or not fdr_col or not rows:
        return {}, "", []
    by_gene, parsed = {}, []
    for r in rows:
        gid = r.get("Gene_ID", "")
        if not gid:
            continue
        try:
            fdr, lfc = float(r[fdr_col]), float(r[lfc_col])
        except (KeyError, TypeError, ValueError):
            continue
        rec = {"gene_id": gid, "logfc": lfc, "fdr": fdr}
        prev = by_gene.get(gid)
        if prev is None or fdr < prev["fdr"]:
            by_gene[gid] = rec
        parsed.append(rec)
    sig = [r for r in parsed if r["fdr"] < fdr_cut]
    return by_gene, lfc_col[len("logFC"):], sig


def pick_degs(sig, per_sense):
    picked = []
    for sense, keep in (("UP", lambda x: x["logfc"] > 0), ("DOWN", lambda x: x["logfc"] < 0)):
        chosen = sorted([r for r in sig if keep(r)], key=lambda x: x["fdr"])[:per_sense]
        for r in chosen:
            picked.append({**r, "sense": sense})
    return picked


######### human symbol resolution

def load_orthologs(path):
    """query gene -> human symbol, and human symbol -> [query genes], from gene_to_human.tsv."""
    fwd, rev, ident = {}, collections.defaultdict(list), {}
    if not path or not os.path.isfile(path):
        return fwd, rev, ident
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        q, h = (r.get("query_gene") or "").strip(), (r.get("human_symbol") or "").strip()
        if not q or not h:
            continue
        fwd[q] = h
        fwd.setdefault(q.upper(), h)
        if q not in rev[h.upper()]:
            rev[h.upper()].append(q)
        ident[q] = (r.get("perc_identity") or "").strip()
    log(f"  orthologs: {len(fwd)} query gene identifier(s), {len(rev)} human symbol(s)")
    return fwd, rev, ident


def resolve_symbols(picked, is_human, symbols, fwd):
    """Attach a human_symbol to each picked DEG; return the resolved ones and the failures."""
    ok, bad = [], []
    for g in picked:
        gid = g["gene_id"]
        if is_human:
            sym = symbols.get(gid) or gid
        else:
            sym = fwd.get(gid) or fwd.get(gid.upper()) or ""
        sym = (sym or "").strip()
        if sym and not ENSEMBL_RE.match(sym):
            ok.append({**g, "human_symbol": sym})
        else:
            bad.append(g)
    return ok, bad


######### ENCODE rE2G annotations

def search_re2g(biosample=None, limit=200):
    params = [("type", "Annotation"), ("annotation_type", RE2G_TYPE), ("format", "json"),
              ("limit", str(limit)), ("field", "accession"),
              ("field", "biosample_ontology.term_name"), ("field", "description"),
              ("field", "assembly")]
    if biosample:
        params.append(("biosample_ontology.term_name", biosample))
    url = f"{ENCODE}/search/?" + urllib.parse.urlencode(params)
    d = get_json(url) or {}
    out = []
    for x in d.get("@graph", []):
        out.append({"accession": x.get("accession", ""),
                    "biosample": embedded(x, "biosample_ontology").get("term_name") or "",
                    "assembly": ",".join(as_list(x.get("assembly"))),
                    "description": x.get("description", "")})
    return out, d.get("total", len(out))


def pick_link_file(accession, prefer=None):
    ann = get_json(f"{ENCODE}/annotations/{accession}/?format=json")
    if not ann:
        log(f"  [warn] annotation {accession} could not be read")
        return None, None
    cands = [f for f in ann.get("files", []) if f.get("output_type") == LINK_OUTPUT_TYPE]
    if not cands:
        log(f"  [warn] {accession} has no '{LINK_OUTPUT_TYPE}' file")
        return None, None
    chosen = None
    if prefer:
        chosen = next((f for f in cands if f.get("accession") == prefer), None)
        if not chosen:
            log(f"  [warn] {prefer} is not a link file of {accession}; falling back to the smallest")
    if not chosen:
        chosen = min(cands, key=lambda f: f.get("file_size") or 0)
    return ann.get("description", ""), chosen


def fetch_link_file(meta, cache_dir):
    local = os.path.join(cache_dir, f"{meta['accession']}.bed.gz")
    if os.path.isfile(local) and os.path.getsize(local) > 1000:
        return local
    os.makedirs(cache_dir, exist_ok=True)
    url = ENCODE + meta["href"]
    log(f"  downloading {url}")
    tmp = local + ".part"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA["User-Agent"]})
        with urllib.request.urlopen(req, timeout=1800) as r, open(tmp, "wb") as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    except TRANSIENT as exc:
        if os.path.exists(tmp):
            os.remove(tmp)
        log(f"  [warn] download failed: {exc}")
        return None
    os.replace(tmp, local)
    return local


COLUMN_ALIASES = {
    "gene": ("TargetGene", "targetGene", "target_gene", "gene"),
    "name": ("name", "element", "ElementName"),
    "chrom": ("chr", "chrom", "#chr", "chromosome"),
    "start": ("start", "chromStart"),
    "end": ("end", "chromEnd"),
    "class": ("class", "ElementClass", "element_class"),
    "self_promoter": ("isSelfPromoter", "is_self_promoter"),
    "distance": ("distanceToTSS.Feature", "distanceToTSS", "distance", "distanceToTSS.feature"),
    "score": ("Score", "score", "rE2G.Score", "ENCODE-rE2G.Score"),
    "celltype": ("CellType", "cellType", "biosample"),
}


def resolve_columns(header):
    idx = {name: i for i, name in enumerate(header)}
    out = {}
    for key, names in COLUMN_ALIASES.items():
        for n in names:
            if n in idx:
                out[key] = idx[n]
                break
    return out


def load_elements(local, wanted, accession):
    """gene symbol -> element records, for the symbols in `wanted`."""
    out = collections.defaultdict(list)
    with gzip.open(local, "rt", errors="replace") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        col = resolve_columns(header)
        missing = [k for k in ("gene", "chrom", "start", "end") if k not in col]
        if missing:
            log(f"  [warn] {accession}: link file lacks the column(s) {missing}; skipped")
            return out
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) <= col["end"]:
                continue
            gene = c[col["gene"]]
            if gene not in wanted:
                continue
            get = lambda k, d="": c[col[k]] if k in col and col[k] < len(c) else d
            out[gene].append({
                "element": get("name", f"{c[col['chrom']]}:{c[col['start']]}-{c[col['end']]}"),
                "chrom": c[col["chrom"]], "start": c[col["start"]], "end": c[col["end"]],
                "class": get("class"), "isSelfPromoter": get("self_promoter"),
                "distanceToTSS": get("distance"), "rE2G_score": get("score", "0"),
                "cell_type": get("celltype"), "re2g_accession": accession,
            })
    return out


######### ENCODE region search

def region_experiments(region, assay_titles, page):
    """{assay_title: {target: {biosamples}}} for the ENCODE experiments overlapping `region`."""
    out = collections.defaultdict(lambda: collections.defaultdict(set))
    for assay in assay_titles:
        frm, total = 0, None
        while True:
            params = [("region", region), ("genome", "GRCh38"), ("assay_title", assay),
                      ("limit", str(page)), ("from", str(frm)), ("format", "json"),
                      ("field", "target.label"), ("field", "biosample_ontology.term_name"),
                      ("field", "assay_title")]
            d = get_json(f"{ENCODE}/region-search/?" + urllib.parse.urlencode(params))
            if d is None:
                break
            if total is None:
                total = d.get("total", 0)
            graph = d.get("@graph", [])
            for x in graph:
                title = x.get("assay_title") or assay
                if isinstance(title, (list, tuple)):
                    title = title[0] if title else assay
                label = embedded(x, "target").get("label") or title
                term = embedded(x, "biosample_ontology").get("term_name") or ""
                if term:
                    out[title][label].add(term)
                else:
                    out[title][label]
            frm += len(graph)
            if not graph or len(graph) < page or frm >= (total or 0):
                break
    return {t: {k: sorted(v) for k, v in d.items()} for t, d in out.items()}


def cached_regions(cache_path):
    if cache_path and os.path.isfile(cache_path):
        try:
            with open(cache_path, encoding="utf-8") as fh:
                return json.load(fh)
        except (OSError, ValueError):
            pass
    return {}


def save_cache(cache_path, cache):
    if not cache_path:
        return
    try:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        tmp = cache_path + ".part"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(cache, fh)
        os.replace(tmp, cache_path)
    except OSError:
        pass


######### per-comparison assembly

def float_or(v, default):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def tf_deg_status(tf, is_human, rev, ident, de_by_gene, de_ci, symbols_inv, fdr_cut):
    """Is the gene encoding `tf` itself differentially expressed in this comparison?"""
    rec = {"TF_is_DEG": "no", "TF_gene_id": "", "TF_query_gene": "", "TF_logFC": "",
           "TF_FDR": "", "TF_in_DE_table": "no", "TF_ortholog_identity": ""}
    if is_human:
        candidates = [tf] + symbols_inv.get(tf.upper(), [])
    else:
        candidates = list(rev.get(tf.upper(), []))
    hit = None
    for cand in candidates:
        d = de_by_gene.get(cand) or de_ci.get(cand.upper())
        if d and (hit is None or d["fdr"] < hit[1]["fdr"]):
            hit = (cand, d)
    if not candidates:
        return rec
    rec["TF_query_gene"] = ",".join(candidates[:5])
    rec["TF_ortholog_identity"] = ident.get(candidates[0], "")
    if hit:
        cand, d = hit
        rec.update({"TF_in_DE_table": "yes", "TF_gene_id": d["gene_id"],
                    "TF_query_gene": cand, "TF_logFC": f"{d['logfc']:+.3f}",
                    "TF_FDR": f"{d['fdr']:.3e}",
                    "TF_is_DEG": "yes" if d["fdr"] < fdr_cut else "no"})
    return rec


def main():
    p = argparse.ArgumentParser(description="Look the DEGs up in ENCODE for regulatory elements and TF binding.")
    p.add_argument("--dge-dir", help="directory holding DGE_analysis_comp*.txt")
    p.add_argument("--out-dir", help="where to write the ENCODE tables")
    p.add_argument("--organism", default="", help="scientific name of the analysed organism")
    p.add_argument("--orthologs-tsv", default="",
                   help="gene_to_human.tsv from ortho_human.py; required unless the organism is human")
    p.add_argument("--cache-dir", default="", help="reused across runs for the rE2G files and region searches")
    p.add_argument("--biosample", default="K562",
                   help="comma-separated ENCODE biosample term name(s) whose rE2G predictions to use "
                        "(default K562); see --list-biosamples")
    p.add_argument("--re2g-accession", default="",
                   help="comma-separated ENCODE annotation accession(s), used instead of --biosample")
    p.add_argument("--re2g-file", default="", help="specific ENCFF of the thresholded element-gene links")
    p.add_argument("--list-biosamples", action="store_true",
                   help="print the ENCODE rE2G annotations available and exit")
    p.add_argument("--genes-per-sense", type=int, default=25,
                   help="top N DEGs per direction, ranked by FDR (default 25)")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff for calling a gene a DEG")
    p.add_argument("--max-elements-per-gene", type=int, default=10,
                   help="keep only the N highest-scoring elements per gene (default 10)")
    p.add_argument("--max-regions", type=int, default=1500,
                   help="hard cap on distinct genomic regions queried per run (default 1500)")
    p.add_argument("--extra-assays", default="DNase-seq,ATAC-seq,Histone ChIP-seq",
                   help="comma-separated extra ENCODE assay titles to record over the same regions "
                        "(empty string to only do TF ChIP-seq)")
    p.add_argument("--pause", type=float, default=0.34, help="seconds between ENCODE region searches")
    p.add_argument("--page", type=int, default=100, help="ENCODE region-search page size")
    a = p.parse_args()

    if a.list_biosamples:
        rows, total = search_re2g(limit=500)
        log(f"{total} ENCODE rE2G annotation(s); showing {len(rows)}")
        print("accession\tbiosample\tassembly\tdescription")
        for r in sorted(rows, key=lambda x: (x["biosample"], x["accession"])):
            print("\t".join([r["accession"], r["biosample"], r["assembly"], r["description"]]))
        return 0

    if not a.dge_dir or not a.out_dir:
        log("ERROR: --dge-dir and --out-dir are required unless --list-biosamples is given.")
        return 4

    organism = a.organism.replace("_", " ").strip()
    is_human = organism.lower() in ("homo sapiens", "human")
    a.out_dir = os.path.abspath(a.out_dir)
    cache_dir = os.path.abspath(a.cache_dir) if a.cache_dir else os.path.join(a.out_dir, "cache")
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)
    extra_assays = [x.strip() for x in a.extra_assays.split(",") if x.strip()]

    tables = sorted(f for f in os.listdir(a.dge_dir)
                    if re.fullmatch(r"DGE_analysis_comp\d+\.txt", f))
    if not tables:
        log("No DGE_analysis_comp*.txt table found; nothing to look up.")
        return 3

    fwd, rev, ident = ({}, {}, {})
    if not is_human:
        if not a.orthologs_tsv or not os.path.isfile(a.orthologs_tsv):
            log("ERROR: the organism is not human, so a gene_to_human.tsv from ortho_human.py is "
                "required (pass --orthologs-tsv). Enable ortho_detection to produce it.")
            return 4
        fwd, rev, ident = load_orthologs(a.orthologs_tsv)
        if not fwd:
            log("ERROR: no usable ortholog mapping in " + a.orthologs_tsv)
            return 3

    ######### which DEGs, and their human symbols
    per_comp, unmapped = {}, []
    for tbl in tables:
        comp = tbl[:-4]
        de_by_gene, contrast, sig = parse_de(os.path.join(a.dge_dir, tbl), a.fdr)
        picked = pick_degs(sig, a.genes_per_sense)
        symbols = symbol_map(os.path.join(a.dge_dir, comp + "_annotation.txt"),
                             [g["gene_id"] for g in picked]) if is_human else {}
        symbols_inv = collections.defaultdict(list)
        for gid, sym in symbols.items():
            symbols_inv[sym.upper()].append(gid)
        resolved, bad = resolve_symbols(picked, is_human, symbols, fwd)
        for g in bad:
            unmapped.append({"comparison": comp, "gene_id": g["gene_id"], "sense": g["sense"]})
        de_ci = {}
        for gid, rec in de_by_gene.items():
            prev = de_ci.get(gid.upper())
            if prev is None or rec["fdr"] < prev["fdr"]:
                de_ci[gid.upper()] = rec
        per_comp[comp] = {"contrast": contrast, "de_by_gene": de_by_gene, "de_ci": de_ci,
                          "n_significant": len(sig), "resolved": resolved,
                          "symbols_inv": symbols_inv}
        log(f"{comp}: {len(sig)} gene(s) at FDR < {a.fdr}; {len(picked)} taken, "
            f"{len(resolved)} with a human symbol")

    wanted = sorted({g["human_symbol"] for v in per_comp.values() for g in v["resolved"]})
    if not wanted:
        log("No DEG could be translated to a human gene symbol.")
        return 3
    log(f"\n{len(wanted)} distinct human symbol(s) to look up in ENCODE")

    ######### rE2G annotations
    annotations = []
    if a.re2g_accession:
        for acc in (x.strip() for x in a.re2g_accession.split(",") if x.strip()):
            annotations.append({"accession": acc, "biosample": "", "description": "", "assembly": ""})
    else:
        for bs in (x.strip() for x in a.biosample.split(",") if x.strip()):
            found, total = search_re2g(biosample=bs, limit=50)
            if not found:
                log(f"  [warn] no rE2G annotation for biosample '{bs}'")
                continue
            log(f"  biosample '{bs}': {total} annotation(s), using {found[0]['accession']}")
            annotations.append(found[0])
    if not annotations:
        log("ERROR: no ENCODE rE2G annotation could be resolved. Run with --list-biosamples to see "
            "the biosample term names available, or give --re2g-accession directly.")
        return 4

    used, elements_by_gene = [], collections.defaultdict(list)
    want_set = set(wanted)
    for ann in annotations:
        desc, meta = pick_link_file(ann["accession"], a.re2g_file or None)
        if not meta:
            continue
        local = fetch_link_file(meta, cache_dir)
        if not local:
            continue
        got = load_elements(local, want_set, ann["accession"])
        n_el = sum(len(v) for v in got.values())
        log(f"  {ann['accession']} ({ann.get('biosample') or desc}): {len(got)} gene(s), {n_el} element(s)")
        for gene, els in got.items():
            elements_by_gene[gene].extend(els)
        used.append({"accession": ann["accession"], "biosample": ann.get("biosample", ""),
                     "file": meta.get("accession", ""), "description": desc or ann.get("description", ""),
                     "genes": len(got), "elements": n_el})
    write_tsv(os.path.join(a.out_dir, "biosamples.tsv"),
              ["accession", "biosample", "file", "description", "genes", "elements"], used)
    if not elements_by_gene:
        log("None of the requested genes has predicted regulatory elements in the chosen biosample(s).")
        return 3

    for gene, els in elements_by_gene.items():
        els.sort(key=lambda e: -float_or(e["rE2G_score"], 0.0))
        if a.max_elements_per_gene > 0:
            del els[a.max_elements_per_gene:]

    ######### TF ChIP-seq (and other assays) over those elements
    regions = []
    for gene in sorted(elements_by_gene):
        for el in elements_by_gene[gene]:
            regions.append((float_or(el["rE2G_score"], 0.0), f"{el['chrom']}:{el['start']}-{el['end']}"))
    ordered, seen = [], set()
    for _, reg in sorted(regions, key=lambda x: -x[0]):
        if reg not in seen:
            seen.add(reg)
            ordered.append(reg)
    truncated = False
    if a.max_regions > 0 and len(ordered) > a.max_regions:
        log(f"  {len(ordered)} distinct regions exceeds --max-regions {a.max_regions}; keeping the "
            "highest-scoring ones")
        ordered = ordered[:a.max_regions]
        truncated = True
    keep_regions = set(ordered)

    assay_titles = ["TF ChIP-seq"] + extra_assays
    cache_key = "regions_" + re.sub(r"[^A-Za-z0-9]+", "_", "_".join(assay_titles)) + ".json"
    cache_path = os.path.join(cache_dir, cache_key)
    cache = cached_regions(cache_path)
    log(f"  region cache: {len(cache)} region(s) already done ({cache_path})")

    failed, done = [], 0
    for n, region in enumerate(ordered, 1):
        if region in cache:
            continue
        try:
            cache[region] = region_experiments(region, assay_titles, a.page)
        except TRANSIENT as exc:
            log(f"  FAILED {region}: {type(exc).__name__} — not cached, re-run to retry")
            failed.append((region, type(exc).__name__))
            continue
        done += 1
        if done % 20 == 0:
            save_cache(cache_path, cache)
            log(f"  [{n}/{len(ordered)}] {region}: "
                f"{len(cache[region].get('TF ChIP-seq', {}))} TF(s)")
        time.sleep(a.pause)
    save_cache(cache_path, cache)
    if failed:
        write_tsv(os.path.join(a.out_dir, "failed_regions.tsv"), ["region", "error"],
                  [{"region": r, "error": e} for r, e in failed])
        log(f"INCOMPLETE: {len(failed)} region(s) could not be fetched; re-run to retry only those.")

    ######### per-comparison tables
    statuses = []
    for comp in sorted(per_comp):
        v = per_comp[comp]
        cdir = os.path.join(a.out_dir, comp)
        element_rows, tf_rows, assay_rows = [], [], []
        assay_seen = set()
        for g in sorted(v["resolved"], key=lambda x: x["fdr"]):
            sym = g["human_symbol"]
            for el in elements_by_gene.get(sym, []):
                region = f"{el['chrom']}:{el['start']}-{el['end']}"
                element_rows.append({**el, "gene_id": g["gene_id"], "human_symbol": sym,
                                     "sense": g["sense"], "logFC": f"{g['logfc']:+.3f}",
                                     "FDR": f"{g['fdr']:.3e}"})
                if region not in keep_regions or region not in cache:
                    continue
                per_assay = cache[region]
                for tf, biosamples in sorted(per_assay.get("TF ChIP-seq", {}).items(),
                                             key=lambda kv: (-len(kv[1]), kv[0])):
                    tf_rows.append({"gene_id": g["gene_id"], "human_symbol": sym,
                                    "element": el["element"], "class": el["class"],
                                    "rE2G_score": el["rE2G_score"],
                                    "distanceToTSS": el["distanceToTSS"], "TF": tf,
                                    "n_biosamples": len(biosamples),
                                    "biosamples": ",".join(biosamples)})
                for title in extra_assays:
                    for target, biosamples in sorted(per_assay.get(title, {}).items()):
                        if (el["element"], title, target) in assay_seen:
                            continue
                        assay_seen.add((el["element"], title, target))
                        assay_rows.append({"element": el["element"], "chrom": el["chrom"],
                                           "start": el["start"], "end": el["end"],
                                           "assay_title": title, "target": target,
                                           "n_experiments": len(biosamples),
                                           "biosamples": ",".join(biosamples)})

        if not element_rows:
            statuses.append({"comparison": comp, "contrast": v["contrast"],
                             "n_significant": v["n_significant"], "genes": len(v["resolved"]),
                             "genes_with_elements": 0, "elements": 0, "tfs": 0, "tf_is_deg": 0,
                             "usable": False})
            continue

        seen_tf = collections.Counter((r["human_symbol"], r["TF"]) for r in tf_rows)
        for r in tf_rows:
            r["n_elements_with_this_TF"] = seen_tf[(r["human_symbol"], r["TF"])]

        targets, at_promoter, n_elements = (collections.defaultdict(set),
                                            collections.defaultdict(set),
                                            collections.defaultdict(set))
        for r in tf_rows:
            targets[r["TF"]].add(r["human_symbol"])
            n_elements[r["TF"]].add(r["element"])
            if (r["class"] or "").lower().startswith("promoter"):
                at_promoter[r["TF"]].add(r["human_symbol"])
        summary = []
        for tf in sorted(targets):
            rec = {"TF": tf, "n_target_genes": len(targets[tf]),
                   "n_at_promoter": len(at_promoter.get(tf, ())),
                   "n_elements": len(n_elements[tf]),
                   "target_genes": ",".join(sorted(targets[tf]))}
            rec.update(tf_deg_status(tf, is_human, rev, ident, v["de_by_gene"], v["de_ci"],
                                     v["symbols_inv"], a.fdr))
            summary.append(rec)
        summary.sort(key=lambda r: (r["TF_is_DEG"] != "yes",
                                    float_or(r["TF_FDR"], 9.0), -r["n_target_genes"]))

        write_tsv(os.path.join(cdir, "elements.tsv"), ELEMENT_COLS, element_rows)
        write_tsv(os.path.join(cdir, "tfs.tsv"), TF_COLS, tf_rows)
        write_tsv(os.path.join(cdir, "tf_summary.tsv"), TF_SUMMARY_COLS, summary)
        if assay_rows:
            write_tsv(os.path.join(cdir, "assays.tsv"), ASSAY_COLS, assay_rows)

        n_deg_tf = sum(1 for r in summary if r["TF_is_DEG"] == "yes")
        statuses.append({"comparison": comp, "contrast": v["contrast"],
                         "n_significant": v["n_significant"], "genes": len(v["resolved"]),
                         "genes_with_elements": len({r["human_symbol"] for r in element_rows}),
                         "elements": len({r["element"] for r in element_rows}),
                         "tfs": len(targets), "tf_is_deg": n_deg_tf,
                         "other_assay_rows": len(assay_rows), "usable": True})
        log(f"{comp}: {len({r['human_symbol'] for r in element_rows})} gene(s) with element(s), "
            f"{len(targets)} TF(s), {n_deg_tf} of them differentially expressed themselves")

    if unmapped:
        write_tsv(os.path.join(a.out_dir, "unmapped.tsv"),
                  ["comparison", "gene_id", "sense"], unmapped)

    usable = [s for s in statuses if s["usable"]]
    with open(os.path.join(a.out_dir, "status.json"), "w", encoding="utf-8") as fh:
        json.dump({"organism": organism, "is_human": is_human,
                   "genes_per_sense": a.genes_per_sense, "fdr": a.fdr,
                   "max_elements_per_gene": a.max_elements_per_gene,
                   "extra_assays": extra_assays, "regions_queried": len(ordered),
                   "regions_truncated": truncated, "regions_failed": len(failed),
                   "biosamples": used, "comparisons": statuses,
                   "usable_comparisons": [s["comparison"] for s in usable]}, fh, indent=1)

    log(f"\n{len(usable)}/{len(statuses)} comparison(s) produced usable ENCODE results.")
    return 0 if usable else 3


if __name__ == "__main__":
    sys.exit(main())
