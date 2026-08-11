#!/usr/bin/env python3
"""Gather NCBI Gene + PubMed context for the top DEGs of each comparison.

Reads the base DE tables (DGE_analysis_comp*.txt), takes the strongest genes in
each direction, resolves each to an NCBI GeneID, and saves the curated RefSeq
summary plus the most recent linked papers.

Exit codes:
    0   at least one comparison produced usable context
    3   ran cleanly but nothing usable was gathered (no significant genes, no
        gene resolved, or no literature linked)
    4   NCBI E-utilities unreachable
"""

import argparse, json, os, re, sys, time, urllib.error, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "reanalyzerGSE"
ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d{6,}", re.IGNORECASE)
RETRY_BACKOFFS = [5, 15]
RETRYABLE = {429, 500, 502, 503, 504}


def log(msg):
    print(msg, file=sys.stderr, flush=True)


class Eutils:
    """Minimal E-utilities client: throttled, retried, stdlib only."""

    def __init__(self, email=None, api_key=None):
        self.email = email or os.environ.get("NCBI_EMAIL", "")
        self.api_key = api_key or os.environ.get("NCBI_API_KEY", "")
        # NCBI allows 3 req/s without a key, 10 with one
        self.min_interval = 0.11 if self.api_key else 0.35
        self._last = 0.0
        self.calls = 0

    def _wait(self):
        gap = time.time() - self._last
        if gap < self.min_interval:
            time.sleep(self.min_interval - gap)
        self._last = time.time()

    def get(self, endpoint, params, parse="json", timeout=30):
        p = dict(params)
        p["tool"] = TOOL
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(p)}"
        for attempt in range(1 + len(RETRY_BACKOFFS)):
            self._wait()
            self.calls += 1
            try:
                with urllib.request.urlopen(url, timeout=timeout) as r:
                    raw = r.read()
                return json.loads(raw) if parse == "json" else ET.fromstring(raw)
            except urllib.error.HTTPError as e:
                if e.code not in RETRYABLE or attempt == len(RETRY_BACKOFFS):
                    log(f"  [warn] {endpoint} HTTP {e.code}")
                    return None
            except Exception as e:
                if attempt == len(RETRY_BACKOFFS):
                    log(f"  [warn] {endpoint} failed: {e}")
                    return None
            time.sleep(RETRY_BACKOFFS[attempt])
        return None

    def reachable(self):
        return self.get("esearch.fcgi",
                        {"db": "pubmed", "term": "rna-seq", "retmax": 1, "retmode": "json"}) is not None


def resolve_gene(api, symbol, organism):
    d = api.get("esearch.fcgi", {"db": "gene", "retmode": "json",
                                 "term": f'{symbol}[sym] AND "{organism}"[orgn]'})
    ids = ((d or {}).get("esearchresult") or {}).get("idlist") or []
    return ids[0] if ids else None


def gene_summary(api, gid):
    d = api.get("esummary.fcgi", {"db": "gene", "id": gid, "retmode": "json"})
    rec = ((d or {}).get("result") or {}).get(gid) or {}
    return {"name": rec.get("name", ""), "description": rec.get("description", ""),
            "summary": (rec.get("summary") or "").strip()}


def linked_pmids(api, gid, n):
    d = api.get("elink.fcgi", {"dbfrom": "gene", "db": "pubmed", "id": gid, "retmode": "json"})
    try:
        links = d["linksets"][0]["linksetdbs"][0]["links"]
    except (KeyError, IndexError, TypeError):
        return [], 0
    ordered = sorted(links, key=int, reverse=True)
    return ordered[:n], len(ordered)


def fetch_papers(api, pmids):
    if not pmids:
        return []
    root = api.get("efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids),
                                   "retmode": "xml", "rettype": "abstract"}, parse="xml")
    if root is None:
        return []
    out = []
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//MedlineCitation/PMID", "")
        # every AbstractText element, so structured abstracts survive whole
        parts = []
        for node in art.findall(".//Abstract/AbstractText"):
            label = node.get("Label")
            txt = "".join(node.itertext()).strip()
            if txt:
                parts.append(f"{label}: {txt}" if label else txt)
        out.append({
            "pmid": pmid,
            "year": art.findtext(".//PubDate/Year") or art.findtext(".//PubDate/MedlineDate") or "",
            "journal": art.findtext(".//Journal/ISOAbbreviation") or art.findtext(".//Journal/Title") or "",
            "title": "".join(art.find(".//ArticleTitle").itertext()).strip()
                     if art.find(".//ArticleTitle") is not None else "",
            "abstract": " ".join(parts),
        })
    order = {p: i for i, p in enumerate(pmids)}
    out.sort(key=lambda r: order.get(r["pmid"], 999))
    return out


def read_table(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]
    return header, rows


def pick_genes(path, per_sense, annot_path=None):
    """Top `per_sense` genes each way: FDR < 0.05, ranked by FDR ascending."""
    header, rows = read_table(path)
    lfc_col = next((c for c in header if c.startswith("logFC")), None)
    fdr_col = next((c for c in header if c.lower() in ("fdr", "padj", "adj.p.val", "qvalue")), None)
    if not lfc_col or not fdr_col or not rows:
        return [], "", 0

    symbol_of = {}
    if annot_path and os.path.isfile(annot_path) and ENSEMBL_RE.match(rows[0].get("Gene_ID", "")):
        ah, arows = read_table(annot_path)
        if "SYMBOL" in ah:
            for r in arows:
                sym = (r.get("SYMBOL") or "").strip()
                if sym and sym.upper() not in ("NA", ""):
                    symbol_of[r.get("Gene_ID", "")] = sym.split(",")[0].strip()

    parsed = []
    for r in rows:
        try:
            fdr, lfc = float(r[fdr_col]), float(r[lfc_col])
        except (ValueError, KeyError, TypeError):
            continue
        parsed.append((fdr, lfc, r.get("Gene_ID", "")))
    sig = [p for p in parsed if p[0] < 0.05]
    picked = []
    for sense, keep in (("UP", lambda x: x[1] > 0), ("DOWN", lambda x: x[1] < 0)):
        for fdr, lfc, gid in sorted([p for p in sig if keep(p)], key=lambda x: x[0])[:per_sense]:
            picked.append({"gene_id": gid, "symbol": symbol_of.get(gid, gid),
                           "sense": sense, "logfc": lfc, "fdr": fdr})
    return picked, lfc_col[len("logFC"):], len(sig)


def cache_path(cache_dir, organism, symbol):
    key = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{organism}__{symbol}")
    return os.path.join(cache_dir, key + ".json")


def gather_gene(api, symbol, organism, n_papers, cache_dir):
    cf = cache_path(cache_dir, organism, symbol) if cache_dir else None
    if cf and os.path.isfile(cf):
        try:
            with open(cf, encoding="utf-8") as fh:
                rec = json.load(fh)
            if len(rec.get("papers", [])) >= n_papers or rec.get("geneid") is None:
                rec["papers"] = rec.get("papers", [])[:n_papers]
                rec["cached"] = True
                return rec
        except (OSError, ValueError):
            pass

    gid = resolve_gene(api, symbol, organism)
    rec = {"symbol": symbol, "organism": organism, "geneid": gid,
           "gene": {}, "papers": [], "n_linked": 0, "cached": False}
    if gid:
        rec["gene"] = gene_summary(api, gid)
        pmids, total = linked_pmids(api, gid, n_papers)
        rec["n_linked"] = total
        rec["papers"] = fetch_papers(api, pmids)
    if cf:
        try:
            os.makedirs(cache_dir, exist_ok=True)
            with open(cf, "w", encoding="utf-8") as fh:
                json.dump(rec, fh)
        except OSError:
            pass
    return rec


def write_comparison(out_dir, comp, contrast, records, n_sig):
    cdir = os.path.join(out_dir, comp)
    os.makedirs(cdir, exist_ok=True)

    with open(os.path.join(cdir, "genes.tsv"), "w", encoding="utf-8") as fh:
        fh.write("Gene\tSense\tlogFC\tFDR\tGeneID\tPapers_shown\tPapers_linked\tPubMed_URL\n")
        for r in records:
            gid = r["rec"].get("geneid") or ""
            url = f"https://www.ncbi.nlm.nih.gov/gene/{gid}" if gid else ""
            fh.write("\t".join([r["symbol"], r["sense"], f"{r['logfc']:+.3f}", f"{r['fdr']:.3e}",
                                gid or "-", str(len(r["rec"]["papers"])),
                                str(r["rec"]["n_linked"]), url]) + "\n")

    for r in records:
        if not r["rec"]["papers"]:
            continue
        with open(os.path.join(cdir, f"{r['symbol']}.abstracts.tsv"), "w", encoding="utf-8") as fh:
            fh.write("PMID\tYear\tJournal\tTitle\tAbstract\n")
            for p in r["rec"]["papers"]:
                fh.write("\t".join(x.replace("\t", " ").replace("\n", " ") for x in
                                   (p["pmid"], p["year"], p["journal"], p["title"], p["abstract"])) + "\n")

    blob, pmids = [], []
    blob.append(f"# Comparison {comp}" + (f" ({contrast})" if contrast else ""))
    blob.append(f"# {n_sig} genes at FDR < 0.05; showing the strongest {len(records)} in both directions.\n")
    for r in records:
        rec = r["rec"]
        head = f"## {r['symbol']}  [{r['sense']}, logFC {r['logfc']:+.2f}, FDR {r['fdr']:.2e}]"
        if rec.get("geneid"):
            head += f"  (NCBI GeneID {rec['geneid']})"
        blob.append(head)
        g = rec.get("gene") or {}
        if g.get("description"):
            blob.append(f"Description: {g['description']}")
        if g.get("summary"):
            blob.append(f"Curated RefSeq summary: {g['summary']}")
        if rec["papers"]:
            blob.append(f"Linked papers ({rec['n_linked']} curated in total; {len(rec['papers'])} most recent shown). "
                        "Cite only these PMIDs:")
            for p in rec["papers"]:
                pmids.append(p["pmid"])
                blob.append(f"- PMID:{p['pmid']} ({p['year']}, {p['journal']}) {p['title']}")
                if p["abstract"]:
                    blob.append(f"  {p['abstract']}")
        elif not rec.get("geneid"):
            blob.append("No NCBI Gene record for this identifier.")
        else:
            blob.append("No curated literature linked to this gene.")
        blob.append("")

    with open(os.path.join(cdir, "for_llm.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(blob))
    with open(os.path.join(cdir, "pmids.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(sorted(set(pmids))) + ("\n" if pmids else ""))

    resolved = sum(1 for r in records if r["rec"].get("geneid"))
    with_lit = sum(1 for r in records if r["rec"]["papers"])
    return {"comparison": comp, "contrast": contrast, "n_significant": n_sig,
            "genes": len(records), "resolved": resolved, "with_literature": with_lit,
            "papers": len(pmids), "usable": bool(with_lit or
                                                 any((r["rec"].get("gene") or {}).get("summary") for r in records))}


def main():
    p = argparse.ArgumentParser(description="Gather NCBI Gene + PubMed context for the top DEGs.")
    p.add_argument("--dge-dir", required=True, help="directory holding DGE_analysis_comp*.txt")
    p.add_argument("--organism", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--cache-dir", default=None)
    p.add_argument("--genes-per-sense", type=int, default=5)
    p.add_argument("--papers-per-gene", type=int, default=5)
    p.add_argument("--email", default=None)
    a = p.parse_args()

    organism = a.organism.replace("_", " ").strip()
    tables = sorted(f for f in os.listdir(a.dge_dir)
                    if re.fullmatch(r"DGE_analysis_comp\d+\.txt", f))
    if not tables:
        log("No DGE_analysis_comp*.txt tables found; nothing to gather.")
        return 3

    api = Eutils(email=a.email)
    if not api.reachable():
        log("NCBI E-utilities unreachable; skipping literature gathering.")
        return 4

    os.makedirs(a.out_dir, exist_ok=True)
    statuses = []
    for tbl in tables:
        comp = tbl[:-4]
        path = os.path.join(a.dge_dir, tbl)
        picked, contrast, n_sig = pick_genes(path, a.genes_per_sense,
                                             os.path.join(a.dge_dir, comp + "_annotation.txt"))
        if not picked:
            log(f"{comp}: no genes at FDR < 0.05; skipped.")
            statuses.append({"comparison": comp, "contrast": contrast, "n_significant": n_sig,
                             "genes": 0, "resolved": 0, "with_literature": 0,
                             "papers": 0, "usable": False})
            continue
        log(f"{comp}: {n_sig} significant; looking up {len(picked)} genes ...")
        records = []
        for g in picked:
            rec = gather_gene(api, g["symbol"], organism, a.papers_per_gene, a.cache_dir)
            records.append({**g, "rec": rec})
            mark = "cached" if rec.get("cached") else f"{len(rec['papers'])} paper(s)"
            log(f"  {g['symbol']:<20} {g['sense']:<5} GeneID={rec.get('geneid') or '-':<10} {mark}")
        st = write_comparison(a.out_dir, comp, contrast, records, n_sig)
        statuses.append(st)

    unresolved = []
    for st in statuses:
        cdir = os.path.join(a.out_dir, st["comparison"])
        gt = os.path.join(cdir, "genes.tsv")
        if os.path.isfile(gt):
            with open(gt, encoding="utf-8") as fh:
                next(fh, None)
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) > 4 and f[4] == "-":
                        unresolved.append(f"{st['comparison']}\t{f[0]}")
    if unresolved:
        with open(os.path.join(a.out_dir, "unresolved.tsv"), "w", encoding="utf-8") as fh:
            fh.write("Comparison\tGene\n" + "\n".join(unresolved) + "\n")

    usable = [s for s in statuses if s["usable"]]
    with open(os.path.join(a.out_dir, "status.json"), "w", encoding="utf-8") as fh:
        json.dump({"organism": organism, "api_calls": api.calls,
                   "genes_per_sense": a.genes_per_sense, "papers_per_gene": a.papers_per_gene,
                   "comparisons": statuses, "usable_comparisons": [s["comparison"] for s in usable]},
                  fh, indent=1)

    log(f"\n{len(usable)}/{len(statuses)} comparison(s) produced usable context "
        f"({api.calls} E-utilities calls).")
    return 0 if usable else 3


if __name__ == "__main__":
    sys.exit(main())
