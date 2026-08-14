#!/usr/bin/env python3
"""Parse functional annotations (GAF, GMT, GTF, GFF, or 2-column TSV/TXT)
and resolve protein/transcript/symbol identifiers to reference gene IDs.

Extracts:
  - Gene-to-GO associations (source_id <tab> Computed_GO_Process_IDs)
  - Gene-to-KEGG associations (source_id <tab> KEGG_KO_ID)
"""

import argparse
import gzip
import os
import re
import sys


def open_file(path):
    """Open plain text or gzip-compressed files transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def strip_version(val):
    """Strip version suffix from identifiers (e.g., WP_002386188.1 -> WP_002386188)."""
    if not val:
        return val
    if "." in val and not val.startswith("."):
        return val.rsplit(".", 1)[0]
    return val


def parse_reference_annotations(ref_path):
    """Build protein_id/transcript_id/locus/symbol -> gene_id mapping from reference GTF/GFF."""
    prot_to_gene = {}
    known_genes = set()
    if not ref_path or not os.path.isfile(ref_path):
        return prot_to_gene, known_genes

    with open_file(ref_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            gid = pid = locus = name = tid = ""
            for item in attrs.split(";"):
                item = item.strip()
                if not item:
                    continue
                if item.startswith("gene_id"):
                    gid = item.split('"')[1] if '"' in item else item.split("=")[1].strip('"') if "=" in item else item.split()[1]
                elif item.startswith("protein_id"):
                    pid = item.split('"')[1] if '"' in item else item.split("=")[1].strip('"') if "=" in item else item.split()[1]
                elif item.startswith("transcript_id"):
                    tid = item.split('"')[1] if '"' in item else item.split("=")[1].strip('"') if "=" in item else item.split()[1]
                elif item.startswith("locus_tag"):
                    locus = item.split('"')[1] if '"' in item else item.split("=")[1].strip('"') if "=" in item else item.split()[1]
                elif item.startswith("gene ") or item.startswith("gene="):
                    name = item.split('"')[1] if '"' in item else item.split("=")[1].strip('"') if "=" in item else item.split()[1]
                elif item.startswith("ID=") and not gid:
                    gid = item.split("=")[1].strip('"')
                elif item.startswith("Parent=") and not gid:
                    gid = item.split("=")[1].strip('"')

            gene_key = gid or locus or name
            if gene_key:
                known_genes.add(gene_key)
                known_genes.add(gene_key.upper())
                for id_val in [pid, tid, locus, name]:
                    if id_val and id_val != gene_key:
                        prot_to_gene[id_val] = gene_key
                        prot_to_gene[id_val.upper()] = gene_key
                        base_val = strip_version(id_val)
                        if base_val != id_val:
                            prot_to_gene[base_val] = gene_key
                            prot_to_gene[base_val.upper()] = gene_key

    return prot_to_gene, known_genes


def parse_functional_file(input_path, ref_path, out_go, out_kegg):
    """Extract GO and KEGG terms from input annotation file and resolve IDs."""
    prot_to_gene, known_genes = parse_reference_annotations(ref_path)

    go_pairs = set()
    kegg_pairs = set()

    is_gaf = False
    is_gmt = False
    is_gtf_gff = False

    fname = os.path.basename(input_path).lower()
    if fname.endswith(".gmt") or fname.endswith(".gmt.gz"):
        is_gmt = True
    elif fname.endswith(".gtf") or fname.endswith(".gtf.gz") or fname.endswith(".gff") or fname.endswith(".gff.gz") or fname.endswith(".gff3") or fname.endswith(".gff3.gz"):
        is_gtf_gff = True
    else:
        with open_file(input_path) as f:
            for line in f:
                line_s = line.strip()
                if not line_s or line_s.startswith("#"):
                    continue
                if line_s.startswith("!gaf-version") or line_s.startswith("!#DB"):
                    is_gaf = True
                    break
                if line_s.startswith("!"):
                    continue
                parts = line_s.split("\t")
                if len(parts) >= 15:
                    is_gaf = True
                elif len(parts) <= 2:
                    is_gaf = False  # 2-col TSV
                elif len(parts) >= 9 and ("gene_id" in line_s or "ID=" in line_s or "Parent=" in line_s):
                    is_gtf_gff = True
                else:
                    is_gaf = True
                break

    mapped_prot = 0
    mapped_sym = 0
    mapped_direct = 0

    with open_file(input_path) as f:
        for line in f:
            line_s = line.strip()
            if not line_s or line_s.startswith("#") or line_s.startswith("!"):
                continue
            parts = line_s.split("\t")

            if is_gaf:
                if len(parts) < 5:
                    continue
                col2_id = parts[1].strip()
                col3_sym = parts[2].strip() if len(parts) > 2 else ""
                goid = parts[4].strip() if len(parts) > 4 else ""

                target_gid = ""
                if col2_id in prot_to_gene:
                    target_gid = prot_to_gene[col2_id]
                    mapped_prot += 1
                elif col2_id.upper() in prot_to_gene:
                    target_gid = prot_to_gene[col2_id.upper()]
                    mapped_prot += 1
                elif strip_version(col2_id) in prot_to_gene:
                    target_gid = prot_to_gene[strip_version(col2_id)]
                    mapped_prot += 1
                elif col3_sym and col3_sym in known_genes:
                    target_gid = col3_sym
                    mapped_sym += 1
                elif col3_sym and col3_sym.upper() in known_genes:
                    target_gid = col3_sym
                    mapped_sym += 1
                elif col2_id in known_genes:
                    target_gid = col2_id
                    mapped_direct += 1
                elif col3_sym:
                    target_gid = col3_sym
                elif col2_id:
                    target_gid = col2_id

                if target_gid and goid.startswith("GO:"):
                    go_pairs.add((target_gid, goid))

                if target_gid:
                    for item in parts:
                        for km in re.findall(r'(?:ko:)?(K\d{5})', item):
                            kegg_pairs.add((target_gid, km))
                        for pm in re.findall(r'(?:ko|map)(\d{5})', item):
                            kegg_pairs.add((target_gid, f"map{pm}"))

            elif is_gtf_gff:
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                gid = ""
                m = re.search(r'\bgene_id\s+"([^"]+)"', attrs)
                if m:
                    gid = m.group(1)
                if not gid:
                    m = re.search(r'(?:^|[; ])gene_id=([^;]+)', attrs)
                    if m:
                        gid = m.group(1).strip('"')
                if not gid:
                    m = re.search(r'(?:^|[; ])ID=([^;]+)', attrs)
                    if m:
                        gid = m.group(1).strip('"')
                if not gid:
                    m = re.search(r'(?:^|[; ])Parent=([^;]+)', attrs)
                    if m:
                        gid = m.group(1).strip('"')
                if not gid:
                    m = re.search(r'\blocus_tag\s+"([^"]+)"', attrs)
                    if m:
                        gid = m.group(1)
                if not gid:
                    continue

                for go in re.findall(r'GO:\d{7}', attrs):
                    go_pairs.add((gid, go))
                for p in re.findall(r'\|(\d{7})\|', attrs):
                    go_pairs.add((gid, f"GO:{p}"))

                for km in re.findall(r'(?:ko:)?(K\d{5})', attrs):
                    kegg_pairs.add((gid, km))
                for pm in re.findall(r'(?:ko|map)(\d{5})', attrs):
                    kegg_pairs.add((gid, f"map{pm}"))

            elif is_gmt:
                if len(parts) < 3:
                    continue
                term = parts[0].strip()
                for i in range(2, len(parts)):
                    gene = parts[i].strip()
                    if not gene:
                        continue
                    target_g = prot_to_gene.get(gene, prot_to_gene.get(gene.upper(), gene))
                    if term.startswith("GO:"):
                        go_pairs.add((target_g, term))
                    elif re.match(r'^(?:ko:|map)?K?\d{5}$', term):
                        k_clean = re.sub(r'^ko:', '', term)
                        kegg_pairs.add((target_g, k_clean))
                    else:
                        go_pairs.add((target_g, term))

            else:  # 2-column TSV
                col1 = parts[0].strip()
                col2 = parts[1].strip() if len(parts) > 1 else ""
                if not col1:
                    continue
                if col1.startswith("GO:") or re.match(r'^(?:ko|map)\d{5}$', col1) or re.match(r'^(?:ko:)?K\d{5}$', col1):
                    term_str, gene_str = col1, col2
                else:
                    gene_str, term_str = col1, col2

                target_g = prot_to_gene.get(gene_str, prot_to_gene.get(gene_str.upper(), gene_str))
                for term in re.split(r'[,; ]+', term_str):
                    term = term.strip()
                    if not term:
                        continue
                    if term.startswith("GO:"):
                        go_pairs.add((target_g, term))
                    elif re.match(r'^(?:ko:|map)?(?:K\d{5}|\d{5})$', term):
                        k_clean = re.sub(r'^ko:', '', term)
                        kegg_pairs.add((target_g, k_clean))
                    else:
                        go_pairs.add((target_g, term))

    if out_go:
        os.makedirs(os.path.dirname(os.path.abspath(out_go)), exist_ok=True)
        with open(out_go, "w", encoding="utf-8") as f:
            f.write("source_id\tComputed_GO_Process_IDs\n")
            for gid, go in sorted(go_pairs):
                f.write(f"{gid}\t{go}\n")

    if out_kegg:
        if kegg_pairs:
            os.makedirs(os.path.dirname(os.path.abspath(out_kegg)), exist_ok=True)
            with open(out_kegg, "w", encoding="utf-8") as f:
                f.write("source_id\tKEGG_KO_ID\n")
                for gid, kid in sorted(kegg_pairs):
                    f.write(f"{gid}\t{kid}\n")
        elif os.path.exists(out_kegg):
            os.remove(out_kegg)

    print(f"Extracted {len(go_pairs)} gene-GO associations and {len(kegg_pairs)} gene-KEGG associations.")
    if is_gaf and mapped_prot > 0:
        print(f"  (Resolved {mapped_prot} GAF protein accessions -> reference gene IDs)")


def main():
    parser = argparse.ArgumentParser(description="Parse functional annotations and map IDs to reference GTF/GFF.")
    parser.add_argument("-i", "--input", required=True, help="Input annotation file (GAF, GMT, GTF, GFF, TSV).")
    parser.add_argument("-r", "--reference-gtf", default=None, help="Reference GTF/GFF for protein/transcript->gene mapping.")
    parser.add_argument("-o", "--output-go", required=True, help="Output path for 2-column gene->GO mapping file.")
    parser.add_argument("-k", "--output-kegg", default=None, help="Output path for 2-column gene->KEGG mapping file.")
    args = parser.parse_args()

    parse_functional_file(args.input, args.reference_gtf, args.output_go, args.output_kegg)


if __name__ == "__main__":
    main()
