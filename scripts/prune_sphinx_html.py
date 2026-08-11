#!/usr/bin/env python3
"""Shrink a built Sphinx html/ tree.

Sphinx's html_extra_path copies whole result directories into html/, so the tree ends up
holding a second copy of most of the analysis, and the report links only a fraction of it.

Two passes, neither of which can lose data: a file is only ever touched when the same
file still exists elsewhere in the run folder.

  1. delete copies that no built page references
  2. replace each referenced copy with a hard link to its original

Files with no counterpart outside html/ are always kept, whether linked or not.
"""

import argparse
import os
import re
import sys
from urllib.parse import unquote

GENERATED_DIRS = {"_static", "_downloads", "_sources", "_images", "_modules", "_panels_static"}
GENERATED_NAMES = {"objects.inv", "searchindex.js", ".nojekyll", ".buildinfo"}
KEEP_EXT = {
    ".html", ".htm", ".css", ".js", ".map", ".json", ".svg", ".ico",
    ".woff", ".woff2", ".ttf", ".eot", ".otf", ".inv",
}
SCAN_EXT = {".html", ".htm", ".js", ".css"}

ATTR_RE = re.compile(
    rb"""(?:href|src|data-src|data-href|srcset|action|content)\s*=\s*["']([^"'>]+)["']""",
    re.IGNORECASE,
)
URL_RE = re.compile(rb"""url\(\s*['"]?([^'")]+)""", re.IGNORECASE)
EXT_RE = re.compile(
    rb"""\.(?:txt|tsv|csv|xlsx|xls|pdf|png|jpg|jpeg|gif|svg|zip|gz|tar|qs2|rds|bw|bed|bam|cram|bai|csi|html|htm|json|log|out|xml|parquet|npz)(?![\w])""",
    re.IGNORECASE,
)
NAME_BYTES = frozenset(
    b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-+")
MAX_NAME = 200


def collect_references(root, extra_files):
    refs = set()

    def add(raw):
        try:
            token = unquote(raw.decode("utf-8", "replace")).strip()
        except Exception:
            return
        token = token.split("#", 1)[0].split("?", 1)[0].strip()
        if not token or "://" in token or token.startswith(("mailto:", "data:", "javascript:")):
            return
        refs.add(token.lstrip("./"))
        refs.add(os.path.basename(token))

    targets = list(extra_files)
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            if os.path.splitext(fn)[1].lower() in SCAN_EXT:
                targets.append(os.path.join(dirpath, fn))

    for path in targets:
        try:
            with open(path, "rb") as fh:
                blob = fh.read()
        except OSError:
            continue
        for m in ATTR_RE.finditer(blob):
            add(m.group(1))
        for m in URL_RE.finditer(blob):
            add(m.group(1))
        for m in EXT_RE.finditer(blob):
            start = m.start()
            left = max(0, start - MAX_NAME)
            i = start
            while i > left and blob[i - 1] in NAME_BYTES:
                i -= 1
            if i < start:
                add(blob[i:m.end()])
    return refs


def is_generated(rel):
    return rel.split(os.sep)[0] in GENERATED_DIRS or os.path.basename(rel) in GENERATED_NAMES


def is_referenced(rel, name, refs):
    return name in refs or rel in refs or rel.replace(os.sep, "/") in refs


def index_originals(source_root, html_root):
    index = {}
    for dirpath, dirnames, filenames in os.walk(source_root):
        dirnames[:] = [d for d in dirnames
                       if os.path.join(dirpath, d) != html_root and d != "sphinx_report"]
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            try:
                if os.path.islink(full):
                    continue
                index.setdefault(fn, []).append((full, os.path.getsize(full)))
            except OSError:
                continue
    return index


def same_content(a, b, sample=8192):
    try:
        with open(a, "rb") as fa, open(b, "rb") as fb:
            if fa.read(sample) != fb.read(sample):
                return False
            size = os.path.getsize(a)
            if size > 2 * sample:
                fa.seek(-sample, os.SEEK_END)
                fb.seek(-sample, os.SEEK_END)
                if fa.read(sample) != fb.read(sample):
                    return False
    except OSError:
        return False
    return True


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return "%.1f %s" % (n, unit)
        n /= 1024.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("html_root", help="the built sphinx_report/html directory")
    ap.add_argument("--source-root", default=None,
                    help="run folder holding the originals (default: two levels above html_root)")
    ap.add_argument("--extra-html", action="append", default=[],
                    help="additional page to scan for references, e.g. final_report.html")
    ap.add_argument("--no-link", action="store_true",
                    help="do not hard-link surviving copies to their originals")
    ap.add_argument("--dry-run", action="store_true", help="report only, change nothing")
    ap.add_argument("--verbose", action="store_true", help="list every file acted on")
    ap.add_argument("--force", action="store_true", help="skip the sphinx-tree sanity check")
    args = ap.parse_args()

    html_root = os.path.realpath(args.html_root)
    if not os.path.isdir(html_root):
        print("Report tidying: %s is not a directory, nothing to do" % args.html_root)
        return 0
    if not args.force and (os.path.basename(html_root) != "html" or not os.path.isfile(
            os.path.join(os.path.dirname(html_root), "conf.py"))):
        print("Report tidying: %s does not look like a Sphinx build tree, skipping" % html_root)
        return 0

    source_root = os.path.realpath(
        args.source_root or os.path.dirname(os.path.dirname(html_root)))

    refs = collect_references(html_root, [p for p in args.extra_html if os.path.isfile(p)])
    originals = index_originals(source_root, html_root)

    removed = removed_bytes = linked = linked_bytes = unique = 0
    for dirpath, _, filenames in os.walk(html_root):
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, html_root)
            try:
                st = os.lstat(full)
            except OSError:
                continue
            if not os.path.isfile(full) or os.path.islink(full):
                continue

            twin = None
            for cand, size in originals.get(fn, ()):
                if size == st.st_size and os.path.realpath(cand) != full and same_content(full, cand):
                    twin = cand
                    break
            if twin is None:
                unique += 1
                continue

            keep = (is_generated(rel) or os.path.splitext(fn)[1].lower() in KEEP_EXT
                    or is_referenced(rel, fn, refs))
            if not keep:
                removed += 1
                removed_bytes += st.st_size
                if args.verbose:
                    print("  remove %s" % rel)
                if not args.dry_run:
                    try:
                        os.remove(full)
                    except OSError as exc:
                        print("Report tidying: could not remove %s (%s)" % (rel, exc))
                        removed -= 1
                        removed_bytes -= st.st_size
                continue

            if args.no_link or st.st_nlink > 1 or st.st_size == 0:
                continue
            linked += 1
            linked_bytes += st.st_size
            if args.verbose:
                print("  link   %s" % rel)
            if not args.dry_run:
                tmp = full + ".rgse_relink"
                try:
                    os.link(twin, tmp)
                    os.replace(tmp, full)
                except OSError:
                    if os.path.lexists(tmp):
                        try:
                            os.remove(tmp)
                        except OSError:
                            pass
                    linked -= 1
                    linked_bytes -= st.st_size

    if not args.dry_run:
        for dirpath, dirnames, filenames in os.walk(html_root, topdown=False):
            if dirpath != html_root and not dirnames and not filenames:
                try:
                    os.rmdir(dirpath)
                except OSError:
                    pass

    label = "would free" if args.dry_run else "freed"
    print("Report tidying: removed %d unreferenced duplicate(s), %s %s; hard-linked %d "
          "referenced duplicate(s), %s a further %s; kept %d file(s) held only by the report."
          % (removed, label, human(removed_bytes), linked, label, human(linked_bytes), unique))
    return 0


if __name__ == "__main__":
    sys.exit(main())
