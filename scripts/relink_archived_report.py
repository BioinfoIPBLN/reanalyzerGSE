#!/usr/bin/env python3
import html
import os
import re
import sys
import urllib.parse

if len(sys.argv) < 4:
    sys.exit("usage: relink_archived_report.py REPORT_DIR ARCHIVE ARCHIVED_PATH...")

report_dir, archive = sys.argv[1], sys.argv[2]
archived = {p[2:] if p.startswith("./") else p for p in sys.argv[3:]}
archived = {p.split("/", 1)[0] for p in archived}
link = re.compile(r'<a href="([^"]+)" download="[^"]*"([^>]*)>(.*?)</a>', re.S)

for fname in sorted(os.listdir(report_dir)):
    if not fname.endswith(".html"):
        continue
    path = os.path.join(report_dir, fname)
    with open(path, encoding="utf-8", errors="surrogateescape") as fh:
        text = fh.read()
    changed = 0

    def relink(m):
        global changed
        href, attrs, label = m.groups()
        target = urllib.parse.unquote(html.unescape(href))
        if target.split("/", 1)[0] not in archived or os.path.exists(os.path.join(report_dir, target)):
            return m.group(0)
        changed += 1
        return (f'<a href="{archive}" download="{archive}"{attrs} title="{html.escape(target)}">'
                f"{label} (in {archive})</a>")

    new = link.sub(relink, text)
    if changed:
        tmp = path + ".relink_tmp"
        with open(tmp, "w", encoding="utf-8", errors="surrogateescape") as fh:
            fh.write(new)
        os.replace(tmp, path)
        print(f"{fname}: {changed} download link(s) now point to {archive}")
