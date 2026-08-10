#!/bin/bash

path=$1
project_name=$2

final_dir_name=$(basename $(find $path -maxdepth 1 -type d \( -name "final_results_reanalysis0_*" -o -name "final_results_reanalysis0" \) | head -1) 2>/dev/null)

if [ -z "$final_dir_name" ]; then
	echo "WARNING: Could not find final_results_reanalysis0 directory in $path. Sphinx report may be incomplete."
	final_dir_name="final_results_reanalysis0"
fi

mkdir -p $path/sphinx_report; cd $path/sphinx_report 


######### Sphinx quickstart
sphinx-quickstart -a $(echo $(whoami)_reanalyzerGSE) -l en -p $project_name -r "" --no-sep &> sphinx.log


######### Modify conf.py
sed -i "s,html_theme =.*,html_theme = 'sphinxdoc',g" conf.py

rnaseqqc_links=""
bamqc_links=""
for aligner_dir in "$path/miARma_out0"/*_results; do
    if [ -d "$aligner_dir/rnaseqqc_results" ]; then
        for qrep in "$aligner_dir/rnaseqqc_results"/*/qualimapReport.html; do
            if [ -f "$qrep" ]; then
                sample_name=$(basename $(dirname "$qrep"))
                display_name=$(echo "$sample_name" | sed -E 's/_(hisat2|STAR)\.bam//')
                rnaseqqc_links="${rnaseqqc_links}   <a href=\"sphinx_report/html/rnaseqqc_results/${sample_name}/qualimapReport.html\" target=\"_blank\">Click to open RNASeq QC for ${display_name}</a><br>
"
            fi
        done
    fi
    if [ -d "$aligner_dir/bamqc_results" ]; then
        for qrep in "$aligner_dir/bamqc_results"/*/qualimapReport.html; do
            if [ -f "$qrep" ]; then
                sample_name=$(basename $(dirname "$qrep"))
                display_name=$(echo "$sample_name" | sed -E 's/_(hisat2|STAR)\.bam//')
                bamqc_links="${bamqc_links}   <a href=\"sphinx_report/html/bamqc_results/${sample_name}/qualimapReport.html\" target=\"_blank\">Click to open BAM QC for ${display_name}</a><br>
"
            fi
        done
    fi
done

echo -e "\n
import glob
import os

html_extra_path = [
    \"../multiqc_out\",
    \"../${final_dir_name}/QC_and_others\",
    \"../${final_dir_name}/QC_and_others/qualimap_rnaseqqc_results\",
    \"../${final_dir_name}/DGE\",
    \"../${final_dir_name}\"
]
html_extra_path.extend(glob.glob(\"../miARma_out0/*_results\"))
html_extra_path.extend(glob.glob(\"../miARma_out0/*_results/multibamqc_results\"))
html_extra_path.extend(glob.glob(\"../miARma_out0/*_results/rnaseqqc_results\"))
html_extra_path.extend(glob.glob(\"../miARma_out0/*_results/bamqc_results\"))
html_extra_path.extend(glob.glob(\"../preliminar_rrna_qc\"))
html_extra_path = [p for p in html_extra_path if os.path.exists(p)]
" >> conf.py

echo -e "\n
from docutils import nodes
from sphinx.util.docutils import SphinxDirective
import os
import fnmatch
import pandas as pd  # Requires pandas for easy CSV processing


class IncludeMatchingFiles(SphinxDirective):
    required_arguments = 1  # The file pattern is required
    optional_arguments = 2  # Directory and optional mode (e.g., \"degs\")
    has_content = False

    def run(self):
        # Handle arguments
        pattern = self.arguments[0]
        directory = self.arguments[1] if len(self.arguments) > 1 else \".\"
        mode = self.arguments[2] if len(self.arguments) > 2 else None

        # Resolve the path relative to the sphinx_report build directory, not the absolute system path
        directory = os.path.normpath(os.path.join(os.getcwd(), directory))

        # Check if the directory exists
        if not os.path.isdir(directory):
            error_node = nodes.paragraph(text=f\"Directory does not exist: {directory}\")
            return [error_node]

        # Search for files matching the pattern
        all_files = os.listdir(directory)
        matched_files = [file for file in all_files if fnmatch.fnmatch(file, pattern)]

        # If matching files are empty, str error without absolute path
        if not matched_files:
            original_dir = self.arguments[1] if len(self.arguments) > 1 else \".\"
            warning_node = nodes.paragraph(text=f\"No files match the pattern: {pattern} in {original_dir}\")
            return [warning_node]

        file_nodes = []
        for file in matched_files:
            # Ignore files containing \"annotation\" or \"Gene_IDs\" or \"fdr\" or \"merged\"
            if \"annotation\" in file or \"Gene_IDs\" in file or \"fdr\" in file or \"merged\" in file:
                continue

            file_path = os.path.join(directory, file)
            if pattern.endswith(\".pdf\"):
                file_nodes.extend(self.process_pdf(file_path, file))
            elif mode == \"degs\":
                file_nodes.extend(self.process_degs(file_path, file))
            else:
                file_nodes.extend(self.process_text(file_path, file))

        return file_nodes

    def process_pdf(self, file_path, file_name):
        \"\"\"Provide a clickable link to open/download the PDF file.\"\"\"
        try:
            caption_node = nodes.paragraph(text=f\"Found {file_name}:\")
            link_html = f\"<a href=\\\"{file_name}\\\" target=\\\"_blank\\\">Open PDF</a>\"
            raw_node = nodes.raw(\"\", link_html, format=\"html\")
            return [caption_node, raw_node]
        except Exception as e:
            return [nodes.paragraph(text=f\"Error processing PDF file {file_name}: {e}\")]

    def process_degs(self, file_path, file_name):
        \"\"\"Process the file in DEGs mode.\"\"\"
        deg_nodes = []

        try:
            # Read the file as a tab-separated file and use the first row as header
            data = pd.read_csv(file_path, sep=\"\\\t\", header=0)

            # Ensure numeric conversion for columns 3 and 6 (index 2 and 5)
            # Log2FoldChange is usually column 3, padj or FDR is usually column 6
            data.iloc[:, 2] = pd.to_numeric(data.iloc[:, 2], errors=\"coerce\")
            data.iloc[:, 5] = pd.to_numeric(data.iloc[:, 5], errors=\"coerce\")

            # Drop rows with invalid numeric values in columns 3 or 6
            data = data.dropna(subset=[data.columns[2], data.columns[5]])

            # Filter rows where the 6th column < 0.05
            degs_filtered = data[data.iloc[:, 5] < 0.05]

            # Up-regulated DEGs (3rd column > 0)
            degs_up = degs_filtered[degs_filtered.iloc[:, 2] > 0]
            num_degs_up = len(degs_up)

            # Down-regulated DEGs (3rd column < 0)
            degs_down = degs_filtered[degs_filtered.iloc[:, 2] < 0]
            num_degs_down = len(degs_down)

            # Top 10 in EACH sense, taken from the separated up/down sets so an
            # empty sense stays empty. (Previously head(10)+tail(10) of the
            # combined set overlapped when there were < 20 DEGs, duplicating the
            # same rows for both senses, e.g. 2 up / 0 down shown as 2+2.)
            top_down = degs_down.sort_values(by=data.columns[2]).head(10)
            top_up   = degs_up.sort_values(by=data.columns[2]).tail(10)
            head_tail_degs = pd.concat([top_down, top_up])

            # Add information to nodes
            caption_node = nodes.paragraph()
            strong_node = nodes.strong(text=f\"Contents of {file_name}:\")
            caption_node += strong_node
            deg_nodes.append(caption_node)

            deg_nodes.append(nodes.paragraph(text=f\"{num_degs_up} DEGs up\"))
            deg_nodes.append(nodes.paragraph(text=f\"{num_degs_down} DEGs down\"))
            deg_nodes.append(nodes.paragraph(text=f\"Total number of DEGs: {num_degs_up + num_degs_down}\"))
            deg_nodes.append(nodes.paragraph(text=f\"Top 10 DEGs in each sense:\"))

            literal_node = nodes.literal_block()
            literal_node['language'] = \"text\"
            literal_node += nodes.Text(head_tail_degs.to_string(index=False, header=True))
            deg_nodes.append(literal_node)

            # Add a download link (HTML)
            xlsx_name = os.path.splitext(file_name)[0] + \".xlsx\"
            xlsx_path = os.path.join(os.path.dirname(file_path), xlsx_name)
            xlsx_link = f' (<a href=\"{xlsx_name}\" download=\"{xlsx_name}\">xlsx version</a>)' if os.path.exists(xlsx_path) else \"\"
            download_html = f'<p><a href=\"{file_name}\" download=\"{file_name}\" class=\"btn btn-primary\">Download {file_name}</a>{xlsx_link}</p>'
            raw_node = nodes.raw(\"\", download_html, format=\"html\")
            deg_nodes.append(raw_node)

            # --- Extract top 10 from merged annotated file if exists ---
            base_name = os.path.splitext(file_name)[0]
            # Assumes format: DGE_analysis_compX_merged_RPKM.txt or similar
            merged_pattern = os.path.join(os.path.dirname(file_path), f\"{base_name}_merged_*.txt\")
            merged_files = glob.glob(merged_pattern)
            
            if merged_files:
                # Prefer RPKM for the table preview, if multiple, else just take the first
                preview_file_path = next((f for f in merged_files if \"RPKM\" in f), merged_files[0])
                try:
                    merged_data = pd.read_csv(preview_file_path, sep=\"\\\t\", header=0)
                    
                    # Usually logFC is col 2, FDR is col 5, but let's dynamically find them if possible
                    logfc_col_idx = 2
                    padj_col_idx = 5
                    for i, col in enumerate(merged_data.columns):
                        if \"logFC\" in col or \"log2FoldChange\" in col:
                            logfc_col_idx = i
                        elif \"FDR\" in col or \"padj\" in col or \"P.Value\" in col or \"pvalue\" in col:
                            padj_col_idx = i

                    merged_data.iloc[:, logfc_col_idx] = pd.to_numeric(merged_data.iloc[:, logfc_col_idx], errors=\"coerce\")
                    merged_data.iloc[:, padj_col_idx] = pd.to_numeric(merged_data.iloc[:, padj_col_idx], errors=\"coerce\")
                    merged_data = merged_data.dropna(subset=[merged_data.columns[logfc_col_idx], merged_data.columns[padj_col_idx]])
                    
                    merged_degs_filtered = merged_data[merged_data.iloc[:, padj_col_idx] < 0.05]
                    merged_lfc = merged_data.columns[logfc_col_idx]
                    # Top 10 per sense from the separated up/down sets (see note above).
                    merged_down = merged_degs_filtered[merged_degs_filtered.iloc[:, logfc_col_idx] < 0].sort_values(by=merged_lfc).head(10)
                    merged_up   = merged_degs_filtered[merged_degs_filtered.iloc[:, logfc_col_idx] > 0].sort_values(by=merged_lfc).tail(10)
                    merged_head_tail = pd.concat([merged_down, merged_up])

                    # Determine which expression type is being previewed
                    preview_basename = os.path.basename(preview_file_path)
                    if \"RPKM\" in preview_basename:
                        expr_type = \"RPKM\"
                    elif \"TPM\" in preview_basename:
                        expr_type = \"TPM\"
                    elif \"CPM\" in preview_basename:
                        expr_type = \"CPM\"
                    else:
                        expr_type = \"expression\"
                    deg_nodes.append(nodes.paragraph(text=f\"Top 10 DEGs in each sense (with annotation, {expr_type} expression values):\"))
                    merged_literal_node = nodes.literal_block()
                    merged_literal_node['language'] = \"text\"
                    merged_literal_node += nodes.Text(merged_head_tail.to_string(index=False, header=True))
                    deg_nodes.append(merged_literal_node)

                    # Add download link for all merged files
                    for m_file_path in sorted(merged_files):
                        m_file_name = os.path.basename(m_file_path)
                        m_xlsx_name = os.path.splitext(m_file_name)[0] + \".xlsx\"
                        m_xlsx_path = os.path.join(os.path.dirname(m_file_path), m_xlsx_name)
                        m_xlsx_link = f' (<a href=\"{m_xlsx_name}\" download=\"{m_xlsx_name}\">xlsx version</a>)' if os.path.exists(m_xlsx_path) else \"\"
                        m_download_html = f'<p><a href=\"{m_file_name}\" download=\"{m_file_name}\" class=\"btn btn-primary\">Download {m_file_name}</a>{m_xlsx_link}</p>'
                        deg_nodes.append(nodes.raw(\"\", m_download_html, format=\"html\"))

                except Exception as e:
                    deg_nodes.append(nodes.paragraph(text=f\"Error processing merged DEGs from {preview_file_path}: {e}\"))

        except Exception as e:
            error_node = nodes.paragraph(text=f\"Error processing DEGs in {file_name}: {e}\")
            deg_nodes.append(error_node)

        # Per-comparison AI insight box, shown at the END of THIS comparison's
        # subsection so each 'Contents of ...' block gets its own box (rather than
        # all boxes collected together at the end of the DEGs section). The box was
        # pre-rendered to <comp>.ai_insight.html by sphinx_report.sh (reusing the
        # same ai_box() styling as the design/counts boxes).
        try:
            ai_html = os.path.join(os.path.dirname(file_path),
                                   os.path.splitext(file_name)[0] + \".ai_insight.html\")
            if os.path.isfile(ai_html):
                with open(ai_html, encoding=\"utf-8\", errors=\"replace\") as _fh:
                    _box = _fh.read().strip()
                if _box:
                    deg_nodes.append(nodes.raw(\"\", _box, format=\"html\"))
        except Exception:
            pass

        return deg_nodes

    def process_text(self, file_path, file_name):
        \"\"\"Process text files.\"\"\"
        try:
            literal_node = nodes.literal_block()
            literal_node['language'] = \"text\"
            with open(file_path, \"r\") as f:
                literal_node += nodes.Text(f.read())

            caption_node = nodes.paragraph()
            strong_node = nodes.strong(text=f\"Contents of {file_name}:\")
            caption_node += strong_node

            return [caption_node, literal_node]
        except Exception as e:
            return [nodes.paragraph(text=f\"Error processing text file {file_name}: {e}\")]


def setup(app):
    app.add_directive(\"include_matching_files\", IncludeMatchingFiles)
" >> conf.py


######### Optional AI insight boxes (scripts/llm_insight.py wrote *.ai_insight.md
######### into the DGE dir when -llm_endpoint + ai_insights were enabled). Convert
######### each to a single-line HTML <div> and wrap it in a `.. raw:: html` block;
######### the variables stay empty (nothing shown) when no insight files exist.
ai_box() {   # $1 = *.ai_insight.md -> single-line HTML div on stdout (nothing if absent/empty)
	[ -f "$1" ] || return 0
	python3 - "$1" <<'PYEOF'
import sys, html, re
try:
    lines=[l.rstrip("\n") for l in open(sys.argv[1],encoding="utf-8",errors="replace")]
except Exception:
    sys.exit(0)
lines=[l for l in lines if l.strip()]
if not lines: sys.exit(0)
def conv(s):
    s=html.escape(s)                       # escapes & < > " ' -> no quote survives into the RST
    s=re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
    s=re.sub(r"`(.+?)`", r"<code>\1</code>", s)
    return s
stamp=conv(lines[0]); body="".join("<p style='margin:4px 0;'>%s</p>"%conv(l) for l in lines[1:])
print("<div style='background:#eef6ff;border-left:5px solid #4a90d9;padding:10px 14px;margin:12px 0;"
      "border-radius:4px;font-size:13px;'><div style='color:#5a6b7b;font-size:11px;margin-bottom:6px;'>"
      "\U0001f916 %s</div>%s</div>" % (stamp, body))
PYEOF
}
ai_rst() { local box; box=$(ai_box "$1"); [ -n "$box" ] || return 0; printf '\n.. raw:: html\n\n   %s\n' "$box"; }
ai_design_rst=$(ai_rst "$path/$final_dir_name/DGE/study_design.ai_insight.md")
ai_counts_rst=$(ai_rst "$path/$final_dir_name/DGE/counts.ai_insight.md")
# Per-comparison DEG insight boxes are injected per-subsection inside process_degs()
# (conf.py, below), so each appears at the end of its own "Contents of ..." block
# rather than all together at the end of the DEGs section. Pre-render each one to a
# ready HTML file that conf.py embeds verbatim (same styling as the boxes above).
for aif in "$path/$final_dir_name/DGE/"DGE_analysis_comp*.ai_insight.md; do
	[ -e "$aif" ] || continue
	case "$aif" in *".enrichment_ai_insight.md") continue ;; esac   # those go in the enrichment report
	aibox=$(ai_box "$aif")
	[ -n "$aibox" ] && printf '%s\n' "$aibox" > "${aif%.md}.html"
done


######### Pipeline timing: turn the raw start/end epoch rows in step_times.tsv into a
######### per-step DURATION table (hh:mm:ss) + a total, indented 3 spaces for an RST
######### code-block. Empty if the file is absent; a step with a start but no end yet
######### (e.g. the report step itself) shows "(in progress)".
timing_block=""
if [ -f "$path/step_times.tsv" ]; then
	timing_block=$(awk -F'\t' -v now="$(date +%s)" '
		$3=="start"{ s[$1]=$2; if(!($1 in seen)){ order[++n]=$1; seen[$1]=1 } }
		$3=="end"  { e[$1]=$2 }
		function hms(d,  h,m,sec){ if(d<0) d=0; h=int(d/3600); m=int((d%3600)/60); sec=d%60; return sprintf("%02d:%02d:%02d",h,m,sec) }
		END{
			w=8; for(i=1;i<=n;i++) if(length(order[i])>w) w=length(order[i])
			printf "   %-*s  %s\n", w, "Step", "Duration (hh:mm:ss)"
			printf "   %-*s  %s\n", w, "----", "-------------------"
			first=""; last=0
			for(i=1;i<=n;i++){
				k=order[i]
				if(first=="" && (k in s)) first=s[k]
				if((k in e) && (k in s)){ d=e[k]-s[k]; if(e[k]>last) last=e[k]; dur=hms(d) }
				else if(k in s){ d=now-s[k]; if(now>last) last=now; dur=hms(d) }
				else dur="N/A"
				printf "   %-*s  %s\n", w, k, dur
			}
			if(first!="" && last>0) printf "   %-*s  %s\n", w, "TOTAL", hms(last-first)
		}' "$path/step_times.tsv")
fi


# Ensure layout and strand info exist at root level
if [ ! -f "$path/library_layout_info.txt" ] && [ -f "$path/reads_study_info/library_layout_info.txt" ]; then
    cp "$path/reads_study_info/library_layout_info.txt" "$path/library_layout_info.txt" 2>/dev/null || true
fi
if [ ! -f "$path/strand_info.txt" ] && [ -f "$path/reads_study_info/strand_info.txt" ]; then
    cp "$path/reads_study_info/strand_info.txt" "$path/strand_info.txt" 2>/dev/null || true
fi

layout_inc=""
if [ -f "$path/library_layout_info.txt" ]; then
    layout_inc=".. literalinclude:: ../library_layout_info.txt"
elif [ -f "$path/reads_study_info/library_layout_info.txt" ]; then
    layout_inc=".. literalinclude:: ../reads_study_info/library_layout_info.txt"
fi

strand_inc=""
if [ -f "$path/strand_info.txt" ]; then
    strand_inc=".. literalinclude:: ../strand_info.txt"
elif [ -f "$path/reads_study_info/strand_info.txt" ]; then
    strand_inc=".. literalinclude:: ../reads_study_info/strand_info.txt"
fi

batch_vec_inc=""
[ -f "$path/reads_study_info/batch_vector.txt" ] && batch_vec_inc=".. literalinclude:: ../reads_study_info/batch_vector.txt"

batch_bio_inc=""
if [ -f "$path/reads_study_info/batch_biological_variables.txt" ]; then
    batch_bio_inc=".. literalinclude:: ../reads_study_info/batch_biological_variables.txt"
elif [ -f "$path/reads_study_info/batch_biological_variable.txt" ]; then
    batch_bio_inc=".. literalinclude:: ../reads_study_info/batch_biological_variable.txt"
fi

covar_inc=""
[ -f "$path/reads_study_info/covariables.txt" ] && covar_inc=".. literalinclude:: ../reads_study_info/covariables.txt"

######### Modify index.rst
echo "
Welcome to $project_name report!
####################################################################################

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Sections
==================



Summary of samples and experimental conditions
------------------------------------------------------------------------------------
.. literalinclude:: ../reads_study_info/samples_info.txt
.. literalinclude:: ../reads_study_info/organism.txt
.. literalinclude:: ../$final_dir_name/QC_and_others/reads_numbers.txt
.. literalinclude:: ../miARma_out0/Pre_fastqc_results/list_of_files.txt
.. index:: Samples



Summary of comparisons and covariables
------------------------------------------------------------------------------------
.. include_matching_files:: design_possible_*.txt ../reads_study_info/

Covariables or potential batch effect:

$batch_vec_inc
$batch_bio_inc
$covar_inc
.. index:: Comparisons



Layout and strandedness
------------------------------------------------------------------------------------
$layout_inc
$strand_inc
$ai_design_rst
.. index:: Layout



Counts
------------------------------------------------------------------------------------
Please use the following links:

.. raw:: html
   
   <a href=\"sphinx_report/html/TPM_counts_genes_log2_0.1_categ.txt\" target=\"_blank\">Click to get TPM counts (log2 + 0.1)</a>$(if [ -f "$path/$final_dir_name/TPM_counts_genes_log2_0.1_categ.xlsx" ]; then echo ' (<a href="sphinx_report/html/TPM_counts_genes_log2_0.1_categ.xlsx" target="_blank">xlsx version</a>)'; fi)<br>
   <a href=\"sphinx_report/html/RPKM_counts_genes_log2_0.1_categ.txt\" target=\"_blank\">Click to get RPKM counts (log2 + 0.1)</a>$(if [ -f "$path/$final_dir_name/RPKM_counts_genes_log2_0.1_categ.xlsx" ]; then echo ' (<a href="sphinx_report/html/RPKM_counts_genes_log2_0.1_categ.xlsx" target="_blank">xlsx version</a>)'; fi)<br>
   <a href=\"sphinx_report/html/CPM_counts_genes_log2_0.1_categ.txt\" target=\"_blank\">Click to get CPM counts (log2 + 0.1)</a>$(if [ -f "$path/$final_dir_name/CPM_counts_genes_log2_0.1_categ.xlsx" ]; then echo ' (<a href="sphinx_report/html/CPM_counts_genes_log2_0.1_categ.xlsx" target="_blank">xlsx version</a>)'; fi)

If requested, please go to \"$project_name/$final_dir_name/violin\" to check out the figures showing the transcriptional profiles of genes of interest. You may also find the tables \"_annotation.txt\" including the gene annotation available.$(if [ -f "$path/$final_dir_name/DGE/deResults.qs2" ]; then echo " The exploreLocalDE app may be also used (see below)."; fi)
$ai_counts_rst
.. index:: Counts



Differentally Expressed Genes
------------------------------------------------------------------------------------
List of comparisons: 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include_matching_files:: list_comp.txt ../$final_dir_name/DGE/

DEGs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include_matching_files:: DGE_analysis_comp*.txt ../$final_dir_name/DGE/ degs
.. index:: DEGs



Volcano plots
------------------------------------------------------------------------------------

.. include_matching_files:: Volcano_plot_*.pdf ../$final_dir_name/DGE/

.. index:: Volcano


Functional enrichment analyses
------------------------------------------------------------------------------------
Please use the following links:

$(if [ -f "$path/$final_dir_name/DGE/functional_enrichment_report.html" ]; then echo "
.. raw:: html

   <a href=\"functional_enrichment_report.html\" target=\"_blank\">Click to open Functional Enrichment HTML Report</a><br>
"; fi)
$(if [ -f "$path/$final_dir_name/DGE/funct_enrichment_analyses.tar.gz" ]; then echo "You can also download all raw enrichment result files: :download:\`tar.gz archive <../$final_dir_name/DGE/funct_enrichment_analyses.tar.gz>\`"; fi)
$(if [ ! -f "$path/$final_dir_name/DGE/functional_enrichment_report.html" ] && [ ! -f "$path/$final_dir_name/DGE/funct_enrichment_analyses.tar.gz" ]; then echo "Functional enrichment not requested or not available"; fi)

.. index:: Funct_enrich




Housekeeping genes
------------------------------------------------------------------------------------
.. include_matching_files:: HK_genes_*_10.txt ../$final_dir_name/DGE

.. index:: Housekeeping genes



QC analyses
------------------------------------------------------------------------------------
Please use the following links:

.. raw:: html
   
   <a href=\"sphinx_report/html/multiqc_report.html\" target=\"_blank\">Click to open report by MultiQC</a><br>
${rnaseqqc_links}${bamqc_links}$(if [ -f "$path/miARma_out0/star_results/multibamqc_results/multisampleBamQcReport.html" ] || [ -f "$path/$final_dir_name/multisampleBamQcReport.html" ] || [ $(find "$path/miARma_out0" -name "multisampleBamQcReport.html" 2>/dev/null | wc -l) -gt 0 ]; then echo "   <a href=\"sphinx_report/html/multisampleBamQcReport.html\" target=\"_blank\">Click to open Multi-sample BAM QC by Qualimap</a><br>"; fi)$(if [ -f "$path/$final_dir_name/QC_and_others/${project_name}_norm_QC_commented.pdf" ]; then echo "   <a href=\"sphinx_report/html/${project_name}_norm_QC_commented.pdf\" target=\"_blank\">Click to open PDF with multiple QC figures (with AI commentary)</a><br>"; else echo "   <a href=\"sphinx_report/html/${project_name}_norm_QC.pdf\" target=\"_blank\">Click to open PDF with multiple QC figures</a>$(if [ -f "$path/$final_dir_name/QC_and_others/QC_AI_commentary_slides.pdf" ]; then echo " - <a href=\"sphinx_report/html/QC_AI_commentary_slides.pdf\" target=\"_blank\">AI commentary slides</a>"; fi)<br>"; fi)
$(if [ -f "$path/$final_dir_name/QC_and_others/${project_name}_adjusted_QC_commented.pdf" ]; then echo "   <a href=\"sphinx_report/html/${project_name}_adjusted_QC_commented.pdf\" target=\"_blank\">Click to open PDF with multiple QC figures if batch correction/adjusted counts (with AI commentary)</a><br>"; elif [ -f "$path/$final_dir_name/QC_and_others/${project_name}_adjusted_QC.pdf" ]; then echo "   <a href=\"sphinx_report/html/${project_name}_adjusted_QC.pdf\" target=\"_blank\">Click to open PDF with multiple QC figures if batch correction/adjusted counts</a><br>"; fi)

.. index:: QC analyses



$(if [ -d "$path/preliminar_rrna_qc" ] && [ -f "$path/preliminar_rrna_qc/rRNA_mapping_summary_R1.tsv" ]; then echo "
Preliminary rRNA QC
------------------------------------------------------------------------------------
rRNA contamination was assessed by mapping R1 reads against rRNA reference databases using Bowtie2.

.. raw:: html

   <a href=\"sphinx_report/html/rRNA_mapping_barplot.html\" target=\"_blank\">Click to open interactive rRNA mapping barplot</a><br>

.. literalinclude:: ../preliminar_rrna_qc/rRNA_mapping_summary_R1.tsv

.. index:: rRNA QC
"; fi)

$(if [ -f "$path/step_times.tsv" ]; then echo "
Pipeline timing
------------------------------------------------------------------------------------
Wall-clock time per pipeline step:

.. code-block:: text

$timing_block
$(if [ -f "$path/$final_dir_name/QC_and_others/pipeline_gantt_chart.pdf" ]; then echo "
A Gantt chart of the same timeline is also available:

.. raw:: html

   <a href=\"sphinx_report/html/pipeline_gantt_chart.pdf\" target=\"_blank\">Click to open Pipeline Gantt Chart (PDF)</a><br>
"; fi)

.. index:: Pipeline timing
"; fi)


Genome browser visualization (IGV)
------------------------------------------------------------------------------------
Alignment files (.bam) and coverage tracks (.bw) generated by reanalyzerGSE can be visualized using a genome browser, for example the Integrative Genomics Viewer (IGV). You can use the desktop application or the online version:

.. raw:: html

   <a href=\"https://igv.org/app/\" target=\"_blank\">Open IGV online (igv.org/app)</a><br>
   <a href=\"https://igv.org/doc/desktop/\" target=\"_blank\">Download IGV desktop</a>

A local R Shiny app (``igvShinyApp.R``) has also been placed in your results folder, with the paths to the reference genome, annotation and BigWig coverage tracks pre-filled from the pipeline. Launch it from an R console with:

.. code-block:: r

   shiny::runApp(\"$final_dir_name/igvShinyApp.R\")

.. index:: IGV



$(if [ -f "$path/$final_dir_name/DGE/deResults.qs2" ]; then echo "
Interactive exploration
------------------------------------------------------------------------------------
A SummarizedExperiment object has been generated for this analysis, so you can interactively explore your differential expression results, pathway analyses and more in the exploreLocalDE Shiny app.

.. raw:: html

   <a href=\"https://shiny-public.fgcz.uzh.ch/app/exploreLocalDE\" target=\"_blank\">Open exploreLocalDE Shiny app</a>

Load the following object into the app: :download:\`deResults.qs2 <../$final_dir_name/DGE/deResults.qs2>\` (also available in the :file:\`${final_dir_name}/DGE/\` results folder).

.. index:: exploreLocalDE
"; fi)
" > index.rst

######### Build
sphinx-build -M html . . &>> sphinx.log

######### Compress the report:
# tar cf - $path/sphinx_report/ | pigz --best > $path/sphinx_report.tar.gz && rm -rf $path/sphinx_report/

# Make a report in the main folder, correcting the paths:
if [ -f "$path/sphinx_report/html/index.html" ]; then
	cp $path/sphinx_report/html/index.html $path/final_report.html
	sed -i 's,href="_static/,href="sphinx_report/html/_static/,g' $path/final_report.html
	sed -i 's,href="_source,href="sphinx_report/html/_source,g' $path/final_report.html
	sed -i 's,href="_downloads,href="sphinx_report/html/_downloads,g' $path/final_report.html
	sed -i 's,href="genindex,href="sphinx_report/html/genindex,g' $path/final_report.html
	sed -i 's,href="search,href="sphinx_report/html/search,g' $path/final_report.html
	sed -i "s,href=\"$project_name,href=\"sphinx_report/html/$project_name,g" $path/final_report.html
	sed -i 's,href="multi,href="sphinx_report/html/multi,g' $path/final_report.html
	sed -i 's,href="Volcano,href="sphinx_report/html/Volcano,g' $path/final_report.html
	sed -i 's,href="TPM,href="sphinx_report/html/TPM,g' $path/final_report.html
	sed -i 's,href="RPKM,href="sphinx_report/html/RPKM,g' $path/final_report.html
	sed -i 's,href="CPM,href="sphinx_report/html/CPM,g' $path/final_report.html
	sed -i 's,href="DGE_analysis_comp,href="sphinx_report/html/DGE_analysis_comp,g' $path/final_report.html
	sed -i 's,href="functional_enrichment_report,href="sphinx_report/html/functional_enrichment_report,g' $path/final_report.html
	sed -i 's,src="_static/,src="sphinx_report/html/_static/,g' $path/final_report.html
else
	echo "WARNING: Sphinx build did not produce index.html. Check sphinx_report/sphinx.log for details."
fi

echo -e "\n\nfinal_report.html contains the report with the main results generated by Sphinx!\n\n"
