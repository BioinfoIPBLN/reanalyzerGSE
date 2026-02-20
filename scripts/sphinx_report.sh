#!/bin/bash

path=$1
project_name=$2

mkdir -p $path/sphinx_report; cd $path/sphinx_report 


######### Sphinx quickstart
sphinx-quickstart -a $(echo $(whoami)_reanalyzerGSE) -l en -p $project_name -r "" --no-sep &> sphinx.log


######### Modify conf.py
sed -i "s,html_theme =.*,html_theme = 'sphinxdoc',g" conf.py

echo -e "\n
html_extra_path = [
    \"../multiqc_out\",
    \"../miARma_out0/star_results/multibamqc_results\",
    \"../final_results_reanalysis0_$project_name/QC_and_others\",
    \"../final_results_reanalysis0_$project_name/DGE\",
    \"../final_results_reanalysis0_$project_name\"
]
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

        # Resolve the absolute path of the directory
        directory = os.path.abspath(directory)

        # Check if the directory exists
        if not os.path.isdir(directory):
            error_node = nodes.paragraph(text=f\"Directory does not exist: {directory}\")
            return [error_node]

        # Search for files matching the pattern
        all_files = os.listdir(directory)
        matched_files = [file for file in all_files if fnmatch.fnmatch(file, pattern)]

        # If matching files are empty, return a warning node
        if not matched_files:
            warning_node = nodes.paragraph(text=f\"No files match the pattern: {pattern} in {directory}\")
            return [warning_node]

        file_nodes = []
        for file in matched_files:
            # Ignore files containing \"annotation\" or \"Gene_IDs\" or \"fdr\"
            if \"annotation\" in file or \"Gene_IDs\" in file or \"fdr\" in file:
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
            # Read the file as a tab-separated file (assume .txt is tab-delimited)
            data = pd.read_csv(file_path, sep=\"\\\t\", header=None)

            # Ensure numeric conversion for columns 3 and 6
            for col in [2, 5]:
                data[col] = pd.to_numeric(data[col], errors=\"coerce\")

            # Drop rows with invalid numeric values in columns 3 or 6
            data = data.dropna(subset=[2, 5])

            # Filter rows where the 6th column < 0.05
            degs_filtered = data[data[5] < 0.05]

            # Up-regulated DEGs (3rd column > 0)
            degs_up = degs_filtered[degs_filtered[2] > 0]
            num_degs_up = len(degs_up)

            # Down-regulated DEGs (3rd column < 0)
            degs_down = degs_filtered[degs_filtered[2] < 0]
            num_degs_down = len(degs_down)

            # Sort by 3rd column, take top 10 and bottom 10 rows
            sorted_degs = degs_filtered.sort_values(by=2)
            head_tail_degs = pd.concat([sorted_degs.head(10), sorted_degs.tail(10)])

            # Add information to nodes
            caption_node = nodes.paragraph(text=f\"Contents of {file_name}:\")            
            deg_nodes.append(nodes.paragraph(text=f\"Contents of {file_name}:\"))
            deg_nodes.append(nodes.paragraph(text=f\"{num_degs_up} DEGs up\"))
            deg_nodes.append(nodes.paragraph(text=f\"{num_degs_down} DEGs down\"))
            deg_nodes.append(nodes.paragraph(text=f\"Total number of DEGs: {num_degs_up + num_degs_down}\"))
            deg_nodes.append(nodes.paragraph(text=f\"Top 10 DEGs in each sense:\"))

            literal_node = nodes.literal_block()
            literal_node['language'] = \"text\"
            literal_node += nodes.Text(head_tail_degs.to_string(index=False, header=False))
            deg_nodes.append(literal_node)

            # Add a download link (HTML)
            download_html = f'<p><a href=\"{file_name}\" download=\"{file_name}\" class=\"btn btn-primary\">Download {file_name}</a></p>'
            raw_node = nodes.raw(\"\", download_html, format=\"html\")
            deg_nodes.append(raw_node)

        except Exception as e:
            error_node = nodes.paragraph(text=f\"Error processing DEGs in {file_name}: {e}\")
            deg_nodes.append(error_node)

        return deg_nodes

    def process_text(self, file_path, file_name):
        \"\"\"Process text files.\"\"\"
        try:
            literal_node = nodes.literal_block()
            literal_node['language'] = \"text\"
            with open(file_path, \"r\") as f:
                literal_node += nodes.Text(f.read())
            caption_node = nodes.paragraph(text=f\"Contents of {file_name}:\")
            return [caption_node, literal_node]
        except Exception as e:
            return [nodes.paragraph(text=f\"Error processing text file {file_name}: {e}\")]


def setup(app):
    app.add_directive(\"include_matching_files\", IncludeMatchingFiles)
" >> conf.py


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
.. literalinclude:: ../GEO_info/samples_info.txt
.. literalinclude:: ../GEO_info/organism.txt
.. literalinclude:: ../final_results_reanalysis0_$project_name/QC_and_others/reads_numbers.txt
.. literalinclude:: ../miARma_out0/Pre_fastqc_results/list_of_files.txt

.. index:: Samples



Summary of comparisons and covariables
------------------------------------------------------------------------------------
.. include_matching_files:: design_possible_*.txt ../GEO_info/

Covariables or potential batch effect:

.. literalinclude:: ../GEO_info/batch_vector.txt
.. literalinclude:: ../GEO_info/batch_biological_variable.txt
.. literalinclude:: ../GEO_info/covariables.txt
.. index:: Comparisons



Layout and strandedness
------------------------------------------------------------------------------------
.. literalinclude:: ../library_layout_info.txt
.. literalinclude:: ../strand_info.txt
.. index:: Layout



Counts
------------------------------------------------------------------------------------
Please use the following links:

.. raw:: html
   
   <a href=\"RPKM_counts_genes_log2_0.1_categ.txt\" target=\"_blank\">Click to get RPKM counts (log2 + 0.1)</a><br>
   <a href=\"TPM_counts_genes_log2_0.1.txt\" target=\"_blank\">Click to get TPM counts (log2 + 0.1)</a>

If requested, please go to \"$project_name/final_results_reanalysis0_$project_name/violin\" to check out the figures showing the transcriptional profiles of genes of interest. You may also find the tables \"_annotation.txt\" including the gene annotation available

.. index:: Counts



Differentally Expressed Genes
------------------------------------------------------------------------------------
List of comparisons: 
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../final_results_reanalysis0_$project_name/DGE/list_comp.txt

DEGs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include_matching_files:: DGE_analysis_comp*.txt ../final_results_reanalysis0_$project_name/DGE/ degs

.. index:: DEGs



Volcano plots
------------------------------------------------------------------------------------

.. include_matching_files:: Volcano_plot_*.pdf ../final_results_reanalysis0_$project_name/DGE/

.. index:: Volcano



Functional enrichment analyses
------------------------------------------------------------------------------------
If requested, please use the following :download:\`link <../final_results_reanalysis0_$project_name/DGE/funct_enrichment_analyses.tar.gz>\`

.. index:: Funct_enrich



Housekeeping genes
------------------------------------------------------------------------------------
.. include_matching_files:: HK_genes_*_10.txt ../final_results_reanalysis0_$project_name/DGE

.. index:: Housekeeping genes



QC analyses
------------------------------------------------------------------------------------
Please use the following links:

.. raw:: html
   
   <a href=\"multiqc_report.html\" target=\"_blank\">Click to open report by MultiQC</a><br>
   <a href=\"multisampleBamQcReport.html\" target=\"_blank\">Click to open Multi-sample BAM QC by Qualimap</a><br>
   <a href=\"${project_name}_norm_QC.pdf\" target=\"_blank\">Click to open PDF with multiple QC figures</a><br>
   <a href=\"${project_name}_adjusted_QC.pdf\" target=\"_blank\">Click to open PDF with multiple QC figures if batch correction/adjusted counts</a>

.. index:: QC analyses



Genome browser visualization (IGV)
------------------------------------------------------------------------------------
Alignment files (.bam) and coverage tracks (.bw) generated by reanalyzerGSE can be visualized using a genome browser, for example the Integrative Genomics Viewer (IGV). You can use the desktop application or the online version:

.. raw:: html

   <a href=\"https://igv.org/app/\" target=\"_blank\">Open IGV online (igv.org/app)</a><br>
   <a href=\"https://igv.org/doc/desktop/\" target=\"_blank\">Download IGV desktop</a>

A local R Shiny app (``igvShinyApp.R``) has also been placed in your results folder, with the paths to the reference genome, annotation and BigWig coverage tracks pre-filled from the pipeline. Launch it from an R console with:

.. code-block:: r

   shiny::runApp(\"path/to/final_results_*/igvShinyApp.R\")

Requirements: R packages ``shiny``, ``igvShiny``, ``rtracklayer``, ``shinyjs``. Update the ``genome_name`` variable inside the app if the organism name is not recognised by igvShiny (e.g. use ``mm39`` instead of ``Mus_musculus``).

.. index:: IGV



Interactive exploration with exploreDE
------------------------------------------------------------------------------------
If you chose the option \`\`-eDe yes\`\` to create a SummarizedExperiment object compatible with \`exploreDE <https://zenodo.org/records/13927692>\`_, you can use it to interactively explore your differential expression results.

.. raw:: html

   <a href=\"https://shiny-public.fgcz.uzh.ch/app/exploreLocalDE\" target=\"_blank\">Open exploreLocalDE (Shiny app)</a>

The generated .qs2 file can be found in the DGE results folder and can be loaded into exploreDE/exploreLocalDE for interactive visualization of DE results, pathway analyses, and more.

.. index:: exploreDE" > index.rst

######### Build
sphinx-build -M html . . &>> sphinx.log

######### Compress the report:
# tar cf - $path/sphinx_report/ | pigz --best > $path/sphinx_report.tar.gz && rm -rf $path/sphinx_report/

# Make a report in the main folder, correcting the paths:
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
sed -i 's,href="DGE_analysis_comp,href="sphinx_report/html/DGE_analysis_comp,g' $path/final_report.html
sed -i 's,src="_static/,src="sphinx_report/html/_static/,g' $path/final_report.html

echo -e "\n\nfinal_report.html contains the report with the main results generated by Sphinx!\n\n"
