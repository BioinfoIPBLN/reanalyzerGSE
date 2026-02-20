library(shiny)
library(igvShiny)
library(rtracklayer)
library(shinyjs)

# igvShinyApp.R — Interactive IGV genome browser for reanalyzerGSE output
#
# Default paths below are set automatically by the reanalyzerGSE pipeline
# when this script is copied to the final_results folder. You can also
# change them manually in the UI text inputs after launching the app.
#
# NOTE: 'genome_name' should match an igvShiny-recognised short name
# (e.g. "hg38", "mm39", "mm10", "hg19", "rn7"). The pipeline sets a
# placeholder; update it if the pre-filled value is not recognised by igvShiny.

# ---- Paths pre-filled by reanalyzerGSE (replaced at deploy time) ----
default_fasta   <- "/path/to/reference_genome.fa"
default_gtf     <- "/path/to/annotation.gtf"
default_bw_dir  <- "/path/to/bigwig_folder/"
genome_name     <- "GENOME_NAME_PLACEHOLDER"   # e.g. "mm39", "hg38"
# ---------------------------------------------------------------------

# Derive a sensible initial locus from a random contig in the GTF (first 1000 bp)
get_initial_locus <- function(gtf_path) {
  tryCatch({
    if (!file.exists(gtf_path)) return("all")
    # Read non-comment lines and collect unique contig names
    con <- file(gtf_path, "r")
    on.exit(close(con))
    contigs <- character(0)
    while (length(line <- readLines(con, n = 1)) > 0) {
      if (!startsWith(line, "#")) {
        contig <- strsplit(line, "\t")[[1]][1]
        if (!is.na(contig) && nchar(contig) > 0) contigs <- c(contigs, contig)
      }
      # Stop early once we have enough unique contigs (no need to read entire GTF)
      if (length(unique(contigs)) >= 50) break
    }
    contigs <- unique(contigs)
    if (length(contigs) == 0) return("all")
    chosen <- sample(contigs, 1)
    paste0(chosen, ":1-1000")
  }, error = function(e) "all")
}

default_locus <- get_initial_locus(default_gtf)

ui <- fluidPage(
  useShinyjs(),
  tags$style(HTML("
    .full-width { width: 100% !important; }
    .partial-width { width: 75% !important; }
    #sidebar { transition: all 0.3s; }
  ")),
  titlePanel("IGV Genome Browser"),
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      textInput("fasta_path", "FASTA file path:",
                value = default_fasta),
      textInput("gtf_path", "GTF file path:",
                value = default_gtf),
      textInput("bw_folder", "BigWig folder:",
                value = default_bw_dir),
      textInput("locus", "Locus:", value = default_locus),
      actionButton("load", "Load Files", class = "btn-primary"),
      hr(),
      uiOutput("bw_files_found")
    ),
    mainPanel(
      id = "main",
      actionButton("toggle_sidebar", "☰ Settings",
                   class = "btn-sm btn-default",
                   style = "margin-bottom: 10px;"),
      igvShinyOutput("igvShiny_0", height = "700px"),
      width = 12
    )
  )
)

server <- function(input, output, session) {
  
  loaded <- reactiveValues(
    fasta = NULL,
    fai   = NULL,
    gtf   = NULL,
    bws   = NULL,
    locus = NULL
  )
  
  # Toggle sidebar
  observeEvent(input$toggle_sidebar, {
    toggle("sidebar")
    runjs("
      if ($('#sidebar').is(':visible')) {
        $('#main').css('width', '75%');
      } else {
        $('#main').css('width', '100%');
      }
    ")
  })
  
  bw_files <- reactive({
    req(input$bw_folder)
    folder <- input$bw_folder
    if (!dir.exists(folder)) return(NULL)
    list.files(folder, pattern = "\\.bw$", full.names = TRUE)
  })
  
  output$bw_files_found <- renderUI({
    files <- bw_files()
    if (is.null(files) || length(files) == 0) {
      helpText("No .bw files found in folder.")
    } else {
      helpText(paste0(length(files), " BigWig file(s) found: ",
                      paste(basename(files), collapse = ", ")))
    }
  })
  
  index_fasta <- function(fasta_path) {
    fai_path <- paste0(fasta_path, ".fai")
    if (!file.exists(fai_path)) {
      message("Indexing FASTA...")
      system2("samtools", args = c("faidx", fasta_path))
    }
    if (!file.exists(fai_path)) stop("Failed to create FASTA index")
    return(fai_path)
  }
  
  read_bw_as_bedgraph <- function(bw_path) {
    bw <- as.data.frame(rtracklayer::import(bw_path))
    data.frame(
      chrom = as.character(bw$seqnames),
      start = bw$start - 1,
      end   = bw$end,
      score = bw$score
    )
  }
  
  read_gtf_as_gff3 <- function(gtf_path, locus = NULL) {
    message("Reading GTF for locus: ", locus)
    
    if (!is.null(locus)) {
      locus_clean <- gsub(",", "", locus)
      chrom  <- sub(":.*", "", locus_clean)
      coords <- sub(".*:", "", locus_clean)
      start  <- as.integer(sub("-.*", "", coords))
      end    <- as.integer(sub(".*-", "", coords))
      
      message("  Parsed - chrom: ", chrom, " start: ", start, " end: ", end)
      
      if (is.na(start) || is.na(end)) {
        message("  Locus parse failed, loading without region filter")
        gtf <- as.data.frame(rtracklayer::import(gtf_path))
      } else {
        which <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))
        gtf   <- as.data.frame(rtracklayer::import(gtf_path, which = which))
      }
    } else {
      gtf <- as.data.frame(rtracklayer::import(gtf_path))
    }
    
    message("  GTF rows loaded: ", nrow(gtf))
    if (nrow(gtf) == 0) return(NULL)
    
    data.frame(
      seqname   = as.character(gtf$seqnames),
      source    = ".",
      feature   = as.character(gtf$type),
      start     = gtf$start,
      end       = gtf$end,
      score     = ".",
      strand    = as.character(gtf$strand),
      frame     = ".",
      attribute = paste0(
        "gene_id=",    ifelse(is.na(gtf$gene_id),   ".", gtf$gene_id),
        ";gene_name=", ifelse(is.na(gtf$gene_name), ".", gtf$gene_name)
      ),
      stringsAsFactors = FALSE
    )
  }
  
  load_annotation <- function(locus) {
    tryCatch({
      showNotification("Loading gene annotation...", type = "message", duration = 3)
      tbl.gff3 <- read_gtf_as_gff3(loaded$gtf, locus = locus)
      
      if (!is.null(tbl.gff3) && nrow(tbl.gff3) > 0) {
        loadGFF3TrackFromLocalData(
          session,
          id                     = "igvShiny_0",
          trackName              = "Genes",
          tbl.gff3               = tbl.gff3,
          color                  = "blue",
          colorTable             = list(),
          colorByAttribute       = NA_character_,
          displayMode            = "EXPANDED",
          trackHeight            = 150,
          visibilityWindow       = 500000,
          deleteTracksOfSameName = TRUE
        )
        message("Annotation loaded: ", nrow(tbl.gff3), " features")
      } else {
        message("No annotation features in this region")
        showNotification("No genes found in this region", type = "warning", duration = 3)
      }
    }, error = function(e) {
      message("GTF load ERROR: ", e$message)
      showNotification(paste("GTF load failed:", e$message), type = "error")
    })
  }
  
  observeEvent(input$load, {
    req(input$fasta_path, input$gtf_path, input$bw_folder)
    
    fasta <- input$fasta_path
    gtf   <- input$gtf_path
    bws   <- bw_files()
    locus <- input$locus
    
    for (f in c(fasta, gtf)) {
      if (!file.exists(f)) {
        showNotification(paste("File not found:", f), type = "error")
        return()
      }
    }
    if (is.null(bws) || length(bws) == 0) {
      showNotification("No BigWig files found in folder.", type = "error")
      return()
    }
    
    fai <- tryCatch(
      index_fasta(fasta),
      error = function(e) {
        showNotification(paste("FASTA indexing failed:", e$message), type = "error")
        return(NULL)
      }
    )
    if (is.null(fai)) return()
    
    loaded$fasta <- fasta
    loaded$fai   <- fai
    loaded$gtf   <- gtf
    loaded$bws   <- bws
    loaded$locus <- locus
    
    output$igvShiny_0 <- renderIgvShiny({
      igvShiny(
        genomeOptions = parseAndValidateGenomeSpec(
          genomeName   = genome_name,
          initialLocus = locus,
          stockGenome  = FALSE,
          dataMode     = "localFiles",
          fasta        = fasta,
          fastaIndex   = fai
        ),
        height = "700px"
      )
    })
    
    # Hide sidebar and expand main to full width
    hide("sidebar")
    runjs("$('#main').css('width', '100%');")
  })
  
  observeEvent(input$igvReady, {
    req(loaded$bws)
    message("IGV ready - loading tracks")
    
    load_annotation(loaded$locus)
    
    bws    <- loaded$bws
    colors <- colorRampPalette(c("#2166AC", "#F46D43", "#4DAC26",
                                 "#8B008B", "#FF8C00"))(length(bws))
    
    for (i in seq_along(bws)) {
      tryCatch({
        message("Loading BigWig: ", basename(bws[i]))
        tbl.bw <- read_bw_as_bedgraph(bws[i])
        loadBedGraphTrack(
          session,
          id        = "igvShiny_0",
          trackName = sub("\\.bam\\.bw$|\\.bw$", "", basename(bws[i])),
          tbl       = tbl.bw,
          color     = colors[i],
          autoscale = TRUE
        )
        message("  OK: ", nrow(tbl.bw), " rows")
      }, error = function(e) {
        message("  ERROR: ", e$message)
        showNotification(paste("BigWig load failed:", basename(bws[i]), e$message),
                         type = "error")
      })
    }
    
    showNotification(
      paste0("Loaded annotation + ", length(bws), " BigWig file(s)!"),
      type = "message", duration = 5
    )
  })
  
  observeEvent(input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]], {
    req(loaded$gtf)
    newLoc <- input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]]
    loaded$locus <- newLoc
    message("Region changed to: ", newLoc, " - reloading annotation")
    load_annotation(newLoc)
  })
}

shinyApp(ui, server)
