library(shiny)
library(ggplot2)
library(plotly)
library(data.table)
library(dplyr)
library(gtools)
library(ggpubr)
library(ggrepel)
library(reshape)

ui <- fluidPage(
  titlePanel("Expression Visualization App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file_upload", "Upload a file (txt or tab-separated, log2 scale is going to be shown so please provide untransformed counts)"),
      uiOutput("select_columns_ui"),
      textInput("groups_names", "Enter comma-separated (no spaces) group names for each one of the selected columns", value = ""),
      textInput("gene_name", "Enter Gene ID (just one)", value = ""),
      actionButton("plot_button", "Generate Plot"),
      downloadButton("download_plot_png", "Download PNG"),
      downloadButton("download_plot_pdf", "Download PDF"),
      verbatimTextOutput("text_content"),
    ),
    mainPanel(
      tabsetPanel(
        id = 'my_tabsetPanel',
        type = 'tabs',
        tabPanel("Plot", plotlyOutput("scatter_plot")),
        tabPanel("Plot_ttest", plotOutput("scatter_plot_ttest")),
        tabPanel("Plot_wilcoxtest", plotOutput("scatter_plot_wtest")),
        tabPanel("Plot_kruskaltest", plotOutput("scatter_plot_ktest")),
        tabPanel("Plot_anova", plotOutput("scatter_plot_anova"))
      )
    )
  )
  
)

server <- function(input, output) {
  
  # Reactive expression to read the uploaded file
  data <- reactive({
    req(input$file_upload)
    fread(input$file_upload$datapath, header = TRUE, sep = "\t")
  })
  
  # Dynamically generate the UI for column selection based on the uploaded file
  output$select_columns_ui <- renderUI({
    req(data())
    selectInput(inputId = "columns", label = "Select the columns of each group to be included", choices = colnames(data()), multiple = TRUE, selectize = TRUE)
  })
  
  plot_output <- reactive({
    # Create ggplot object
    gene_counts_rpkm_to_plot <- data()
    if (input$gene_name != "" & input$gene_name %in% gene_counts_rpkm_to_plot$Gene_ID){
      gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot %>% filter(Gene_ID == input$gene_name)
    }
    
    a <- as.data.frame(gene_counts_rpkm_to_plot)[,mixedsort(input$columns)]
    df <- data.frame(sample=names(a),
                     expr_RPKM=as.numeric(unname(a)),
                     condition=mixedsort(unlist(strsplit(input$groups_names,","))))
    df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
    
    df2 <- reshape::melt(df)
    df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])        
    df2 <- df2 %>%
      add_count(condition, name = "condition_n")
    df2$sample2 <- as.numeric(1:dim(df2)[1])
    list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
    df2$i <- 1
    labs_total <- aggregate(i~condition,df2,sum)
    labs_total$value <- NA
    labs_total$i <- paste0("Total points: ",labs_total$i)
    p <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 1,label=i)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0("GENE SHOWN: ", input$gene_name))
    return(p)
    
  })
  plot_output_ttest <- reactive({
    # Create ggplot object
    gene_counts_rpkm_to_plot <- data()
    if (input$gene_name != "" & input$gene_name %in% gene_counts_rpkm_to_plot$Gene_ID){
      gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot %>% filter(Gene_ID == input$gene_name)
    }
    
    a <- as.data.frame(gene_counts_rpkm_to_plot)[,mixedsort(input$columns)]
    df <- data.frame(sample=names(a),
                     expr_RPKM=as.numeric(unname(a)),
                     condition=mixedsort(unlist(strsplit(input$groups_names,","))))
    df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
    
    df2 <- reshape::melt(df)
    df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])        
    df2 <- df2 %>%
      add_count(condition, name = "condition_n")
    df2$sample2 <- as.numeric(1:dim(df2)[1])
    list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
    df2$i <- 1
    labs_total <- aggregate(i~condition,df2,sum)
    labs_total$value <- NA
    labs_total$i <- paste0("Total points: ",labs_total$i)
    q <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      stat_compare_means(method = "t.test", comparisons = list_combinations) +
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 1,label=i)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0("GENE SHOWN: ", input$gene_name))
    return(q)
    
  })
  plot_output_wtest <- reactive({
    # Create ggplot object
    gene_counts_rpkm_to_plot <- data()
    if (input$gene_name != "" & input$gene_name %in% gene_counts_rpkm_to_plot$Gene_ID){
      gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot %>% filter(Gene_ID == input$gene_name)
    }
    
    a <- as.data.frame(gene_counts_rpkm_to_plot)[,mixedsort(input$columns)]
    df <- data.frame(sample=names(a),
                     expr_RPKM=as.numeric(unname(a)),
                     condition=mixedsort(unlist(strsplit(input$groups_names,","))))
    df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
    
    df2 <- reshape::melt(df)
    df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])        
    df2 <- df2 %>%
      add_count(condition, name = "condition_n")
    df2$sample2 <- as.numeric(1:dim(df2)[1])
    list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
    df2$i <- 1
    labs_total <- aggregate(i~condition,df2,sum)
    labs_total$value <- NA
    labs_total$i <- paste0("Total points: ",labs_total$i)
    r <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      stat_compare_means(method = "wilcox.test", comparisons = list_combinations) +
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 1,label=i)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0("GENE SHOWN: ", input$gene_name))
    return(r)
    
  })
  plot_output_ktest <- reactive({
    # Create ggplot object
    gene_counts_rpkm_to_plot <- data()
    if (input$gene_name != "" & input$gene_name %in% gene_counts_rpkm_to_plot$Gene_ID){
      gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot %>% filter(Gene_ID == input$gene_name)
    }
    
    a <- as.data.frame(gene_counts_rpkm_to_plot)[,mixedsort(input$columns)]
    df <- data.frame(sample=names(a),
                     expr_RPKM=as.numeric(unname(a)),
                     condition=mixedsort(unlist(strsplit(input$groups_names,","))))
    df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
    
    df2 <- reshape::melt(df)
    df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])        
    df2 <- df2 %>%
      add_count(condition, name = "condition_n")
    df2$sample2 <- as.numeric(1:dim(df2)[1])
    list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
    df2$i <- 1
    labs_total <- aggregate(i~condition,df2,sum)
    labs_total$value <- NA
    labs_total$i <- paste0("Total points: ",labs_total$i)
    s <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      stat_compare_means(method = "kruskal.test", comparisons = list_combinations) +
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 1,label=i)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0("GENE SHOWN: ", input$gene_name))
    return(s)
    
  })
  plot_output_anova <- reactive({
    # Create ggplot object
    gene_counts_rpkm_to_plot <- data()
    if (input$gene_name != "" & input$gene_name %in% gene_counts_rpkm_to_plot$Gene_ID){
      gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot %>% filter(Gene_ID == input$gene_name)
    }
    
    a <- as.data.frame(gene_counts_rpkm_to_plot)[,mixedsort(input$columns)]
    df <- data.frame(sample=names(a),
                     expr_RPKM=as.numeric(unname(a)),
                     condition=mixedsort(unlist(strsplit(input$groups_names,","))))
    df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)
    
    df2 <- reshape::melt(df)
    df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])        
    df2 <- df2 %>%
      add_count(condition, name = "condition_n")
    df2$sample2 <- as.numeric(1:dim(df2)[1])
    list_combinations <- strsplit(unique(unlist(lapply(strsplit(apply(expand.grid(unique(df2$condition), unique(df2$condition)),1,function(x){paste(x,collapse="*****")}),"*****",fixed=T),function(x){if (x[1] != x[2]){paste(sort(x),collapse="*****")}}))),"*****",fixed=T) # Used ***** only to make sure is a separator that's not going to be used in the conditions name
    df2$i <- 1
    labs_total <- aggregate(i~condition,df2,sum)
    labs_total$value <- NA
    labs_total$i <- paste0("Total points: ",labs_total$i)
    t <- ggplot(df2, aes(x=condition, y=value,color=condition,label = sample2)) +
      geom_violin(trim=T) +
      stat_compare_means(method = "anova", comparisons = list_combinations) +
      geom_boxplot(width=0.1) +
      geom_jitter(shape=16, position=position_jitter(0.2,seed = 1)) +
      ggrepel::geom_label_repel(position = position_jitter(0.2,seed = 1)) +
      geom_point(data = dplyr::filter(df2, condition_n == 1)) +
      geom_text(data=labs_total,aes(x=condition,y= max(df2$value) + 1,label=i)) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=9)) +
      labs(x="",y="Expr_RPKM_log2_01",color="Condition", title=paste0("GENE SHOWN: ", input$gene_name))
    return(t)
    
  })

  observeEvent(input$plot_button, {
    output$text_content <- renderText({
      gene_counts_rpkm_to_plot <- data()
      if (input$gene_name != "" & input$gene_name %in% gene_counts_rpkm_to_plot$Gene_ID){
        gene_counts_rpkm_to_plot <- gene_counts_rpkm_to_plot %>% filter(Gene_ID == input$gene_name)
      }

      a <- as.data.frame(gene_counts_rpkm_to_plot)[,mixedsort(input$columns)]
      df <- data.frame(sample=names(a),
                       expr_RPKM=as.numeric(unname(a)),
                       condition=mixedsort(unlist(strsplit(input$groups_names,","))))
      df$"Expr_RPKM_log2_01" <- round(log2(df$expr_RPKM + 0.1),2)

      df2 <- reshape::melt(df)
      df2 <- unique(df2[df2$variable=="Expr_RPKM_log2_01",c("condition","value","sample")])
      df2 <- df2 %>%
        add_count(condition, name = "condition_n")
      df2$sample2 <- as.numeric(1:dim(df2)[1])
      lab_list <- paste0("Note: If the output of a statistical test is not shown, then it yielded some error and it's likely not the most appropriate one\nSamples_numbering:\n",paste0("Number_",1:length(df2$sample),": ",df2$sample,collapse="\n"))
    })
    
    output$scatter_plot <- renderPlotly({
        p <- plot_output()
        ggplotly(p) # Convert ggplot object to plotly object
    })
    
    output$scatter_plot_ttest <- renderPlot({
      q <- plot_output_ttest()
      q
    })
    
    output$scatter_plot_wtest <- renderPlot({
      r <- plot_output_ttest()
      r
    })
    
    output$scatter_plot_ktest <- renderPlot({
      s <- plot_output_ttest()
      s
    })
    
    output$scatter_plot_anova <- renderPlot({
      t <- plot_output_anova()
      t
    })
    
  })
  
  # Download the plots:
  # Download plot in PNG format
  output$download_plot_png <- downloadHandler(
    filename = function() {
      paste("scatter_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Create a ggplot object and save it to a file
      p <- plot_output() # Retrieve the ggplot object
      ggsave(file, plot = p, device = "png", width=30, height=30)
    },
    contentType = "image/png"
  )
  
  # Download plot in PDF format
  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste("scatter_plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      # Create a ggplot object and save it to a file
      p <- plot_output() # Retrieve the ggplot object
      ggsave(file, plot = p, device = "pdf", width = 30, height = 30)
    },
    contentType = "application/pdf"
  )
}
  
shinyApp(ui = ui, server = server)
# shiny::runApp("/home/pepeluisrr/Descargas/")
