### Final project: BF-591-R
### Author: Italo Duran
### Email: duran01@bu.edu

#if (!requireNamespace("BiocManager", quietly= TRUE))
#install.packages("BiocManager")
#library(BiocManager)
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(colourpicker)# you might need to install this.
library(DT)# you might need to install this.
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)# you might need to install this.
#BiocManager::install("DESeq2")
library(DESeq2)# you might need to install this.
library(igraph)# you might need to install this.
#library(bslib)


###############################################################
#####Visual interface - Front End Part#########################
ui = fluidPage(theme = shinythemes::shinytheme("journal"),
  titlePanel("Final Project BF591-R"),
  h3("By: Italo Duran"),
  h4(HTML("<b>Visualization analysis of:</b><br> mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for<br>Huntington's Disease and neurologically normal individuals data set...")),
  h4(HTML("<b>Citation:</b> Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC et al.<br>RNA Sequence Analysis of Human Huntington Disease Brain Reveals an<br>Extensive Increase in Inflammatory and Developmental Gene Expression.<br>PLoS One 2015;10(12):e0143563. PMID: 26636579")),
  tags$h4(HTML("To use this application, download the <b>.csv files</b> from the data folder in my github:")),
  tags$a(href="https://github.com/imd9/Final-Project-BF591-R",target="_blank",rel="noopener noreferrer","Click here for 'My Github' CSV sample files!"), br(),br(),
  #This part is where we build the page layout and how it's going to look
  tabsetPanel(
    tabPanel("Sample", fluid = TRUE,
             h3("Sample Data Analysis"),
             HTML(paste(rep("In this section you will be able to visualize distinct values and distributions of the samples."), collapse = "")),br(),br(),
             sidebarLayout(sidebarPanel(fileInput("meta_data_samples", "Load sample Data", accept = ".csv", placeholder = "Sample_data.csv"),
                                        HTML(paste(rep("Download the <b>Sample1_data.csv</b> from the GitHub link above."), collapse = "")),br(),br(),
                                        HTML(paste(rep("<b>Summary:</b> It shows the different data types, Identifier values and Mean standard deviation.<br>
                                                        <b>Table:</b> Shows the samples that can be sorted by selecting the header names, chose the number of samples entries,
                                                        and use the search bar to look for specific values from the table.<br>
                                                        <b>Plots of Samples:</b> Histogram of the normal samples vs. with Huntington's, with different types of parameters to chose from."), collapse = ""))),
                           mainPanel(tabsetPanel(tabPanel("Summary", fluid = TRUE,tableOutput(outputId = 'sample_summary')),
                                                 tabPanel("Data", fluid = TRUE,dataTableOutput(outputId = "samples_mdata")),
                                                 tabPanel("Plots",plotOutput(outputId = 'p.m.i.'),plotOutput(outputId = 'r.i.n.'),plotOutput(outputId = 'a.o.d.')))))
             ),
    ################################################################################################################################################################################################################
    tabPanel("Counts",fluid = TRUE,
             h3("Counts Matrix Data Analysis"),
             HTML(paste(rep("In this section you will be able to visualize the normalized count data."), collapse = "")),br(),br(),
             sidebarLayout(sidebarPanel(fileInput(inputId = 'sample_counts_file', "Load sample Data", accept = ".csv", placeholder = "norm_counts.csv"),
                                        HTML(paste(rep("Download the <b>norm_counts.csv</b> from the GitHub link above."))),br(),br(),
                                        HTML(paste(rep("<b>Summary:</b> It shows the different data types, Identifier values and Mean standard deviation.<br>
                                                        <b>Diagnostic plot_sampscatter plots:</b> Shows the samples that can be sorted by selecting the header names, chose the number of samples entries,
                                                        and use the search bar to look for specific values from the table.<br>
                                                        <b>Plots of Samples:</b> Histogram of the normal samples vs. with Huntington's, with different types of parameters to chose from."))),br(),br(),
                                        sliderInput(inputId = 'percent_varians', label = 'Select a threshold for percent variance:',min=0,max=100,value=50,step=1),
                                        sliderInput(inputId = 'sample_non_zero', label = 'Select a threshold for number of genes that are non-zero:',min=0,max=100,value=50,step=1),
                                        submitButton("Plot", icon("r-project",class="fab fa-r-project fa-1x"), width = '100%')),
                           mainPanel(tabsetPanel(tabPanel("Summary",tableOutput(outputId='sample_normalz_summry')),
                                                 tabPanel("Diagnostic plot_sampscatter Plots",plotOutput(outputId='sample_median_varians'),plotOutput(outputId='samp_median_zero')),
                                                 tabPanel("Clustered Heatmap",plotOutput(outputId = "samp_heat_map")),
                                                 tabPanel("P.C.A.",selectInput(inputId="var_pca1",label="Select a PCA for the X-axis:",choices=c("PC1","PC2","PC3","PC4","PC5","PC5","PC6","PC7","PC8","PC9","PC10"),selected="PC1"),
                                                          selectInput(inputId="var_pca2",label="Select a PCA for the Y-axis:",choices=c("PC1","PC2","PC3","PC4","PC5","PC5","PC6","PC7","PC8","PC9","PC10"),selected="PC10"),
                                                          plotOutput(outputId="PCAplot")))))           
              ),
    #############################################################################################################################
    tabPanel("Differential Expression",fluid = TRUE,
             h3("Differential Expression Analysis"),
             HTML(paste(rep("In this section you will be able to visualize differential expression data."),collapse = "")),br(),br(),
             sidebarLayout(sidebarPanel(fileInput(inputId ='file1_de',"Load differential expression results",accept = ".csv",placeholder="deseq_res.csv"),
                                        HTML(paste(rep("Download the <b>ndeseq_diff_exp_res.csv</b> from the GitHub link above."))),br(),br(),
                                        HTML(paste(rep("<p>A volcano plot can be generated with <b>'log<sub>2</sub> fold-change'</b> on the x-axis and <b>'p-adjusted'</b> on the y-axis.</p>"),collapse="")),
                                        radioButtons(inputId='xaxis',label='Choose the column for the x-axis',choices=c('baseMean','log2FoldChange', 'lfcSE','stat','pvalue','padj'),selected = 'log2FoldChange'),
                                        radioButtons(inputId='yaxis',label='Choose the column for the x-axis',choices=c('baseMean','log2FoldChange', 'lfcSE','stat','pvalue','padj'),selected = 'padj'),
                                        colourpicker::colourInput(inputId='base_de_color',label="Base point color" ,"#138086"),
                                        colourpicker::colourInput(inputId='high_de_color',label="Highlight point color","#EEB462"),
                                        sliderInput(inputId = 'samp_padjad',min=-300,max=0,label="Select the magnitude of the p adjusted coloring:",value=-16,step=1),
                                        submitButton("Plot", icon("r-project",class="fab fa-r-project fa-1x"), width = '100%')),
                           mainPanel(tabsetPanel(tabPanel('Data Table',dataTableOutput("samp_sumr_dattable")),
                                                 tabPanel('Filtered Table',dataTableOutput("filtered_de_table")),
                                                 tabPanel('Volcano Plot',plotOutput("volcano")))))
              ),
    #####################################################################################################################
    tabPanel("G.S.E.A.",fluid = TRUE,
             h3("Gene Set Enrichment Analysis"),
             HTML(paste(rep("In this section you will be able to visualize<br>Gene Set Enrichment Analysisdata for each gene."),collapse = "")),br(),br(),
             sidebarLayout(sidebarPanel(HTML(paste(rep("Download the <b>Sample1_data.csv</b> from the GitHub link above."), collapse = "")),br(),
                                        fileInput(inputId = 'samp_gsea_data', label = 'Load sample information matrix CSV',accept = ".csv",placeholder="Sample1_data.csv"),
                                        HTML(paste(rep("Download the <b>norm_counts.csv</b> from the GitHub link above."))),br(),
                                        fileInput(inputId = 'samp_count_gsea', label = 'Load normalized counts matrix',accept = ".csv",placeholder="norm_counts.csv"),
                                        selectInput("metachoice", choices = c("PMI"="PMI","RIN"="RIN","Diagnosis"="Diagnosis","Sequence Reads"="Seq_reads","Age of death"="Age_of_death"),label="Select Data Parameter",selected="PMI"),
                                        HTML(paste(rep("E.g. of gene names:<br><b>ENSG00000069011.10, ENSG00000170689.8<br>ENSG00000180818.4, ENSG00000128710.5</b>"))),br(),br(),br(),
                                        textInput("gene", label = "Search for gene:", placeholder = "ENSG00000069011.10"),
                                        submitButton("Plot", icon("r-project",class="fab fa-r-project fa-1x"), width = '100%')),mainPanel(plotOutput("distroplot"))))
    
    
  )
)
##################################server##############################################################################
server = function(input, output, session){
  options(shiny.maxRequestSize=100*1024^2)  #increase file upload limit
  
  #########-SAMPLE-#############################################################
  load_meta = reactive({
    if (!is.null(input$meta_data_samples)){
      sample_metad = read_csv(input$meta_data_samples$datapath)
      return(sample_metad)}
    else{return(NULL)} })
  #function to produce table:
   samp_sumartbl = function(samp_metdat_tb){
    if (!is.null(input$meta_data_samples)){
       totsamp_sumrytbl = tibble("Name of Columns" = colnames(samp_metdat_tb), "Type" = sapply(samp_metdat_tb, class),
                         "Mean" = sapply(samp_metdat_tb, mean, na.rm = TRUE), "Standard Deviation " = sapply(samp_metdat_tb, sd, na.rm = TRUE))
      return( totsamp_sumrytbl)}
    else{return(NULL)} }
  ### PMI histogram:
  samp_pmi_hist = function(samp_metdat_tb){
    if (!is.null(input$meta_data_samples)){
      samp_histogr = ggplot(samp_metdat_tb, aes(PMI))+geom_histogram(bins = 10, color = "black", fill = "#138086")+
        labs(title="PMI Histogram")+xlab("PMI")+ylab("Count")+theme_classic()
      return(samp_histogr)}
    else{return(NULL)} }
  ### RIN histogram:
  samp_rin_hist = function(samp_metdat_tb){
    if (!is.null(input$meta_data_samples)){
      samp_histogr = ggplot(samp_metdat_tb,aes(RIN))+geom_histogram(bins=10,color="black",fill ="#EEB462")+
        labs(title='Histogram of RIN')+xlab('RIN')+ylab('Count')+theme_classic()
      return(samp_histogr)}
    else{return(NULL)} }
  ### AOD histogram:
  samp_aod_hist = function(samp_metdat_tb){
    if (!is.null(input$meta_data_samples)){
      samp_histogr = ggplot(samp_metdat_tb,aes(Age_of_death))+geom_histogram(bins = 10,color="black",fill="pink")+
        labs(title='Histogram of Age of Death')+xlab('Age of Death')+ylab('Count')+theme_classic()
      return(samp_histogr)}
    else{return(NULL)} }
  
  ######### SAMPLE-OUT-PUT #####################################
  output$sample_summary = renderTable({ samp_sumartbl(load_meta())})
  output$samples_mdata = DT::renderDataTable({load_meta()})
  output$a.o.d. = renderPlot({samp_aod_hist(load_meta())})
  output$r.i.n. = renderPlot({samp_rin_hist(load_meta())})
  output$p.m.i. = renderPlot({samp_pmi_hist(load_meta())})
  #_______________________________________________________________________________________________________________________________#
  ###-COUNTS-#################################################################################################  
  load_counts = reactive({
    if (!is.null(input$sample_counts_file)){
      counts <- read_csv(input$sample_counts_file$datapath)
      return(counts)}
    else{return(NULL)} })
  ### Summary table:
  count_sumry_tbl = function(counts_tibs, samp_percvarri, samp_nztib){
    if (!is.null(input$sample_counts_file)){
      samp_totvals = ncol(counts_tibs)-1
      samp_genestotls = nrow(counts_tibs)
      counts_tibs = counts_tibs %>% mutate(variance = apply(counts_tibs[-1], MARGIN = 1, FUN = var))
      sampperc_val = quantile(counts_tibs$variance, probs = samp_percvarri/100) 
      counts_tibs = filter(counts_tibs, variance >= sampperc_val)  
      counts_tibs = na_if(counts_tibs, 0)
      counts_tibs$non_zero = samp_totvals-rowSums(is.na(counts_tibs))  
      counts_tibs = filter(counts_tibs, non_zero >= samp_nztib)  
      samp_filtrgns = nrow(counts_tibs)    
      samp_pass_percgenes = samp_filtrgns/samp_genestotls*100
      sampfailgns = samp_genestotls-samp_filtrgns
      sampperc_failgns = sampfailgns/samp_genestotls*100
      #printing out summary results table: 
       totsamp_sumrytbl = tibble('Measured Parameters:' = c('Total Samples:', 'Total Genes:', 'Genes Passing:', "Genes Passing (%):", 'Genes Not Passing:', 'Genes Not Passing(%):'),
                         'Results:' = c(samp_totvals, samp_genestotls, samp_filtrgns, samp_pass_percgenes, sampfailgns, sampperc_failgns))
      return( totsamp_sumrytbl)}
    else{return(NULL)} }
  #### MEDIAN AND VARIANCE:
  samp_medvar = function(counts_tibs, samp_percvarri){
    if (!is.null(input$sample_counts_file)){
      ### TABLE DATA:
      samp_graphtib = counts_tibs%>%
        mutate(Median = apply(counts_tibs[-1], MARGIN = 1, FUN = median),Variance = apply(counts_tibs[-1], MARGIN = 1, FUN = var))
      sampperc_val = quantile(samp_graphtib$Variance, probs = samp_percvarri/100) 
      samp_graphtib = samp_graphtib %>% mutate(threshold = case_when(Variance >= sampperc_val ~ "TRUE", TRUE ~ "FALSE"))
      ### plot_sampscatter PLOT:
      cols = c("TRUE"="#138086","FALSE"="#EEB462")
      plot_sampscatter = ggplot(samp_graphtib, aes(Median, Variance))+geom_point(aes(color=threshold),alpha=0.75)+scale_color_manual(values = cols)+
        labs(title="Variance VS. Median:")+scale_y_log10()+scale_x_log10()+theme_classic()+theme(legend.position = 'right')
      return(plot_sampscatter)}
    else{return(NULL)} }
  ##### MEDIAN AND VARIANCE PLOT:
  smp_medvsnz <- function(counts_tibs, samp_nztib){
    if (!is.null(input$sample_counts_file)){
      samp_totvals = ncol(counts_tibs)-1  #store original number of samples and genes
      #make a plot tibble
      samp_graphtib = counts_tibs %>%   
        mutate(Median = apply(counts_tibs[-1], MARGIN = 1, FUN = median)) %>% na_if(0)  #calc median, convert 0 to NA
      samp_graphtib$no_zeros = rowSums(is.na(samp_graphtib))  #make new col, with counts.
      samp_graphtib = samp_graphtib %>% mutate(threshold = case_when(no_zeros <= samp_nztib ~ "TRUE", TRUE ~ "FALSE"))
      #plot plot_sampscatter plot
      cols = c("TRUE"="#138086","FALSE"="#EEB462")
      plot_sampscatter = ggplot(samp_graphtib, aes(Median,no_zeros))+geom_point(aes(color=threshold),alpha=0.75)+scale_color_manual(values = cols)+
        scale_x_log10()+labs(title="Non-Zero VS. Median:")+theme_classic()+ylab('Non-Zero Samples')+theme(legend.position='right')
      return(plot_sampscatter)}
    else{return(NULL)} }
  ### HEATMAP:
  plot_heatmap = function(counts_tibs, samp_percvarri){
    if (!is.null(input$sample_counts_file)){
      counts_tibs = log10(counts_tibs[-1])
      samp_graphtib = counts_tibs %>% 
        mutate(variance = apply(counts_tibs, MARGIN = 1, FUN = var))
      sampperc_val = quantile(samp_graphtib$variance, probs = samp_percvarri/100, na.rm = TRUE) 
      samp_graphtib = filter(samp_graphtib, variance >= sampperc_val)
      samp_heat_map = heatmap(as.matrix(samp_graphtib[-ncol(samp_graphtib)]), scale = "row")
      return(samp_heat_map)}
    else{return(NULL)} }
  #### PCA PLOTS:
  plot_pca = function(counts_tibs, samp_percvarri, var_pca1, var_pca2){
    if (!is.null(input$sample_counts_file)){
      filt_tib = counts_tibs %>% 
        mutate(variance = apply(counts_tibs[-1], MARGIN = 1, FUN = var), .after = gene)
      sampperc_val = quantile(filt_tib$variance, probs = samp_percvarri/100, na.rm = TRUE)   #calculate percentile
      filt_tib = filter(filt_tib, variance >= sampperc_val) #filter the tibble
      pca_res = prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) #transpose the data and perform PCA
      variance = summary(pca_res)$importance[2,]
      x = round(variance[var_pca1]*100, 2)
      y = round(variance[var_pca2]*100, 2)
      ### PRINTING PCA GRAPH:
      samp_graphtib = tibble(PC1 = pca_res$x[,var_pca1], PC2=pca_res$x[,var_pca2])
      pca = ggplot(samp_graphtib, aes(PC1, PC2))+geom_point(color="#138086")+labs(title="Princple Component Analysis Visualization:")+xlab(str_c(var_pca1, x,"% Variance",sep=" "))+
                    ylab(str_c(var_pca2,y,"% Variance",sep=" "))+theme_classic()
      return(pca)}
    else{return(NULL)} }
  
  ####### COUNTS-OUT-PUT #####################
  output$sample_normalz_summry = renderTable({count_sumry_tbl(load_counts(), input$percent_varians, input$sample_non_zero)})
  output$sample_median_varians = renderPlot({samp_medvar(load_counts(), input$percent_varians)})
  output$samp_median_zero = renderPlot({smp_medvsnz(load_counts(), input$sample_non_zero)})
  output$samp_heat_map = renderPlot({plot_heatmap(load_counts(), input$percent_varians)})
  output$PCAplot = renderPlot({plot_pca(load_counts(), input$percent_varians, input$var_pca1, input$var_pca2)})
  #_______________________________________________________________________________________________________________________________#
  ###### Differential Expression Analysis #####################################:
  load_de = reactive({
    if (!is.null(input$file1_de)){
      defp = read_csv(input$file1_de$datapath)
      return(defp)}
    else{return(NULL)}  })
  ### SEETING UP VOLCANO PLOT FUNCTION:---------
  volcano_plot = function(dataf, x_name, y_name, slider, color1, color2) {
    if (!is.null(input$file1_de)){
      cols = c("FALSE" = color1, "TRUE" = color2)
      threshold = 10^slider
      plot_data = dataf %>%
        dplyr::mutate(threshold = case_when(padj <= threshold ~ "TRUE", TRUE ~ "FALSE"))
      ### VOLCANO PLOT: (JUST LIKE IN HW7)---
      volc_plot=ggplot(plot_data,aes(x=!!sym(x_name),y=-log10(!!sym(y_name))))+geom_point(aes(color=threshold))+
        labs(color=str_glue('{y_name} 1 x 10^ {slider}'))+scale_color_manual(values=cols)+theme_classic()+
        theme(plot.margin = margin(2, 2, 2, 2, "cm"),legend.position="bottom")
      return(volc_plot)}
    else{return(NULL)} }
  ### HERE WE SET UP AN EXTRA FUNCTION TO SHOW US THE METADATA:
  draw_table = function(dataf, slider) {
    if (!is.null(input$file1_de)){
      df = dataf %>% 
        dplyr::filter(padj <= 10^slider)
      df[ ] = lapply(df,formatC,format="f",digits = 4)
      return(df)}
    else{return(NULL)} }
  ################ DE-OUT-PUT#########################
  output$samp_sumr_dattable = DT::renderDataTable({load_de()})
  output$volcano = renderPlot({volcano_plot(load_de(), input$xaxis, input$yaxis, input$samp_padjad, input$base_de_color, input$high_de_color)})
  output$filtered_de_table = renderDataTable({draw_table(load_de(), input$samp_padjad)})
  #_______________________________________________________________________________________________________________________________#
  ###### G.S.E.A.- Gene Set Enrichment Analysis ###############################:
  gsea_meta = reactive({
    if (!is.null(input$samp_gsea_data)){
      sample_metad = read_csv(input$samp_gsea_data$datapath)
      return(sample_metad)}
    else{return(NULL)} })
  gsea_datcounts = reactive({
    if (!is.null(input$samp_count_gsea)){
      counts = read_csv(input$samp_count_gsea$datapath)
      return(counts)}
    else{return(NULL)} })
  #function to make distribution plots-
  plot_distro = function(counts_tibs, samp_metdat_tb, meta_cat, selectgene){
    if (!is.null(input$samp_gsea_data) & !is.null(input$samp_count_gsea)){
      counts_tibs = column_to_rownames(counts_tibs, var = "gene")
      gene_counts = as.numeric(as.vector(counts_tibs[selectgene,]))
      samp_graphtib = tibble(Gene_Counts = gene_counts, meta_value = pull(samp_metdat_tb, meta_cat))
      if (meta_cat == "Diagnosis"){
        plot = ggplot(samp_graphtib, aes(meta_value))+
          geom_bar()+theme_classic()+labs(title = "Plot of gene counts vs Diagnosis")
        return(plot) }
      else {
        plot = ggplot(samp_graphtib, aes(meta_value,Gene_Counts))+geom_point()+theme_classic()+
                labs(title = str_c("Plot of gene counts vs ", meta_cat))
        return(plot) }}
    else{return(NULL)} }
  #######################G.S.E.A.-OUT-PUT#################################
  output$distroplot=renderPlot({plot_distro(gsea_datcounts(), gsea_meta(), input$metachoice, input$gene)})
  #_______________________________________________________________________________________________________________________________#
  #### LAST BRACKET FOR SERVER:
}
##### RUN THE APP-----------------
shinyApp(ui=ui, server = server)