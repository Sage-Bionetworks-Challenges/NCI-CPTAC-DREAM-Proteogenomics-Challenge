# server.R - Dashboard
library(shiny)
library(DMwR)
library(grid)
library(pheatmap)
library(ggplot2)
library(ggrepel)
# setwd("~/Documents/RWTH_Aachen/DREAM_CPTAC/ProteoExplorer/ProteoExplorer")
source("./functions.R")
load("./data/Testing.Rdata")
load("./data/Prediction.Rdata")
load("./data/Prediction_performance_per_gene.Rdata")

options(shiny.maxRequestSize=30*1024^2)

shinyServer(
  function(input, output, session) {
    
    F1_layer_X_0 = reactive({
      switch(input$F1_layer_X,
             # "mRNA (breast)" = Testing$prospective_breast_RNA_sort_common_gene_15107,
             # "True Proteomics (breast)" = Testing$prospective_breast_proteome_sort_common_gene_10005,
             # "mRNA (ovarian)" = Testing$prospective_ova_rna_seq_sort_common_gene_15121,
             "True Proteomics (ovarian)" = Testing$prospective_ova_proteome_sort_common_gene_7061)
    })
    F1_layer_X = reactive({ F1_layer_X <- convert_ID(input$select_Identifier,F1_layer_X_0() ) })
    
    F1_layer_Y_0 = reactive({
      switch(input$F1_layer_Y,
             "Team Guan, Predicted Proteomics (ovarian)" = Prediction$team_Guan_ovarian,
             "Ensemble, Predicted Proteomics (ovarian)"  = Prediction$ensemble_top4_ovarian)
    })
    F1_layer_Y = reactive({ F1_layer_Y <- convert_ID(input$select_Identifier,F1_layer_Y_0() ) })
    
    F1_layer_X_common = reactive({ 
      F1_layer_X()[ intersect(rownames(F1_layer_X()), rownames(F1_layer_Y())) , ] 
    })
    
    F2_1_layer_0 = reactive({
      switch(input$F2_1_layer,
             "Team Guan (breast)"  = Prediction_performance_per_gene$team_Guan_breast,
             "Team Guan (ovarian)" = Prediction_performance_per_gene$team_Guan_ovarian )
    })
    F2_1_layer = reactive({ F1_2_layer <- convert_ID_col(input$select_Identifier,F2_1_layer_0() ) })
    
    gene_density  = reactive({ input$F2_1_gene })
    
    # Display Omics layer
    output$F1_layer_X_common = DT::renderDataTable(
      F1_layer_X_common(),
      options = list(scrollX = TRUE)
    )
    
    output$F2_1_layer = DT::renderDataTable(
      F2_1_layer(),
      options = list(scrollX = TRUE)
    )
    
    gene_Scatter_plot = reactive({ input$F1_1_gene })
    
    
    # Plot scatter plot
    output$scatter_plot = renderPlot({
      withProgress(message="Draw Scatter plot",value=1,{
        if(gene_Scatter_plot() %not in% rownames(F1_layer_X_common())) {gene <- rownames(F1_layer_X_common())[1]} else {gene <- gene_Scatter_plot()}
        
        df <- cbind( as.numeric(F1_layer_X()[gene, ]) , as.numeric(F1_layer_Y()[gene, ]) ) 
        colnames(df) <- c("a","b") ; df <- data.frame(df)  ; df <- df[complete.cases(df), ]
        scatter_plot(gene=gene,real_gene=gene_Scatter_plot(),df=df,label_X=input$F1_layer_X,label_Y=input$F1_layer_Y,text_size=22,title_size=2 ) 
        
      })
    })
    
    # switch from update to results menu as soon as the go button is pushed
    observeEvent(input$go_scatter_plot, {
      new.menu = switch(input$menu,"Data_Navigation"="Results")
      updateTabItems(session,"menu",new.menu)
    })
    
    # Plot Density
    output$density_plot = renderPlot({
      withProgress(message="Draw Density plot",value=1,{
        df <- as.data.frame(F2_1_layer())
        density_plot(df=df,gene=gene_density())
      })
    })
    
    # switch from update to results menu as soon as the go button is pushed
    observeEvent(input$go_density, {
      new.menu = switch(input$menu,"Data_Navigation"="Results")
      updateTabItems(session,"menu",new.menu)
    })
    
    outputOptions(output,"scatter_plot",suspendWhenHidden=F)
    outputOptions(output,"density_plot",suspendWhenHidden=F)
  }
)

