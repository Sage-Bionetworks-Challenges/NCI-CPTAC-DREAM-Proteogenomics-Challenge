###################################################### PLOTTING ######################################################
library(ggplot2) # for plotting
library(ggsignif)
library(gridExtra) # for arranging plots
library(beeswarm)
library(ggrepel)

scatter_plot <- function(df,title,x_lab,y_lab,switch_anchor=F,text_size=33,title_size=2.5,coeff_x=0.5) {
  # equation, correlation and p value
  df <- as.data.frame(df) ; colnames(df) <- c("a","b")
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(~~italic("r")~"="~r*","~~italic("p")~"="~p,
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(unname(r), digits = 2),
                          p = format(p, digits=4)))
    as.character(as.expression(eq));                 
  }
  g <- ggplot(df, aes(a, b)) + 
    geom_point(shape = 16, size = 10, show.legend = FALSE, alpha = .7, color = "darkred" ) +  geom_smooth(method=lm,se=F,show.legend=F ) + 
    labs(x =x_lab, y=y_lab) + ggtitle(title) + 
    theme(legend.position="bottom",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size=rel(title_size), hjust=0.5),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = min(df$a) + coeff_x*(max(df$a)-min(df$a)), y = min(df$b) + 0*(max(df$b)-min(df$b)), label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 18 )  
  g
}

nice_boxplot <- function (result, title, text_size, title_size, Y_label="Pearson correlation (r)" ,X_label="features") {
  result_transformed <- c()
  for (i in 1:length(colnames(result))) {
    v <- result[ ,i] ; v <- cbind(rep(i,each=length(rownames(result))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("features", "correlation")
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed$features <- factor(result_transformed$features,levels=1:length(colnames(result)), labels=colnames(result)) 
  
  give.n <- function(x){ return(c(y = median(x)*1.1, label = length(x)))  }
  mean.n <- function(x){ return(c(y = median(x)*0.85, label = round(mean(x),2))) }
  
  r <- ggplot(data = result_transformed, aes(y = correlation, x = features, fill = features )) + geom_boxplot(width = 0.35 , colour="black") + 
    geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black"))  + 
    geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange") ) + 
    stat_summary(fun.data = give.n, geom = "text", fun.y = median, colour = "white",size = 5) +
    stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "white",size = 5) + 
    # geom_signif( comparisons = list(c("GEX","progeny11")),test = "t.test",y_position=1) + 
    ggtitle(title) + ylab(Y_label) + xlab(X_label) + 
    theme(legend.position="none",plot.title=element_text(size=text_size+8,hjust=0.5),
          axis.text=element_text(colour="black",size=text_size),
          axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") )  + 
    scale_y_continuous( breaks=c(-0.4,0,0.4, 0.8), limits=c(-0.5, 1)) 
  grid.arrange(r, ncol = 1)
}


nice_boxplot_standard <- function (result, title, text_size, title_size, Y_label="Pearson correlation of observed vs predicted" ,X_label="features") {
  result_transformed <- c()
  for (i in 1:length(colnames(result))) {
    v <- result[ ,i] ; v <- cbind(rep(i,each=length(rownames(result))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("features", "correlation")
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed$features <- factor(result_transformed$features,levels=1:length(colnames(result)), labels=colnames(result)) 
  
  give.n <- function(x){ return(c(y = median(x)*1.1, label = length(x)))  }
  mean.n <- function(x){ return(c(y = median(x)*0.85, label = round(mean(x),2))) }
  
  r <- ggplot(data = result_transformed, aes(y = correlation, x = features, color = features)) + geom_boxplot( width = 0.3 , colour="black", fill="red") + 
    geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black"))  + 
    geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange") ) + 
    stat_summary(fun.data = give.n, geom = "text", fun.y = median, colour = "white",size = 5) +
    stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "white",size = 5) + 
    # geom_signif( comparisons = list(c("mRNA available","mRNA missing")),test = "wilcox.test",y_position=1,color="black" ) + 
    ggtitle(title) + ylab(Y_label) + xlab(X_label) + 
    theme(legend.position="none",plot.title=element_text(size=text_size+8,hjust=0.5),
          axis.text=element_text(colour="black",size=text_size),
          axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") )  + 
    scale_y_continuous( breaks=c(-0.4,0,0.4, 0.8), limits=c(-0.4, 1)) 
  grid.arrange(r, ncol = 1)
}

# nice_boxplot(result_new_cell, title, text_size=18, title_size = 1.8)

nice_boxplot_list_standard <- function (result, title, text_size=20, title_size, Y_label="Prediction performance (r)" ,X_label="features" ) {
  result_transformed <- c()
  for (i in 1:length(result)) {
    v <- na.omit(result[[i]]) ; v <- cbind(rep(i,each=length(na.omit(result[[i]]))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("features", "correlation")
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed$features <- factor(result_transformed$features,levels=1:length(names(result)), labels=names(result)) 
  
  give.n <- function(x){ return(c(y = median(x)*1.1, label = length(x)))  }
  mean.n <- function(x){ return(c(y = median(x)*0.85, label = round(mean(x),2))) }
  
  r <- ggplot(data = result_transformed, aes(y = correlation, x = features, color = features)) + geom_boxplot( width = 0.3 , colour="black", fill="red") + 
    geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black"))  + 
    geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange") ) + 
    stat_summary(fun.data = give.n, geom = "text", fun.y = median, colour = "white",size = 5) +
    stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "white",size = 5) + 
    # geom_signif( comparisons = list(c("mRNA available","mRNA missing")),test = "wilcox.test",y_position=1) + 
    ggtitle(title) + ylab(Y_label) + xlab(X_label) + 
    theme(legend.position="none",plot.title=element_text(size=text_size+8,hjust=0.5),
          axis.text=element_text(colour="black",size=text_size),
          axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") )  + 
    scale_y_continuous( breaks=c(-0.4,0,0.4, 0.8), limits=c(-0.4, 1)) 
  grid.arrange(r, ncol = 1)
}

nice_boxplot_list <- function (result, title, text_size, title_size, Y_label="Pearson correlation (r)" ,X_label="features") {
  result_transformed <- c()
  for (i in 1:length(result)) {
    v <- na.omit(result[[i]]) ; v <- cbind(rep(i,each=length(na.omit(result[[i]]))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("features", "correlation")
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed$features <- factor(result_transformed$features,levels=1:length(names(result)), labels=names(result)) 
  
  give.n <- function(x){ return(c(y = median(x)*1.1, label = length(x)))  }
  mean.n <- function(x){ return(c(y = median(x)*0.85, label = round(mean(x),2))) }
  
  r <- ggplot(data = result_transformed, aes(y = correlation, x = features, fill = features )) + geom_boxplot(width = 0.35 , colour="black") + 
    geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black"))  + 
    geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange") ) + 
    stat_summary(fun.data = give.n, geom = "text", fun.y = median, colour = "black",size = 6) +
    stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "black",size = 6) + 
    # geom_signif( comparisons = list(c("GEX","progeny11")),test = "t.test",map_signif_level=TRUE,y_position=1) + 
    ggtitle(title) + ylab(Y_label) + xlab(X_label) + 
    theme(legend.position="none",plot.title=element_text(size=text_size+8,hjust=0.5),
          axis.text=element_text(colour="black",size=text_size),
          axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") )  + 
    scale_y_continuous( breaks=c(-0.4,0,0.4, 0.8), limits=c(-0.5, 1)) 
  grid.arrange(r, ncol = 1)
}
# nice_boxplot_list(result_new_cell, title, text_size=15 )

nice_density_plot_list <- function (x, title , text_size=22 ) {
  library(ggplot2) 
  result <- x ;
  result_transformed <- c()
  for (i in 1:length(result)) {
    v <- na.omit(result[[i]]) ; v <- cbind(rep(i,each=length(na.omit(result[[i]]))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("groups", "correlation")
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed[,1] <- factor(result_transformed[,1],levels=1:length(result) , labels=names(x) )

  ggplot(result_transformed, aes(x = correlation, fill = groups)) + geom_density(alpha = 0.3 )  + 
    
    theme(legend.position="botttom",plot.title=element_text(size=text_size+7,hjust=0.5),
          axis.text=element_text(colour="black",size=text_size),
          axis.title.x=element_text(colour="black",size=text_size),axis.title.y=element_text(colour="black",size=text_size),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) 
  q  
}

nice_boxplot_flip <- function (result, title, Yaxis="correlation") {
  result_transformed <- c()
  for (i in 1:length(colnames(result))) {
    v <- result[ ,i] ; v <- cbind(rep(i,each=length(rownames(result))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("features", Yaxis)
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed$features <- factor(result_transformed$features,levels=1:length(colnames(result)), labels=colnames(result)) 
  r <- qplot(data = result_transformed, x = features, y = expression,  fill = features, geom = "boxplot", main=title) + aes(group=features)
  r <- r + geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
    geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + coord_flip()
  grid.arrange(r, ncol = 1)
}

#theme(legend.position="right",legend.text=element_text(size=7), axis.text=element_text(size=8)) + 

Violin_plot <- function (result, title) {
  result_transformed <- c()
  for (i in 1:length(colnames(result))) {
    v <- result[ ,i] ; v <- cbind(rep(i,each=length(rownames(result))), v)
    result_transformed <- rbind(result_transformed, v)
    colnames(result_transformed) <- c("features", "correlation")
  }
  result_transformed <- data.frame(result_transformed)
  result_transformed$features <- factor(result_transformed$features,levels=1:length(colnames(result)), labels=colnames(result)) 
  
  ggplot(result_transformed, aes(x = features, y = correlation, fill = features) ) +  
    geom_violin(trim = FALSE) + 
    theme(legend.position="bottom") + labs(x = "features", y = "correlation") + 
    geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
    geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + ggtitle(title)+ 
    theme(legend.position="bottom",axis.text=element_text(size=7)) 
}

##############   r + coord_flip()   ############### to flip the plot horizontally ###################################


############################## subset columns ##############################
subset_col_abs <- function(mat, abs_limit , obs) {
  value <- as.matrix(mat)
  value <- as.numeric(value)
#  value <- value[!is.na(value)]
  th <- quantile(abs(value), abs_limit )
  
  redoheatmap <- c() ;  
  for (i in 1:ncol(mat)) {
    if ( length( which( abs(mat[ ,i]) > th ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, colnames(mat)[i]) 
    }
  }
  r <- mat[ , redoheatmap ] 
  return(r)
}

subset_col_up <- function(mat, up_limit , obs) {
  value <- as.matrix(mat)
  value <- as.numeric(value)
#  value <- value[!is.na(value)]
  th <- quantile(value, up_limit )
  
  redoheatmap <- c() ;  
  for (i in 1:ncol(mat)) {
    if ( length( which( mat[ ,i] > th ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, colnames(mat)[i]) 
    }
  }
  r <- mat[ , redoheatmap ] 
  return(r)
}

subset_col_down <- function(mat, down_limit , obs) {
  value <- as.matrix(mat)
  value <- as.numeric(value)
#  value <- value[!is.na(value)]
  th <- quantile( value , down_limit )
  
  redoheatmap <- c() ;  
  for (i in 1:ncol(mat)) {
    if ( length( which(  mat[ ,i]  < th ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, colnames(mat)[i]) 
    }
  }
  r <- mat[ , redoheatmap ] 
  return(r)
}

############################## subset rows ##############################
subset_row_abs <- function(mat, abs_limit , obs) {
  value <- as.matrix(mat)
  value <- as.numeric(value)
#  value <- value[!is.na(value)]
  th <- quantile(abs(value), abs_limit )
  
  redoheatmap <- c() ;  
  for (i in 1:nrow(mat)) {
    if ( length( which( abs(mat[i,]) > th ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, rownames(mat)[i]) 
    }
  }
  r <- mat[ redoheatmap , ] 
  return(r)
}

subset_row_up <- function(mat, up_limit , obs) {
  value <- as.matrix(mat)
  value <- as.numeric(value)
#  value <- value[!is.na(value)]
  th <- quantile(value, up_limit )
  
  redoheatmap <- c() ;  
  for (i in 1:nrow(mat)) {
    if ( length( which( mat[ i , ] > th ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, rownames(mat)[i]) 
    }
  }
  r <- mat[ redoheatmap , ] 
  return(r)
}

subset_row_up_threshold <- function(mat, threshold , obs) {
  redoheatmap <- c() ;  
  for (i in 1:nrow(mat)) {
    if ( length( which( mat[ i , ] > threshold ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, rownames(mat)[i]) 
    }
  }
  r <- mat[ redoheatmap , ] 
  return(r)
}


subset_row_down_threshold <- function(mat, threshold , obs) {
  redoheatmap <- c() ;  
  for (i in 1:nrow(mat)) {
    if ( length( which( mat[ i , ] < threshold ) ) >= obs ) {
      redoheatmap <- c(redoheatmap, rownames(mat)[i]) 
    }
  }
  r <- mat[ redoheatmap , ] 
  return(r)
}

subset_row_variance <- function(mat, N) {
  v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:N] , ]
  return(mat)
}
subset_col_variance <- function(mat, N) {
  v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[  , names(v)[1:N] ]
  return(mat)
}



#################################### plot matrix ####################################

######## Label rotation pheatmap
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

library(pheatmap) ; library(grid)
plot_pheatmap <- function(mat, row_names, col_names , title ,cluster_rows=T,cluster_cols=T,fontsize=25,fontsize_row=25, fontsize_col=25, scale="none") {
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(mat, main=title, fontsize=fontsize, fontsize_row=fontsize_row,fontsize_col=fontsize_col,cluster_rows = cluster_rows, cluster_cols = cluster_cols,scale=scale,display_numbers=F)
  setHook("grid.newpage", NULL, "replace")
  grid.text(row_names, x=-0.03, rot=90, gp=gpar(fontsize=30))
  grid.text(col_names, y=0.01, gp=gpar(fontsize=30)  )
}


# library(devtools) ; install_github("plotflow", "trinker")
library(plotflow)

#  Basically you need columns as sample(or conditions) and rows as individuals(Proteins, genes, whatever..) that are measured in each condition

magicPlotMaker <- function(df, outpath){
  
  setwd(outpath)
  
  print("Dimension of dataframe :")
  print(dim(df))
  
  ##This part is just to generate the melted dataframe for ggplot
  df$ID <- row.names(df)
  
  melted_df <- melt(df)
  index <- c(1:length(melted_df[,1]))
  
  df <- df[,-length(df[1,])]
  
  #########################################
  ##                ggplots              ##
  #########################################
  
  a <- ggplot(melted_df, aes(x = value, fill = variable)) + geom_density(alpha = 0.3) + theme_minimal()
  b <- ggplot(melted_df, aes(x = index, y = value, group = variable, color = variable)) + geom_point() + theme_minimal()
  c <- ggplot(melted_df, aes(x = index, y = value, group = variable, color = variable)) + geom_boxplot() + theme_minimal()
  
  ggsave("densities.pdf", plot = a, device = pdf)
  ggsave("cloud_point.pdf", plot = b, device = pdf)
  ggsave("box_plot.pdf", plot = c, device = pdf)
  
  #########################################
  ##                PCA                  ##
  #########################################
  
  complete_df <- df[complete.cases(df),]
  print("Dimension of complete cases :")
  print(dim(complete_df))
  t_complete_df <- t(complete_df)
  PCA <- prcomp(t_complete_df ,center = TRUE, scale. = T) 
  
  pdf("boulder_plot.pdf")
  plot(PCA)
  dev.off()
  pdf("PCA.pdf")
  plot(PCA$x[,1],PCA$x[,2], pch = 19, xlab = paste("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"), ylab = paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"), main = "PCA of samples") 
  text(PCA$x[,1],PCA$x[,2], labels = names(PCA$x[,1]), pos = 3)
  dev.off()
  
  #########################################
  ##                heatmap              ##
  #########################################
  
  pheatmap(t_complete_df, show_colnames =  FALSE, filename = "profil_heatmap.pdf")
  pheatmap(cor(complete_df), filename = "correlation_heatmap.pdf")
}


#df = result from limma differnecial expression (with column X as identifier column)
#gtt = gene to term table
#consensus = consensus table of piano analisys (with column X as term column)
#outpath = path to the directory were the figures will be stored
#mt = mapping table between the identifiers of limma analisys and gene to term table
volcano_nice <- function(df,hAss,vAss,FCIndex,pValIndex){
  df <- df[complete.cases(df),]
  
  hAss <- -log(hAss)
  
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  
  xlimAbs <- ceiling(abs(max(df[,FCIndex])))
  ylimAbs <- ceiling(abs(max(-log(df[,pValIndex]))))
  
  xneg <- function(x) abs(hAss-1+x/(x+vAss))
  xpos <- function(x) abs(hAss-1+x/(x-vAss))
  
  test <- function(x,y,vAss){
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss){
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  
  
  df$couleur <- "0"
  df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),as.numeric(x[pValIndex]), vAss))
  df_005 <- subset(df, couleur == "1")
  
  if (length(df_005[,1]) > 0) 
  {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) + 
      geom_point(alpha = 0.5) + 
      #geom_text(aes(label = df_005$X, vjust = -1.5, nudge_y = 0.5, hjust = 0.15, check_overlap = TRUE)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs,-vAss), color = "black", alpha = 0.7) + 
      ylim(c(0,ylimAbs)) + xlim(c(-xlimAbs,xlimAbs)) + 
      stat_function(fun = xpos, xlim = c(vAss,xlimAbs), color = "black", alpha = 0.7) + 
      scale_colour_manual(values = c("grey30","royalblue3")) + 
      theme_minimal() +
      theme(legend.position = "none") 
  }
  else
  {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) + 
      geom_point(alpha = 0.5) + 
      stat_function(fun = xneg, xlim = c(-xlimAbs,-vAss), color = "black", alpha = 0.7) + 
      ylim(c(0,ylimAbs)) + xlim(c(-xlimAbs,xlimAbs)) + 
      stat_function(fun = xpos, xlim = c(vAss,xlimAbs), color = "black", alpha = 0.7) + 
      scale_colour_manual(values = c("grey30","royalblue3")) + 
      theme_minimal() +
      theme(legend.position = "none")
  }
  return(a)
}

volcano_genesets <- function(df, gtt, consensus, outpath, mt = NULL, get_df = 0) {
  gtt <- as.data.frame(apply(gtt, 2, function(x) toupper(x)))
  names(gtt) <- c("gene","term")
  colnames(df)[1] <- "X"
  colnames(consensus)[1] <- "X"
  
  if (!is.null(mt))
  {
    mt <- apply(mt, 2, function(x) toupper(x))
    
    names(mt) <- c("X","gene")
    
    gtt <- merge(gtt,mt, all = FALSE)
    
    names(gtt) <- c("gene","term","X")
  }
  else
  {
    names(gtt) <- c("X","term")
  }
  
  df$X <- toupper(df$X)
  df <- merge(df,gtt, all = FALSE)
  
  setwd(outpath)
  
  terms <- consensus$X
  
  for (term in terms) 
  {
    subDf <- df[df$term == term,]
    
    if (length(subDf[,1]) > 4) 
    {
      print(term)
      print(length(subDf[,1]))
      #normalized <- (subDf$AveExpr-min(subDf$AveExpr))/(max(subDf$AveExpr)-min(subDf$AveExpr))
      a <- volcano_nice(subDf,0.05,0.5,7,6)
      print("yo")
      ggsave(paste(term,".pdf"), a, device = pdf)
    }
  }
  if (get_df == 1){
    return(df)
  }
}

######################################################### LINE PLOT #########################################################

## baseline=corr_HGSC_PROT_RNA ; comparison=pred_HGSC_PROT_by_RNA[common_gene,1]
comparison_plot <- function(baseline,baseline_name,comparison, comparison_name,legend="features",xlab="ordered genes",ylab="expression", text_size= 20 , title_size = 1 )  {
  ###### baseline 
  baseline <- baseline[order(baseline)]
  target <- names(baseline)
  ###### comparison
  comparison <- comparison[match(target, names(comparison))]
  baseline <- as.vector(baseline) ;comparison <- as.vector(comparison)   
  tab <- as.data.frame(cbind(1:length(baseline), comparison, baseline) ); colnames(tab)[1] <- "index"
  rownames(tab) <- target
  
  tab$diffup <- tab[,2] - tab[,3] ; tab$diffdown <- tab[,3] - tab[,2]
  tab_up <- tab[order(-tab$diffup), ]
  tab_up <- tab_up[ 1:10, c(1,2,4) ]
  tab_down <- tab[order(-tab$diffdown), ]
  tab_down <- tab_down[ 1:10, c(1,2,5) ]
  tab_up_down <- rbind(tab_up[ ,1:2],tab_down[ ,1:2])
  label <- target[tab_up_down$index]
  
  ggplot(tab, aes(x=index,y=comparison)) + 
    geom_line(aes(y=comparison, group=2 , colour="comparison"), size=1 ) + 
    geom_line(aes(y=baseline, group=1 , colour="baseline"), size=1.4) + 
    geom_point(data=tab_up_down,aes(x=index,y=comparison),size=2) +  
    ggtitle(paste0("Comparison of ",comparison_name," vs ",baseline_name," (baseline)")) + 
    xlab(xlab) + ylab(ylab) + scale_colour_manual(name=legend, values=c(baseline="violetred3",comparison= "steelblue3")) + 
    theme(legend.position="bottom", text = element_text(size=text_size), axis.text=element_text(size= text_size) , axis.title= element_text(size= text_size), plot.title = element_text(size = title_size, hjust = 0.5 ), panel.background = element_rect(fill='white'), panel.grid.major = element_line(colour = "grey90")) + 
    geom_text_repel(data=tab_up_down ,aes(label=label), size=5 )  
}




