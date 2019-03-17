

group_prediction_in_family <- function(pred) {

  gene_families_HGNC_092016 <- read.delim2("~/Documents/RWTH_Aachen/GENERAL_DATA/gene_families_HGNC_092016.txt")  ## Hugo gene name, families 
  
  result <- gene_families_HGNC_092016[gene_families_HGNC_092016$Approved.Symbol %in% rownames(pred), ]
  result <- result[ ,c("Approved.Symbol","Gene.family.description") ]
  pred <- pred[as.character(result$Approved.Symbol), ]
  result <- cbind(result, pred )
  
  to_remove <- names(which( table( result$Gene.family.description ) < 10 ))
  result <- result[ result$Gene.family.description %not in% to_remove , ]
  result <- result[ ,-1]
  result$Gene.family.description <- as.character(result$Gene.family.description)
  
  average <- aggregate(result, by=list(result$Gene.family.description), FUN="mean")
  average <- average[ ,-2] ; colnames(average) <- c("family","predictability")
  
  sd <- aggregate(result, by=list(result$Gene.family.description), FUN="sd")
  sd <- sd[ ,-2] ; colnames(sd) <- c("family","predictability")
  return(list(average,sd))
}

scatter_plot_pred_by_family <- function(df,inverse_sum_CV,title,x_lab,y_lab ,text_size=33,title_size=3.5) {
  # equation, correlation and p value
  df <- as.data.frame(df) ; colnames(df) <- c("a","b")
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(~~italic("r")~"="~r*","~~italic("p")~"="~p,
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(r, digits = 2),
                          p = format(p, digits=2)))
    as.character(as.expression(eq));                 
  }
  top_hits1  <- df[ which( (df$a>0.65) & (df$b>0.65) ) , ] 
  top_hits2  <- df[ which( (df$a<0.3) & (df$b<0.3) ) , ] 
  top_hits <- rbind(top_hits1,top_hits2)
  label_top_hits <- rownames(top_hits)    ## (df$a + df$b)<0.6 )
  rest  <- df[ which(rownames(df) %not in% label_top_hits) , ] ; label_rest <- rownames(rest)
  g <- ggplot(df, aes(a, b, color=b)) + 
    geom_point(shape = 16, size = 5*inverse_sum_CV, show.legend = FALSE, alpha = .4 ) +  geom_smooth(method=lm,se=F,show.legend=F) + labs(x = x_lab, y=y_lab ) + ggtitle(title) + 
    theme(legend.position="none",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = min(df$a) + 0.57*(max(df$a)-min(df$a)), y = min(df$b) + 0*(max(df$b)-min(df$b)), label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 15) +
    geom_text_repel(data=top_hits ,aes(label=label_top_hits), size=9, color="black") + 
    scale_color_gradient(low = "#990000", high = "#990000" )
  g
}
