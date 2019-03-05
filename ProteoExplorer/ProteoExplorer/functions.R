`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

convert_ID <- function(ID=1,data) {
  swissprot_ID_HGNC <- read.csv("./data/swissprot_ID_HGNC", row.names=1)
  if(ID==2) { 
    common <- intersect( swissprot_ID_HGNC[ ,1] , rownames(data) )
    matrix <- data[common, ]
    x <- swissprot_ID_HGNC[which(swissprot_ID_HGNC[ ,1] %in% common), 2 ]
    rownames(matrix) <- x
    output <- matrix
  } else { output <- data }
  return(output)
}

convert_ID_col <- function(ID=1,data) {
  swissprot_ID_HGNC <- read.csv("./data/swissprot_ID_HGNC", row.names=1)
  if(ID==2) { 
    common <- intersect( swissprot_ID_HGNC[ ,1] , rownames(data) )
    matrix <- data.matrix(data[common, ]) ; colnames(matrix)[1] <- "x"
    x <- swissprot_ID_HGNC[which(swissprot_ID_HGNC[ ,1] %in% common), 2 ]
    rownames(matrix) <- x
    output <- matrix
  } else { output <- data }
  return(output)
}


scatter_plot <- function(gene,real_gene,df,label_X,label_Y,text_size=26,title_size=2.5 ) {
  if(real_gene %not in% c(gene,"")){
    ggplot(df ) + ggtitle(paste0(real_gene," is not in common for those 2 layers"))
    
  } else {
    
    # equation, correlation and p value
    out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
    lm_eqn <- function(df){
      m <- lm(b ~ a, df);
      eq <- substitute(~~italic("r")~"="~r*","~~italic("p")~"="~p,
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r = format(unname(r), digits = 2),
                            p = format(p, digits=2)))
      as.character(as.expression(eq));                 
    }
    
    g <- ggplot(df, aes(a, b, color = "#000033" )) + 
      geom_point(shape = 16, size = 3, show.legend = FALSE, alpha = .8 ) + geom_smooth(method=lm,se=F,show.legend=F) + 
      labs(x = label_X, y=label_Y ) + ggtitle(gene) + 
      theme(legend.position="none",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
            panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
      geom_text(x = min(df$a) + 0.5*(max(df$a)-min(df$a)), y = min(df$b) , label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 7) 
    g
  }
}

density_plot <- function(df,gene,text_size=20,title_size=2) {
  if(gene %not in% c(rownames(df),"")){
    df <- cbind(df, rep(colnames(df)[1], length(rownames(df)))) ; df <- data.frame(df) ; df[ ,1] <- as.numeric(as.character(df[ ,1]))
    df <- df[complete.cases(df) , ];
    ggplot(df ) + ggtitle(paste0("The data is lacking for ",gene))
    
  } else {
    if(gene==""){gene <- rownames(df)[1]}
    gene_performance <- df[gene, ] 
    df <- cbind(df, rep(colnames(df)[1], length(rownames(df)))) ; df <- data.frame(df) ; df[ ,1] <- as.numeric(as.character(df[ ,1]))
    df <- df[complete.cases(df) , ]; df_ordered <- df[ order(-df[ ,1]) , ] ; rank <- which(df_ordered[ ,1] == gene_performance)[1]
    ggplot(df, aes(x=x)) + geom_density(fill="grey95") + ggtitle(paste0(gene)) + xlab("Prediction performance (r)") + 
      theme(legend.position="none",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), 
            plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
            panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour ='white') )  + 
      geom_vline(xintercept=gene_performance, color="red") + geom_text(x = min(df[ ,1]) + 0.5*(max(df[ ,1])-min(df[ ,1])), aes(x=gene_performance, 
    label=paste0("r=",signif(df[gene,1],2),",  Rank ",rank,"/",length(df[ ,1])), y=0.05), size=7 ) 
  }
  
}
