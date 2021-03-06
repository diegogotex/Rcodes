cool_table <- function(enrichment.obj,enrich.type,up,down){
  #this function only works for gene SYMBOL

  #selecting the enrichment type
  # 1 for ClusterProfiler
  # 2 for enrichR
  if(enrich.type == 1){
    sep = "/"
    splited_col <- strsplit(enrichment.obj@result$geneID, sep)
  } else if(enrich.type == 2){
    sep = ";"
    splited_col <- strsplit(enrichment.obj$Genes, sep) 
  } 
  
  #creating table to store the values 
  tab <- data.frame(matrix(nrow = length(splited_col), ncol = 4)) #creating a data.frame
  colnames(tab) <- c("up","down", "up_genes", "down_genes")       #giving names to the 4 columns
   
  #loop to identify which and how many gene SYMBOL are into the vectors
  for (i in 1:length(splited_col)){
    up2 <- subset(up, up %in% splited_col[[i]])
    down2 <- subset(down, down %in% splited_col[[i]])
    tab[i,] <- c(length(up2), length(down2), paste(up2, collapse = sep), paste(down2, collapse = sep))
  }
  
  #creating the new tab by collapsing the data frame with the enrichment and the df with the genes information
  if(enrich.type == 1){
    new_tab <- cbind(enrichment.obj@result, tab)
    #a partir daqui é novo
    new_tab[,4] <- sub("*\\/.*","",new_tab$BgRatio)
    colnames(new_tab)[4] <- "Term_size"
    new_tab <- new_tab[,c(1,2,4,6,10,11,12,13)]
    
  } else if(enrich.type == 2){
    new_tab <- cbind(enrichment.obj, tab)
    #a partir daqui é novo
    new_tab[,2] <- sub(".*\\/","",new_tab$Overlap)
    colnames(new_tab)[2] <- "Term_size"
    new_tab <- new_tab[,c(1,2,4,10,11,12,13)]
  } 
  
  return(new_tab)
  
  
  
}
