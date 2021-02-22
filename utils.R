#set of functions used to:
#1 - run fgsea on group of log2FC results with custom GMT file
#2 - reduce redundancy of fgsea results using gsea_elim function from
#@alserg in https://www.biostars.org/p/215956/
#3 - plot fgsea results
ranks_df="data/All_edgeR_combined_log2fc.txt"
gmt="data/KEGG_2019_Human.txt"
keep_only_common=F
#1 - run fgsea on group of log2FC results with custom GMT file
run_fgsea <- function(ranks_df, #path to ranks df. Ex.: "data/ranks.txt"
                      gmt, #path to gmt. Ex.: "data/KEGG.gmt"
                      keep_only_common=F #keep only genes w/ logFC values in all samples
                      ){
  
  #message("loading packages: clusterProfiler,fgsea,dplyr,tidyverse,data.table")
  require(clusterProfiler,quietly = T)
  require(fgsea,quietly = T)
  require(dplyr,quietly = T)
  require(tidyverse,quietly = T)
  require(data.table,quietly = T)
  #message("packages loaded")
  
  #load ranks file into dataframe
  ranks_df <- data.table::fread(ranks_df)
  
  #load gmt file using read.gmt from clusterProfiler
  gmt <- clusterProfiler::read.gmt(gmt)
  
  #split df gmt into list gmt suited for fgsea function
  gmt_fgsea <- base::split(gmt[,2],gmt[,1])
  
  #if arg T, will keep only genes that have log2FC values across all samples
  if(keep_only_common){
    ranks_mx <- ranks_mx[complete.cases(ranks_mx)]
  }
  
  #transform ranks_df to matrix with GeneSymnol as rownames
  #obs: first column of ranks_df must contain gene symbols/names compatible
  #with gmt file
  colnames(ranks_df)[1] <- "GeneSymbol"
  ranks_mx <- ranks_df %>%
    column_to_rownames("GeneSymbol") %>%
    as.matrix()
  
  #create list of ranks in each sample (columns in ranks_mx) removing NA values if present
  ranks <- apply(ranks_mx,MARGIN = 2,FUN = function(x) x[!is.na(x)])
  
  #run fgsea using lapply
  fgsea_res <- lapply(ranks,fgsea,pathways=gmt_fgsea,nperm=1000)
  
  #return list containing dataframes with each fgsea result
  #each list element has the name of the sample it corresponds to
  return(fgsea_res)
}


#2 - reduce redundancy of fgsea results using gsea_elim function from
#@alserg in https://www.biostars.org/p/215956/
redundance_reduce_fgsea <- function(fgseaRes, #list containing fgsea results
                                    ranks_df, #same used in "run_fgsea"
                                    gmt, #same used in "run_fgsea"
                                    keep_only_common=F, #same used in "run_fgsea"
                                    pval_cutoff=0.05, #enrichment cutoff
                                    gsea_elim #path to "gsea_elim.R" script
                                    ){
  source(gsea_elim)
  
  gmt <- clusterProfiler::read.gmt(gmt)
  gmt_fgsea <- split(gmt[,2],gmt[,1])
  
  if(keep_only_common){
    ranks_mx <- ranks_mx[complete.cases(ranks_mx)]
  }
  
  ranks_df <- fread(ranks_df)
  colnames(ranks_df)[1] <- "GeneSymbol"
  ranks_mx <- ranks_df %>%
    column_to_rownames("GeneSymbol") %>%
    as.matrix()
  ranks <- apply(ranks_mx,MARGIN = 2,FUN = function(x) x[!is.na(x)])
  
  pathway2name <- lapply(fgseaRes,FUN = function(x){x[pval < pval_cutoff][order(pval), setNames(pathway, pathway)]})
  
  #run redundancy eliminatio on each fgsea result
  samples <- names(fgseaRes)
  
  elimRes <- sapply(samples,simplify = F,FUN = function(sample){
    fgsea_sample <- fgseaRes[sample][[1]]
    ranks_sample <- ranks[sample][[1]]
    pathway2name_sample <- pathway2name[sample][[1]]
    elim <- eliminatePathways(universe = names(ranks_sample),
                              pathways = gmt_fgsea,
                              pathway2name = pathway2name_sample,
                              isNonRandomPval = fgseaIsNonRandomPval(gmt_fgsea,ranks_sample,nperm = 5000))
    return(elim)
  })
  return(elimRes)
}

