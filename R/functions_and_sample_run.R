#set of functions used to:
#1 - run fgsea on group of log2FC results with custom GMT file
#2 - reduce redundancy of fgsea results using gsea_elim function from
#@alserg in https://www.biostars.org/p/215956/


#####Necessary functions#####
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
eliminatePathways <- function(universe, pathways, pathway2name, isNonRandomPval, pvalThreshold=1e-3) {
  require(fgsea)
  require(data.table)
  messagef <- function(...) message(sprintf(...))
  pathways <- lapply(pathways, intersect, universe)
  mainPathways <- c()
  
  res <- list()
  
  pp <- names(pathway2name)[3] # for debug
  
  for (pp in names(pathway2name)) {
    messagef("Testing pathway '%s'", pathway2name[pp])
    
    if (length(mainPathways) == 0) {
      # nothing to compare with
      mainPathways <- pp
      next
    }    
    
    interesting <- TRUE
    pp1 <- mainPathways[1] # for debug
    # We want our new pathway to give us additional information compared to all others that we have already
    for (pp1 in mainPathways) {
      notExplainsPval <- min(isNonRandomPval(universe = pathways[[pp1]], 
                                             p = pathways[[pp]]),
                             isNonRandomPval(universe = setdiff(universe, pathways[[pp1]]), 
                                             p= pathways[[pp]]))
      res <- c(res, list(data.frame(
        checking=pp,
        reference=pp1,
        pval=notExplainsPval)))
      
      if (notExplainsPval > pvalThreshold) {
        
        messagef("Not adding pathway '%s', because we already have '%s'", 
                 pathway2name[pp],
                 pathway2name[pp1])
        
        eqPval <- min(
          isNonRandomPval(universe = setdiff(universe, pathways[[pp]]), 
                          pathways[[pp1]]),
          isNonRandomPval(universe = pathways[[pp]], 
                          pathways[[pp1]]))
        
        res <- c(res, list(data.frame(
          checking=pp1,
          reference=pp,
          pval=eqPval)))
        
        if (eqPval <= pvalThreshold) {
          messagef("And they are not equivalent: %s", eqPval)
        } else {
          messagef("They are equivalent: %s", eqPval)
        }
        interesting = FALSE
        break
      }
    }    
    
    if (interesting) {
      messagef("Adding pathway '%s', because we can", 
               pathway2name[pp])
      mainPathways <- c(mainPathways, pp)
      
      # OK it gives use new info, may be now we don't need some of already added pathways?
      
      # Now roles have changed, we test every interesting pathways agains out newly added one
      for (pp1 in setdiff(mainPathways, pp)) {
        notExplainsPval <- min(
          isNonRandomPval(universe = setdiff(universe, pathways[[pp]]), 
                          pathways[[pp1]]),
          isNonRandomPval(universe = pathways[[pp]], 
                          pathways[[pp1]]))
        res <- c(res, list(data.frame(
          checking=pp1,
          reference=pp,
          pval=notExplainsPval)))
        if (notExplainsPval > pvalThreshold) {
          messagef("Removing pathway %s, because we now have %s", 
                   pathway2name[pp1],
                   pathway2name[pp])
          mainPathways <- setdiff(mainPathways, pp1)
        }
      }
    }
  }
  
  list(mainPathways=mainPathways,
       table=rbindlist(res))
}

fgseaIsNonRandomPval <- function(pathways, ranks, nperm) {
  isNonRandomPval <- function(universe, p) {
    pathway <- intersect(p, universe)
    if (length(pathway) == 0 || length(pathway) == length(universe)) {
      return(1)
    }
    fgsea(pathways = list(pathway),
          stats=ranks[names(ranks) %in% universe], 
          nperm = nperm)$pval
  }
}

#function to use the two functions above on a previously obtained fgsea results list
redundance_reduce_fgsea <- function(fgseaRes, #list containing fgsea results
                                    ranks_df, #same used in "run_fgsea"
                                    gmt, #same used in "run_fgsea"
                                    keep_only_common=F, #same used in "run_fgsea"
                                    pval_cutoff=0.05, #enrichment cutoff
                                    gsea_elim #path to "gsea_elim.R" script
){
  source(gsea_elim)
  require(dplyr)
  
  gmt <- clusterProfiler::read.gmt(gmt)
  gmt_fgsea <- base::split(gmt[,2],gmt[,1])
  
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

#####Script to use the functions above#####
#Author: Thomaz Lüscher Dias
#Script to run fgsea on log2FC file using KEGG gmt
#Then, run redundancy reduction with gsea_elim.R script

#set working dir (change in your PC)
wd <- "~/Área de Trabalho/Papers_Helder/Reduce_redundancy_GSEA/"
setwd(wd)

#load functions from R/ folder
source("R/utils.R")
source("R/gsea_elim.R")

#run gsea with KEGG as pathways and edgeR combined log2FC file as input
res <- run_fgsea(ranks_df = "data/All_edgeR_combined_log2fc.txt",
                 gmt = "data/KEGG_2019_Human.txt",keep_only_common = F)

#run redundancy reduction for the results above
#(use same files for "ranks_df" and "gmt")
elim_res <- redundance_reduce_fgsea(fgseaRes = res,
                                    ranks_df = "data/All_edgeR_combined_log2fc.txt",
                                    gmt = "data/KEGG_2019_Human.txt",
                                    keep_only_common = F,
                                    pval_cutoff = 0.05,
                                    gsea_elim = "R/gsea_elim.R")

######STILL WORKING ON NEXT STEPS: filtration, ploting, etc################