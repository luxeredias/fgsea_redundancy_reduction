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