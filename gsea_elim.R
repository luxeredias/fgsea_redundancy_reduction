# eliminating indirect hits based on bayesian networks ideas
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
