###############################################################################
####                    6. PLOT SUBGRAPHS FOR SELECTED CCs                 ####
###############################################################################
### Define function to plot k-partite subgraphs representing connected components 
### of interest, using igraph package

plotCCs <- function(cc_ids){
  
  for (cc_id in 1:length(cc_ids)){
    
    g <- igraph::graph_from_incidence_matrix(cc.subseq.subincM$cc.subincM[[cc_id]])
    V(g)$label.cex <- 0.5
    V(g)$label.color <- "black"
    V(g)$color <- rep("#0072B2", length(names(as.list(V(g))))) # Blue for RNAseq subseq nodes
    V(g)$color[grep("ENST",names(as.list(V(g))))] <- "#D55E00" # Orange for transcript nodes
    V(g)$color[grep("P_",names(as.list(V(g))))] <- "#009E73" # Green for peptide nodes
    opar <- par()$mar; par(mar=rep(0, 4)) #Give the graph lots of room
    #igraph::plot.igraph(g, edge.color="gray30", cex=.5,layout=layout_as_bipartite)
    igraph::plot.igraph(g, edge.color="gray30", cex=.5, vertex.size=5)
    #title(paste0("CC_", CC_id), line = -1)
    title(paste0("CC_", CCs_toPlot[CC_id]), line = -1)
    par(mar=opar)
  }
  dev.off()
}