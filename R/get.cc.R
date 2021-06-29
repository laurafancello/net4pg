###############################################################################
####           4. BUILD GRAPH AND CALCULATE CONNECTED COMPONENTS           ####
###############################################################################
# Define function to build a graph object representing origin-to-origin 
# connections from adjacency matrix. Calculate graph connected components.

get.cc <- function(adjM){
  
  ### Generate graph and calculate components
  print("Calculating connected components")
  g <- graphAM(adjM, edgemode='undirected', values=NA)
  cc.origins <- connComp(as(g, 'graphNEL'))
  
  result <- list(g=g, cc.origins=cc.origins)
  return(result)
}