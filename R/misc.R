# accumulate list of graphs for MAP estimuation

store.map <- function(Map, A, lkl){

  nl <- length(Map)
  for(k in seq_len(nl)){
    x <- Map[[k]]
    if(sum(x$A-A)==0) if(x$lkl < lkl) x$lkl <- lkl
           # new lkl for existing graph found
  }
}
