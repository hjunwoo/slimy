#' Discrete data score in multinomial-dirichlet distribution
#' @param xi Data frame of discrete data with samples/nodes in rows/columns
#' @param dag Graph object of class \code{graphAM}
#' @param nprime Prior count
#' @export
multinom.score <- function(xi, dag, nprime=1){

  p <- ncol(xi)
  nodes <- nodes(dag)
  plist <- inEdges(dag)
  levels <- dag@nodeData@data
  if(length(levels)==0){
    levels <- vector('list',p)
    names(levels) <- nodes
    for(i in seq_len(p))
      levels[[i]] <- levels(factor(xi[,i]))
  }
  nsample <- nrow(xi)

  score <- 0
  for(i in seq_len(p)){
    parents <- plist[[i]]
    np <- length(parents)
    tmp <- list()
    for(j in seq_len(np))
      tmp[[j]] <- levels[[parents[j]]]
    tmp[[np+1]] <- levels[[i]]
    eg <- expand.grid(tmp)     # enumerated states for (parent,i) set
    colnames(eg) <- c(parents,nodes[i])
    x <- xi[,c(parents,nodes[i])]
    if(is.null(dim(x)))
      cx <- as.character(x)
    else
      cx <- apply(x,1,function(x){paste0(x,collapse=',')})
    count <- NULL
    for(m in seq_len(nrow(eg))){
      ijk <- paste0(eg[m,], collapse=',')
      count <- c(count,sum(cx==ijk))
    }
    eg <- cbind(eg,count)
    if(np > 0){
      by <- list()
      for(j in seq_len(np)) by[[j]] <- eg[,j]
      egs <- aggregate(eg$count, by=by, FUN=sum)
      names(egs) <- c(parents, 'count')
      egsc <- egs$count
    }
    else egsc <- sum(eg$count)
    sc <- sum(lgamma(nprime+eg$count)-lgamma(nprime))
    npij <- nprime*nrow(eg)
    sc <- sc + sum(lgamma(npij)-lgamma(npij+egsc))

    score <- score + sc
  }

  return(score)
}
