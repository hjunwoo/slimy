#' Discrete data score in multinomial-dirichlet distribution
#' @param xi Data frame of discrete data with samples/nodes in rows/columns
#' @param dag Graph object of class \code{graphAM}
#' @param nprime Prior count
#' @export
multinom.score <- function(xi, dag, scoring='K2', nprime=1){

  if(scoring=='K2') nprime <- 1

  p <- ncol(xi)
  nodes <- nodes(dag)
  plist <- inEdges(dag)
  levels <- vector('list',p)
  names(levels) <- nodes
  for(i in seq_len(p))
    levels[[i]] <- levels(factor(xi[,i]))
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
    count <- apply(eg,1,function(x){sum(paste0(x,collapse=',')==cx)})
    eg <- cbind(eg,count)
    if(np > 0){
      by <- list()
      for(j in seq_len(np)) by[[j]] <- eg[,j]
      egs <- aggregate(eg$count, by=by, FUN=sum)
      names(egs) <- c(parents, 'count')
      egsc <- egs$count
    }
    else egsc <- sum(eg$count)

    nijk <- nprime
    if(scoring=='BDeu'){
      nijk <- nijk/length(levels[[nodes[i]]])
      if(np > 0) nijk <- nijk/np
    }
    sc <- sum(lgamma(nijk + eg$count)-lgamma(nijk))
    npij <- nijk*nrow(eg)
    sc <- sc + sum(lgamma(npij)-lgamma(npij+egsc))

    score <- score + sc
  }

  return(score)
}

multinom.local.score <- function(xi, node, parents, scoring='K2',
                                 hyper=NULL, nprime=1){

  if(scoring=='K2') nprime <- 1

  nodes <- c(parents,node)
  p <- length(nodes)
  levels <- vector('list',p)
  names(levels) <- nodes
  for(w in nodes)
    levels[[w]] <- levels(factor(xi[,w]))

  nsample <- nrow(xi)

  score <- 0
  np <- length(parents)
  tmp <- list()
  for(j in seq_len(np))
    tmp[[j]] <- levels[[parents[j]]]
  tmp[[np+1]] <- levels[[node]]
  eg <- expand.grid(tmp)     # enumerated states for (parent,i) set
  colnames(eg) <- nodes
  x <- xi[,nodes]
  if(is.null(dim(x)))
    cx <- as.character(x)
  else
    cx <- apply(x,1,function(x){paste0(x,collapse=',')})
  count <- apply(eg,1,function(x){sum(paste0(x,collapse=',')==cx)})
  eg <- cbind(eg,count)
  if(np > 0){
    by <- list()
    for(j in seq_len(np)) by[[j]] <- eg[,j]
    egs <- aggregate(eg$count, by=by, FUN=sum)
    names(egs) <- c(parents, 'count')
    egsc <- egs$count
  }
  else egsc <- sum(eg$count)

  nijk <- nprime
  if(scoring=='BDeu'){
    nijk <- nijk/length(levels[[node]])
    if(np > 0) nijk <- nijk/np
  }
  sc <- sum(lgamma(nijk + eg$count)-lgamma(nijk))
  npij <- nijk*nrow(eg)
  sc <- sc + sum(lgamma(npij)-lgamma(npij+egsc))

  return(sc)
}

multinom.cache.score <- function(dag, ac, cache){

  e <- 0
  nodes <- nodes(dag)
  parents <- inEdges(dag)
  for(w in nodes){
    pa <- parents[[w]]
    npa <- length(pa)
    ipa <- nodes %in% pa
    if(npa>0)
      iac <- which(apply(ac,2,function(x){all(x==ipa)}))
    else
      iac <- 1
    e <- e + cache[w,iac]
  }
  return(e)
}
