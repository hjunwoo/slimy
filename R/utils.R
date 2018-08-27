#' @export
read.bif <- function(file){

  x <- readLines(file)
  nodes <- c()
  edges <- list()

  for(xi in x){
    x1 <- strsplit(xi,split=' ')[[1]]
    if(x1[1]=='variable')
      nodes <- c(nodes,x1[2])
    else if(x1[1]=='probability'){
      progeny <- x1[3]
      x2 <- strsplit(xi,split='[|]')[[1]][2]
      x3 <- strsplit(x2,split=')')[[1]][1]
      parents <- strsplit(x3,split='[,]')[[1]]
      parents <- vapply(parents,function(x){gsub(' ','',x,fixed=TRUE)},
                        character(1))
      for(p in parents){
        if(is.na(p)) next
        if(p %in% names(edges))
          edges[[p]]$edges <- c(edges[[p]]$edges,progeny)
        else{
          edges[[length(edges)+1]] <- list(edges=progeny)
          names(edges)[length(edges)] <- p
        }
      }
    }
  }
  for(v in nodes){
    if(!(v %in% names(edges))){
      edges[[length(edges)+1]] <- list(edges=NULL)
      names(edges)[length(edges)] <- v
    }
  }
  g <- graphNEL(nodes=nodes, edgeL=edges, edgemode='directed')
  return(g)
}

#' @export
distance <- function(m, mref, equivalence=TRUE){

  if(nrow(m)!=nrow(mref) | ncol(m)!=ncol(mref)) stop('distance error')

  if(equivalence){
    m <- m + t(m)
    mref <- mref + t(mref)
    x <- m!=mref
    d <- sum(x[upper.tri(x)])
  }
  else d <- sum(m!=mref)

  return(d)
}

mat2nel <- function(nodes, A){

  p <- nrow(A)
  edL <- list()
  for(i in seq_len(p))
    edL[[i]] <- list(edges=which(A[i,]!=0))
  names(edL) <- nodes
  gR <- graphNEL(nodes=nodes,edgeL=edL,edgemode='directed')

  return(gR)
}

#' Generate random graph with discrete distributions
#' @param nodes Vector of node names
#' @param mean.degree Mean in-degree per node (Poisson lambda)
#' @param max.degree Upper bound for in-degree
#' @param levels Levels for each variable
#' @param alpha Hyperparameter for dirichlet-distributed distributions for
#'        conditional probability
#' @export
rgraph <- function(dag=NULL, nodes=NULL, mean.degree=1, max.degree=Inf,
                   levels=c('1','2'), alpha=NULL){

  if(!is.null(dag)){
    nodes <- nodes(dag)
    A1 <- as(dag,'matrix')
  }
  p <- length(nodes)
  nlevels <- length(levels)
  if(is.null(alpha)) alpha <- rep(1,nlevels)

  if(is.null(dag)){
    A <- matrix(0, nrow=p, ncol=p)
    rownames(A) <- colnames(A) <- nodes

    while(1){
      A1 <- A
      for(k in seq_len(p)){
        while(TRUE){
          L <- rpois(n=1, lambda=mean.degree)
          if(L<=min(max.degree,p-1)) break
        }
        if(L==0) next
        m <- sample(x=nodes[-k],size=L,replace=FALSE)
        A1[m,k] <- 1
      }
      if(isValidGraph(A1,type='dag')) break
    }
  }
  prob <- vector('list',p)
  names(prob) <- nodes

  node_levels <- vector('list',p)
  names(node_levels) <- nodes

  for(k in seq_len(p)){
    np <- sum(A1[,k]!=0)
    if(np==0)
      prob[[k]] <- rdirichlet(n=1,alpha=alpha)  # marginal distribution
    else{
      prob[[k]] <- rdirichlet(n=nlevels^np, alpha=alpha)
    }
    if(np>0){
      tmp <- list()
      for(l in seq_len(np)) tmp[[l]] <- levels
      eg <- expand.grid(tmp)
      if(np==1) ename <- levels
      else
        ename <- apply(eg,1,
                     function(x){z <- paste0(x,collapse=',');
                                 return(paste0('(',z,')',collapse=''))})
      rownames(prob[[k]]) <- ename
    }
    node_levels[[k]] <- levels
  }

  g <- graphAM(adjMat=A1,edgemode='directed')
  g@edgeData@data <- prob
  g@nodeData@data <- node_levels

  return(g)
}

#' Generate simulated data for discrete graph
#' @export
simulate.data <- function(dag, nsample){

  nodes <- nodes(dag)
  parents <- inEdges(dag)
  par <- dag@edgeData@data
  p <- length(nodes)
  levels <- dag@nodeData@data

  xi <- NULL
  for(k in seq_len(nsample)){
    x <- NULL
    for(i in seq_len(p)){
      x <- hfunc(x, i, parents, par, levels)
      if(length(x)==p) break
    }
    x <- x[match(nodes,names(x))]
    xi <- rbind(xi,x)
  }
  xi <- as.data.frame(xi)
  rownames(xi) <- seq_len(nsample)
  colnames(xi) <- nodes
  return(xi)
}

hfunc <- function(x, i, parents, par, levels){

  nodes <- names(parents)
  if(nodes[i] %in% names(x)) return(x)
  np <- length(parents[[i]])
  values <- levels[[i]]
  if(np==0)
    prob <- par[[i]][1,]
  else{
    for(k in seq_len(np)){
      if(!(parents[[i]][k]%in% names(x)))
        x <- hfunc(x, i=which(nodes==parents[[i]][k]),
                 parents, par, levels)
    }
    v <- x[parents[[i]]]
    irow <- paste0(v,collapse=',')
    if(np > 1) irow <- paste0('(',irow,')',collapse='')
    prob <- par[[i]][irow,]
  }
  ix <- sample(x=values, size=1, prob=prob)
  names(ix) <- names(par)[i]
  x <- c(x,ix)
  return(x)

}
