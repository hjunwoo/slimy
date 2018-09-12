#' @export
read.bif <- function(file){

  x <- readLines(file)
  nodes <- c()
  edges <- list()
  node.levels <- list()
  p.tmp <- list()

  nx <- length(x)
  ix <- 1
  while(TRUE){
    xi <- x[[ix]]
    x1 <- strsplit(xi,split=' ')[[1]]
    x1 <- x1[x1!='']
    if(x1[1]=='variable')
      nodes <- c(nodes,x1[2])
    else if(x1[1]=='type'){
      nlevel <- as.numeric(x1[4])
      factors <- x1[seq(7,7+nlevel-1)]
      for(i in seq_len(nlevel)){
        nchar <- nchar(factors[i])
        if(substr(factors[i],start=nchar,stop=nchar)==',')
          factors[i] <- substr(factors[i],start=1,stop=nchar-1)
      }
      node.levels[[length(nodes)]] <- factors
      names(node.levels) <- nodes[seq_len(length(nodes))]
    }
    else if(x1[1]=='probability'){
      progeny <- x1[3]
      x2 <- strsplit(xi,split='[|]')[[1]][2]
      x3 <- strsplit(x2,split=')')[[1]][1]
      parents <- strsplit(x3,split='[,]')[[1]]
      parents <- vapply(parents,function(x){gsub(' ','',x,fixed=TRUE)},
                        character(1))
      if(sum(is.na(parents))==0){
        idx <- list()
        i <- 1
        for(p in parents){
          if(is.na(p)) next
          if(p %in% names(edges))
            edges[[p]]$edges <- c(edges[[p]]$edges,progeny)
          else{
            edges[[length(edges)+1]] <- list(edges=progeny)
            names(edges)[length(edges)] <- p
          }
          idx[[i]] <- node.levels[[p]]
          i <- i+1
        }
        grid <- expand.grid(idx)
        npa <- length(parents)
        names(grid) <- parents
        if(npa > 1){
          idx <- match(parents,nodes)
          grid <- grid[,order(idx)]
        }

        if(npa > 1)
          na <- apply(grid,1,function(x){y <- paste0(x,collapse=',');
                                     paste0('(',y,')')})
        else na <- grid[,1]
        for(k in seq_len(nrow(grid))){
          ix <- ix + 1
          xi <- x[[ix]]
          x1 <- strsplit(xi,split=')')[[1]]
          y <- strsplit(x1[-1],split=' ')[[1]]
          y <- y[y!='']
          x1 <- gsub(';','',y)
          x1 <- gsub(',','',x1)
          x1 <- as.numeric(x1)
          if(k==1) prob <- matrix(x1, nrow=1)
          else prob <- rbind(prob, matrix(x1, nrow=1))
        }
        rownames(prob) <- na
        gprob <- cbind(grid,prob)

        p.tmp[[length(p.tmp)+1]] <- gprob[,seq(npa+1,ncol(gprob))]
        names(p.tmp)[length(p.tmp)] <- progeny
      }
      else{
        ix <- ix + 1
        xi <- x[[ix]]
        y <- strsplit(xi,split=' ')[[1]]
        y <- y[y!=''][-1]
        x1 <- gsub(';','',y)
        x1 <- gsub(',','',x1)
        prob  <- matrix(as.numeric(x1),nrow=1)
        p.tmp[[length(p.tmp)+1]] <- prob
        names(p.tmp)[length(p.tmp)] <- progeny
      }
    }
    ix <- ix + 1
    if(ix > nx) break
  }
  names(node.levels) <- nodes
  for(v in nodes){
    if(!(v %in% names(edges))){
      edges[[length(edges)+1]] <- list(edges=NULL)
      names(edges)[length(edges)] <- v
    }
  }
  g <- graph::graphNEL(nodes=nodes, edgeL=edges, edgemode='directed')
  A <- as(g,'matrix')
  g <- graph::graphAM(adjMat=A, edgemode='directed')
  p.tmp <- p.tmp[match(nodes,names(p.tmp))]
  g@edgeData@data <- p.tmp
  g@nodeData@data <- node.levels

  return(g)
}

#' @export
distance <- function(m, mref, equivalence=TRUE){

  if(class(m)!='matrix') m <- as(m,'matrix')
  if(class(mref)!='matrix') mref <- as(mref,'matrix')
  if(nrow(m)!=nrow(mref) | ncol(m)!=ncol(mref))
    stop('distance error')

  if(equivalence){
#   m <- m + t(m)
#   mref <- mref + t(mref)
    m <- m | t(m)            #  could be bidirectional
    mref <- mref | t(mref)
    x <- m!=mref
    d <- sum(x[upper.tri(x)])
  }
  else d <- sum(m!=mref)

  return(d)
}

#' sensitivity and specificity of adjacency matrix with respect to ref
#' @export
senspec <- function(m, mref, equivalence=TRUE){

  if(class(m)!='matrix') m <- as(m,'matrix')
  if(class(mref)!='matrix') mref <- as(mref,'matrix')
  if(nrow(m)!=nrow(mref) | ncol(m)!=ncol(mref)) stop('distance error')

  if(equivalence){
#   m <- m + t(m)
#   mref <- mref + t(mref)
    m <- m | t(m)
    mref <- mref | t(mref)
    x <- m[upper.tri(m)]
    xref <- mref[upper.tri(mref)]
  }
  else{
    x <- m
    xref <- mref
  }
  tpr <- sum(x[xref==1])/sum(xref==1)
  tnr <- sum(x[xref==0]==0)/sum(xref==0)

  z <- list(sens=tpr,spec=tnr)
  return(z)
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
#' @param discrete Discrete or continuous
#' @param alpha Hyperparameter for dirichlet-distributed distributions for
#'        conditional probability
#' @param brange Lower and upper bound of regression coefficients for edges
#' @export
rgraph <- function(dag=NULL, nodes=NULL, mean.degree=1, max.degree=Inf,
                   levels=c('1','2'), discrete=TRUE, alpha=NULL,
                   brange=c(0.1,1)){

  if(!is.null(dag)){
    nodes <- nodes(dag)
    A1 <- as(dag,'matrix')
  }
  p <- length(nodes)
  nlevels <- length(levels)
  if(discrete) if(is.null(alpha)) alpha <- rep(1,nlevels)

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
      if(is.DAG(A1)) break
    }
  }
  prob <- vector('list',p)
  names(prob) <- nodes

  node_levels <- vector('list',p)
  names(node_levels) <- nodes

  for(k in seq_len(p)){
    np <- sum(A1[,k]!=0)
    if(discrete){
      if(np==0)
        prob[[k]] <- gtools::rdirichlet(n=1,alpha=alpha)
                              # marginal distribution
      else{
        prob[[k]] <- gtools::rdirichlet(n=nlevels^np, alpha=alpha)
      }
    }
    else if(np >0)
      prob[[k]] <- runif(n=np, min=brange[1], max=brange[2])
    if(np>0){
      if(discrete){
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
      else names(prob[[k]]) <- nodes[which(A1[,k]!=0)]
    }
    if(discrete)
      node_levels[[k]] <- levels
  }

  g <- graphAM(adjMat=A1,edgemode='directed')
  g@edgeData@data <- prob
  g@nodeData@data <- node_levels

  return(g)
}

#' Generate simulated data for discrete graph
#' @export
simulate.data <- function(dag, nsample, sd=1.0, progress.bar=FALSE){

  nodes <- nodes(dag)
  parents <- inEdges(dag)
  par <- dag@edgeData@data
  p <- length(nodes)
  levels <- dag@nodeData@data

  xi <- NULL
  if(progress.bar) pb <- txtProgressBar(style=3)
  for(k in seq_len(nsample)){
    x <- NULL
    for(i in seq_len(p)){
      x <- hfunc(x, i, parents, par, levels, sd=sd)
      if(length(x)==p) break
    }
    x <- x[match(nodes,names(x))]
    xi <- rbind(xi,x)
    if(progress.bar) setTxtProgressBar(pb, value=k/nsample)
  }
  if(progress.bar) close(pb)
  xi <- as.data.frame(xi)
  rownames(xi) <- seq_len(nsample)
  colnames(xi) <- nodes
  return(xi)
}

hfunc <- function(x, i, parents, par, levels, sd=1){

  nodes <- names(parents)
  if(nodes[i] %in% names(x)) return(x)
  np <- length(parents[[i]])
  values <- levels[[i]]
  discrete <- !is.null(values)
  pari <- par[[nodes[i]]]
  if(np==0){    # a root node
    if(discrete) prob <- pari[1,]
    else xm <- 0
  }
  else{
    for(k in seq_len(np)){
      if(!(parents[[i]][k]%in% names(x)))
        x <- hfunc(x, i=which(nodes==parents[[i]][k]),
                 parents=parents, par=par, levels=levels, sd=sd)
    }
    v <- x[parents[[i]]]
    if(discrete){
      irow <- paste0(v,collapse=',')
      if(np > 1) irow <- paste0('(',irow,')',collapse='')
      prob <- pari[irow,]
    }
    else{
      xm <- 0
      for(xp in parents[[i]])
        xm <- xm + pari[xp]*v[xp]
    }
  }
  if(discrete)
    ix <- sample(x=values, size=1, prob=prob)
  else
    ix <- rnorm(n=1, mean=xm, sd=sd)
  names(ix) <- nodes[i]

  x <- c(x,ix)
  return(x)

}

#' Generate Gaussian random graph
#' @param p Number of nodes
#' @param prob Edge probability
#' @export
rgraph.gauss <- function(p, prob=0.1){

  r <- r.gauss.pardag(p=p, prob=prob, top.sort=TRUE)
  nodes <- r$.nodes
  A <- as(r,'matrix')
  A <- apply(A,1:2, as.numeric)
  rownames(A) <- colnames(A) <- nodes

  dag <- graphAM(adjMat=A, edgemode='directed')
  dag@nodeData@data <- r$.params

  return(dag)
}

#' Generate simulated gaussian data
simulate.gauss <- function(dag, nsample){

  r <- new('GaussParDAG', nodes=nodes(dag),in.edges=inEdges(dag),
           params=dag@nodeData@data)
  xi <- r$simulate(nsample)

  return(xi)
}

#' Check if graph is DAG
#' @export
#'
is.DAG <- function(obj){

  if(class(obj)=='graphAM')
    obj <- as(obj,'matrix')
  else if(class(obj)=='data.frame')
    obj <- as.matrix(obj)

  if(class(obj)!='matrix') stop('Wrong object class in is.DAG')

  g <- igraph::graph_from_adjacency_matrix(obj)
  flag <- igraph::is_dag(g)

  return(flag)

}

#' Graph from posterior edge probability
#' @export

posterior.graph <- function(edge.prob, cut=0.5){

  nodes <- rownames(edge.prob)
  A <- apply(edge.prob,1:2,function(x){x>=0.5})
  g <- graphAM(adjMat=A,edgemode='directed')

  return(g)
}

#' ROC
#' @export
roc.graph <- function(edge.prob, mref){

  if(class(mref)!='matrix') mref <- as(mref,'matrix')
  a <- edge.prob + t(edge.prob)
  mref <- mref + t(mref)
  x <- a[upper.tri(a)]
  xref <- mref[upper.tri(mref)]
  roc <- pROC::roc(response=xref,predictor=x)

  return(roc)
}
