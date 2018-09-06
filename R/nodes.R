select_edges <- function(A, xi, q,score,nnodes=1,discrete=FALSE){

  qpair <- NULL
  p <- nrow(A)
  pair <- expand.grid(i=seq_len(p),j=seq_len(p))
  pair <- pair[pair[,1]!=pair[,2],]
  idx <- sample(x=p*(p-1)/2, size=q, replace=FALSE)
  qpair <- pair[idx,]

  ir <- list()
  for(i in seq_len(q)) ir[[i]] <- 0:1
  x <- expand.grid(ir)
  nx <- nrow(x)
  sc <- c()
  for(k in seq_len(nx)){
    Ak <- A
    for(l in seq_len(q))
      Ak[qpair[l,1],qpair[l,2]] <- x[k,l]
    if(is.DAG(Ak)){
      if(discrete){
        dag <- graphAM(adjMat=Ak,edgemode='directed')
        sc <- c(sc, multinom.score(xi=xi, dag=dag))
      }
      else{
        dag <- as(Ak,'GaussParDAG')
        sc <- c(sc, score$global.score(dag))
      }
    }
    else
      sc <- c(sc,NA)
  }
  bad <- is.na(sc)
  sc <- sc[!bad]
  x <- x[!bad,]

  if(sum(!bad)>0){
    sc <- exp(sc - max(sc))
    sc <- sc/sum(sc)
    x <- cbind(x,sc)
  }
  else x <- NULL
  z <- list(index=qpair, score=x)

  return(z)
}

select_parents <- function(A, xi, q,ac,score,kappa=3,discrete=FALSE,
                           nnodes=1, verbose=1,progress.bar=progress.bar){

  p <- nrow(A)
  nodes <- sample(seq_len(p),size=q,replace=FALSE)

  nset <- ncol(ac)
  idx <- list()
  for(k in seq_len(q)) idx[[k]] <- seq_len(nset)
  grid <- expand.grid(idx)
  sc <- c()
  if(progress.bar) pb <- txtProgressBar(style=3)
  for(m in seq_len(nrow(grid))){
    Ak <- A
    for(k in seq_len(q)) Ak[,nodes[k]] <- ac[,grid[m,k]]
    if(sum(diag(Ak)) > 0) sc <- c(sc,NA)
    else if(!is.DAG(Ak))
      sc <- c(sc, NA)
    else{
      if(discrete){
        dag <- graphAM(adjMat=Ak,edgemode='directed')
        sc <- c(sc, multinom.score(xi=xi, dag=dag))
      }
      else{
        dag <- as(Ak,'GaussParDAG')
        sc <- c(sc, score$global.score(dag))
      }
    }
    if(progress.bar) setTxtProgressBar(pb, m/nrow(grid))
  }
  if(progress.bar) close(pb)
  bad <- is.na(sc)
  sc <- sc[!bad]
  grid <- grid[!bad,]

  if(sum(!bad)>0){
    sc <- exp(sc - max(sc))
    sc <- sc/sum(sc)
    grid <- cbind(grid,sc)
  }
  else{
    return(list(flag=TRUE))
  }

  id <- sample(seq_len(nrow(grid)),size=1,prob=grid[,q+1])
  Ak <- A
  for(k in seq_len(q)) Ak[,nodes[k]] <- ac[,grid[id,k]]

  return(list(Ak=Ak,flag=FALSE))
}

# returns matrix of columns each representing all possible parent sets

parent.sets <- function(nodes, kappa=3){

  p <- length(nodes)
  ac <- matrix(0,nrow=p,ncol=1)
  rownames(ac) <- nodes
  a <- diag(p)
  if(kappa > 0) ac <- cbind(ac,a)
  if(kappa > 1 & p > 1) for(i in seq(1,p-1)) for(j in seq(i+1,p))
    ac <- cbind(ac, a[,i] | a[,j])
  if(kappa > 2 & p > 2)
      for(i in seq(1,p-2)) for(j in seq(i+1,p-1)) for(k in seq(j+1,p))
        ac <- cbind(ac, a[,i] | a[,j] | a[,k])
  # columns of ac = repertoire of all parent sets
  if(kappa > 3) stop('Maximum in-degree is limited to <=3')
  return(ac)
}

# computes local score conditional to all possible parent sets for each node

local.score <- function(xi, ac, kappa, discrete=TRUE, scoring='ml',
                        score=NULL, g=NULL,
                        hyper=NULL, progress.bar=FALSE, ncores=1){

  nodes <- colnames(xi)
  p <- length(nodes)
  cache <- matrix(0, nrow=p, ncol=ncol(ac))
  rownames(cache) <- nodes

  cat('Computing local scores ...\n')

  bundle <- list(xi=xi, ac=ac, hyper=hyper, discrete=discrete,
                 scoring=scoring, g=g)  # parameter set

  nac <- ncol(ac)  # no. of parent sets
  if(ncores==1){
    if(progress.bar) pb <- txtProgressBar(style=3)
    lcache <- list()
    for(iac in seq_len(nac)){
      lcache[[iac]] <- fill.cache(iac, bundle)
      if(progress.bar) setTxtProgressBar(pb, iac/nac)
    }
    if(progress.bar) close(pb)
  }
  else{            # parallel
#   Rmpi::mpi.spawn.Rslaves(nslaves=ncores)
#   Rmpi::mpi.bcast.cmd(library(slimy))    # these are done outside
    Rmpi::mpi.bcast.Robj2slave(bundle)
    lcache <- Rmpi::mpi.applyLB(seq_len(nac), FUN=fill.cache, bundle)
#   Rmpi::mpi.close.Rslaves()
#   Rmpi::mpi.finalize()   # appears necessary for clean-up
  }

  for(iac in seq_len(nac))
    cache[,iac] <- lcache[[iac]]

  return(cache)
}

fill.cache <- function(iac, bundle){

  xi <- bundle$xi
  ac <- bundle$ac
  hyper <- bundle$hyper
  discrete <- bundle$discrete
  scoring <- bundle$scoring
  g <- bundle$g

  nodes <- colnames(xi)
  p <- length(nodes)

  pa <- nodes[which(ac[,iac]==1)]
  lcache <- c()
  for(i in seq_len(p)){
    w <- nodes[i]
    if(w %in% pa) sc <- NA
    else{
      wpa <- nodes[nodes %in% c(w,pa)]
      nw <- length(wpa)
      A <- matrix(0, nrow=nw, ncol=nw)
      rownames(A) <- colnames(A) <- wpa
      A[pa,w] <- 1
      if(!is.DAG(A)) sc <- NA
      else{
        if(discrete) sc <- multinom.local.score(xi, w, pa)
        else if(scoring=='ml')
          sc <- score$local.score(vertex=which(w==nodes),
                                       parents=match(pa,nodes))
        else if(scoring=='bge'){
          par <- hyper.par(xi=xi, nu=hyper$nu, alpha=hyper$alpha,
                             v=hyper$v, mu0=hyper$mu0, Sig=hyper$Sig)
          A1 <- A
          A1[pa,w] <- 1
          sc <- mvn.score(xi=xi, hyper.par=par, A=A1)
        }
        else if(scoring=='g-score'){
          sc <- g.score(xi=xi, node=w, pa=pa, g=g)
        }
        else stop('Unknown scoring')
      }
    }
    lcache[i] <- sc
  }

  return(lcache)
}
