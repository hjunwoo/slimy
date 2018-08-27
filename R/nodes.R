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
    if(isValidGraph(Ak,type='dag')){
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
    else if(!isValidGraph(Ak,type='dag'))
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

local.score <- function(xi, ac, kappa, discrete=TRUE, score=NULL,
                        progress.bar=FALSE){

  nodes <- colnames(xi)
  p <- length(nodes)
  cache <- matrix(0, nrow=p, ncol=ncol(ac))
  rownames(cache) <- nodes

  cat('Computing local scores ...\n')

  if(progress.bar) pb <- txtProgressBar(style=3)

  for(k in seq_len(ncol(ac))){
    pa <- nodes[which(ac[,k]==1)]
    for(i in seq_len(p)){
      w <- nodes[i]
      if(w %in% pa) sc <- NA
      else{
        wpa <- nodes[nodes %in% c(w,pa)]
        nw <- length(wpa)
        A <- matrix(0, nrow=nw, ncol=nw)
        rownames(A) <- colnames(A) <- wpa
        A[pa,w] <- 1
        if(!isValidGraph(A,type='dag')) sc <- NA
        else{
          if(discrete) sc <- multinom.local.score(xi, w, pa)
          else
            sc <- score$local.score(vertex=which(w==nodes),
                                       parents=match(pa,nodes))
        }
      }
      cache[i,k] <- sc
    }
    if(progress.bar) setTxtProgressBar(pb, k/ncol(ac))
  }
  if(progress.bar) close(pb)
  return(cache)
}
