# Partition parent sets using GM algorithm
#
partition.pset <- function(A, xi, q, ac, path, kappa=3, progress.bar=FALSE,
                           discrete=TRUE, score){

  p <- nrow(A)
  nodes <- rownames(A)
  W <- sample(nodes, size=q)
  Ag <- A
  Cgbar <- path
  for(w in W)
    for(k in nodes[which(Ag[,w]>0)])   # remove edge k -> w
      Cgbar <- Cgbar - outer(Cgbar[,k],Cgbar[w,],'*')
  Ag[,W] <- 0     # matrix of subgraph Gbar with incoming edges into W

  dew <- ndw <- vector('list',q)
  names(dew) <- names(ndw) <- W

  nd <- NULL
  for(w in W){
    dew[[w]] <- nodes[which(Cgbar[w,]>0)]
    ndw[[w]] <- nodes[which(Cgbar[w,]==0)]
    nd <- union(nd,ndw[[w]])
  }

  H <- matrix(0,nrow=q,ncol=q)           # generate all dags for W
  rownames(H) <- colnames(H) <- W
  hac <- parent.sets(nodes=W, kappa=kappa)
  nset <- ncol(hac)
  idx <- list()
  for(k in seq_len(q)) idx[[k]] <- seq_len(nset)
  grid <- expand.grid(idx)
  colnames(grid) <- W

  Pawgh <- list()   # list of parent sets of w partitioned into H
  eta <- 0          # partition count

  KH <- NULL
  ngrid <- nrow(grid)
  if(progress.bar) pb <- txtProgressBar(style=3)
  for(m in seq_len(ngrid)){

    Hk <- H
    for(w in W) Hk[,w] <- hac[,grid[m,w]]
    if(!isValidGraph(Hk,type='dag')) next
    eta <- eta + 1
    Pawgh[[eta]] <- vector('list', q)
    names(Pawgh[[eta]]) <- W

    lkh <- 0
    for(w in W){
      nx <- sum(Hk[,w])
      x <- rownames(Hk)[Hk[,w]==1]   # pa_w^H
      deghw <- NULL
      for(ix in x) deghw <- union(deghw, dew[[ix]])
      y <- rownames(Hk)[Hk[,w]==0]   # y not pa_w^H
      denhw <- NULL
      for(iy in y) denhw <- union(denhw, dew[[iy]])
      qwh <- union(nd, deghw)
      qwh <- qwh[!(qwh%in%denhw)]    # (nd U de_w) \ de_-w

      pawA <- NULL
      for(v in nodes[!(nodes %in% qwh)])
        pawA <- union(pawA, which(ac[w,]==0 & ac[v,]==1)) # Pa_w^{all}(v)

      pawB <- seq_len(ncol(ac))
      for(ix in x){
        rwh <- dew[[ix]]
        rwh <- rwh[!(rwh%in%denhw)]  # R_w,x = de_x \ de_{-w}
        U <- NULL
        for(r in rwh)
          U <- union(U, which(ac[w,]==0 & ac[r,]==1))
        pawB <- intersect(pawB, U)
      }

      if(nx==0){
        pawgh <- which(ac[w,]==0)
        pawgh <- pawgh[!(pawgh %in% pawA)]
      }
      else
        pawgh <- pawB[!(pawB %in% pawA)]  # these are column indices of ac

      Pawgh[[eta]][[w]] <- pawgh
      sc <- NULL
      aw <- A
      plist <- w
      for(i in pawgh)
        plist <- union(plist, rownames(ac)[which(ac[,i]==1)])
      aw <- aw[plist,plist]
      aw[,w] <- 0
      for(i in pawgh){
        pa <- rownames(ac)[which(ac[,i]==1)]
        aw[pa,w] <- 1
        gw <- graphAM(adjMat=aw, edgemode='directed')
        e <- multinom.score(xi=xi[,plist], dag=gw)
        sc <- c(sc, e)
      }
      lsc <- log(sum(exp(sc-max(sc)))) + max(sc)
      lkh <- lkh + lsc
    }
    KH <- c(KH,lkh)
    if(progress.bar) setTxtProgressBar(pb, m/ngrid)
  }
  KH <- KH - max(KH)
  if(progress.bar) close(pb)

  h <- list(PaH=Pawgh, KH=KH)
  return(h)
}

# path count matrix for transitive closure

path.count <- function(dag){

  nodes <- nodes(dag)
  p <- length(nodes)
  desc <- edges(dag)

  C <- diag(p)  # path count matrix
  rownames(C) <- colnames(C) <- nodes
  Clist <- vector('list', p*(p-1))  # path list (Clist_ij = list(i->j))
  ij <- NULL
  for(i in nodes) for(j in nodes){
    if(i==j) next
    ij <- c(ij,paste(i,j,sep='->'))
  }
  names(Clist) <- ij

  for(i in nodes){
    path <- NULL
    Clist <- cfunc(Clist, i, desc, path)
  }

  C <- diag(p)
  rownames(C) <- colnames(C) <- nodes

  for(k in seq_len(p*(p-1))){
    x <- Clist[[k]]
    if(is.null(x)) next
    if(sum(duplicated(x))>0) stop('Error in path.count')
    z <- strsplit(ij[k],split='->')[[1]]
    start <- z[1]
    end <- z[2]
    C[start,end] <- C[start,end] + length(x)
  }

  return(list(C=C, Clist=Clist))
}

cfunc <- function(Clist, i, desc, path){

  path <- c(path, i)
  de <- desc[[i]]
  np <- length(de)
  if(np==0) return(Clist)

  for(k in de){
    x <- paste(path[1],k,sep='->')
    Clist[[x]] <- c(Clist[[x]], paste(c(path,k),collapse=','))
    Clist <- cfunc(Clist, k, desc, path)
  }
  return(Clist)
}

# updates the path count matrix C after change from A0 to A1
path.update <- function(A0, A1, C){

  diff <- which(A0!=A1,arr.ind=TRUE)
  if(is.null(diff)) nd <- 1
  else nd <- nrow(diff)
  C1 <- C
  sn <- NULL
  for(m in seq_len(nd)){
    i <- diff[m,1]
    j <- diff[m,2]
    if(A0[i,j]==0 & A1[i,j]==1) sign <- 1       # edge added
    else if(A0[i,j]==1 & A1[i,j]==0) sign <- -1 # edge removed
    sn <- c(sn,sign)
  }
  diff <- cbind(diff,sn)
  if(nd>1) diff <- diff[order(diff[,3]),]

  for(m in seq_len(nd)){
    i <- diff[m,1]
    j <- diff[m,2]
    sign <- diff[m,3]
    C1 <- C1 + sign*outer(C1[,i],C1[j,],'*')
  }
  if(sum(C1<0)!=0) stop('Error in path update')
  return(C1)
}