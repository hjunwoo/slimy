#' @export
mc.sample <- function(xi, dag=NULL, ref=NULL,
                      discrete=FALSE,
                      nstep=1000, Big=100, verbose=1, method='gibbs',
                      scoring='ml', representation='parent.set',
                      progress.bar=FALSE,
                      hyper=NULL, q=1,npr=100, kappa=3){

  m <- nrow(xi) # no. of samples
  p <- ncol(xi) # no. of nodes

  if(!is.null(ref)){
    Aref <- as(ref,'matrix')
    Aref <- apply(Aref,1:2,as.numeric)
    colnames(Aref) <- colnames(xi)
  }

  if(is.null(dag)){
    if(discrete)
      dag <- rgraph(nodes=nodes,mean.degree=1.0)
    else
      dag <- r.gauss.pardag(p=p, prob=.1)
  }
  A <- as(dag,'matrix')
  A <- apply(A,1:2,as.numeric)
  colnames(A) <- colnames(xi)

  if(progress.bar) pb <- txtProgressBar(style=3)
  if(verbose==3){
    sumd <- 0
    cnt <- 0
  }

  if(discrete)
    s1 <- sc <- multinom.score(xi=xi, dag=dag)
  else{
    if(scoring=='ml') score <- new('GaussL0penObsScore',xi)
    else s1 <- sc <- mvn.score(xi=xi, hyper=hyper, A=A)
  }

  if(method=='gibbs'){
    ac <- parent.sets(nodes=colnames(xi), kappa)
    if(representation=='gm')
      path <- path.count(dag=dag)$C
  }

  istep <- 0
  while(TRUE){

    if(method=='metropolis'){
      nbr <- neighbor(A)
      k <- floor(runif(n=1)*nrow(nbr))+1
      A1 <- A
      i <- nbr[k,1]
      j <- nbr[k,2]
      if(discrete){
        g1 <- graphAM(adjMat=A1,edgemode='directed')
        s0 <- multinom.score(xi=xi, dag=g1)/nrow(nbr)
      }
      else{
        if(scoring=='ml')
          s0 <- score$local.score(vertex=j,parents=which(A1[j]!=0))/nrow(nbr)
        else
          s0 <- mvn.score(xi=xi, hyper=hyper, A=A1)/nrow(nbr)
      }
      if(A1[i,j]==1){
        if(runif(n=1)>0.5){
          A1 <- A
          A1[j,i] <- 1    # flip the edge
          A1[i,j] <- 0
          if(!isValidGraph(A1,type='dag')) next
        }
        else A1[i,j] <- 0
      }
      else A1[i,j] <- 1
      nbp <- neighbor(A1)
      if(discrete){
        g1 <- graphAM(adjMat=A1,edgemode='directed')
        s1 <- multinom.score(xi=xi, dag=g1)/nrow(nbr)
      }else{
        if(scoring=='ml')
          s1 <- score$local.score(vertex=j,parents=which(A1[,j]!=0))
        else
          s1 <- mvn.score(xi=xi, hyper=hyper, A=A1)
      }
      s1 <- s1/nrow(nbp)
      accept <- FALSE
      if(s1>s0) accept <- TRUE
      else{
        e <- exp(s1-s0)
        if(runif(n=1)<e) accept <- TRUE
      }
      if(accept) A <- A1
      else s1 <- s0
      if(nrow(nbp)>0) istep <- istep + 1
    }
    else if(method=='gibbs'){
      if(representation=='edge.set'){
        qtable <- select_edges(A=A, xi=xi, q=q, score=score,discrete=discrete)
        if(is.null(qtable$score)) next
        iq <- sample(nrow(qtable$score), size=1, prob=qtable$score[,q+1])
        for(l in seq_len(q)){
          i <- qtable$index[l,1]
          j <- qtable$index[l,2]
          A[i,j] <- qtable$score[iq,l]
        }
      }
      else if(representation=='parent.set'){
        Ak <- select_parents(A=A, xi=xi, q=q, ac=ac, score=score,
                             discrete=discrete,
                             kappa=kappa, verbose=verbose,
                             progress.bar=progress.bar)
        if(is.na(Ak$flag)) next
        A <- Ak$Ak
      }
      else if(representation=='gm'){    # goudie-mukherjee

        Pawgh <- partition.pset(A=A, xi=xi, q=q, ac=ac, path=path, kappa=kappa,
                                progress.bar=progress.bar)
        nh <- length(Pawgh$PaH)
        prob <- exp(Pawgh$KH)
        prob <- prob/sum(prob)
        h <- sample(seq_len(nh), size=1, prob=prob)  # Parition index

        pah <- Pawgh$PaH[[h]]
        W <- names(pah)
        aw <- A
        A1 <- A
        A1[,W] <- 0
        for(w in W){
          prob <- NULL
          np <- length(pah[[w]])
          paw <- list()
          for(j in seq_len(np)){
            pa <- pah[[w]][j]
            paw[[j]] <- nodes[which(ac[,pa]==1)]
            aw[,w] <- 0
            aw[paw[[j]],w] <- 1
            g <- graphAM(adjMat=aw, edgemode='directed')
            sc <- multinom.score(xi=xi, dag=g)
            prob <- c(prob,sc)
          }
          prob <- prob-max(prob)
          prob <- exp(prob)
          prob <- prob/sum(prob)
          hi <- sample(np, size=1,prob=prob)
          A1[paw[[hi]],w] <- 1
        }
        path <- path.update(A0=A, A1=A1, path)
        A <- A1
      }
      istep <- istep + 1
    }
    else stop('Unknown method')

    if(progress.bar) setTxtProgressBar(pb, istep/nstep)
    if(verbose==3){
      d <- distance(m=A, mref=Aref)
      cnt <- cnt + 1
      sumd <- sumd + d
      if(cnt==npr){
        cat('istep = ',istep,', mean distance = ',sumd/cnt,'\n',sep='')
        plot(graphAM(adjMat=Aref,edgemode='directed'), main='True')
        plot(graphAM(adjMat=A,edgemode='directed'), main=paste0('dist=',d))
        sumd <- cnt <- 0
      }
    }
    if(istep>=nstep) break
  }
  if(progress.bar) close(pb)

# dag <- mat2nel(nodes=colnames(xi),A)
  colnames(A) <- colnames(xi)
  dag <- graphAM(adjMat=A, edgemode='directed')

  return(dag)
}

neighbor <- function(A){

  p <- nrow(A)
  nbr <- NULL
  for(i in seq(1,p)) for(j in seq(1,p)){
    if(i==j) next
    Ap <- A
    Ap[i,j] <- 1-A[i,j]
    if(isValidGraph(Ap,type='dag')) nbr <- rbind(nbr,c(i,j))
  }

  return(nbr)
}
