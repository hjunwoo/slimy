#' Main MCMC sampler
#' @param xi Data frame of continuous or discrete data with samples in rows and
#'           nodes in columns
#' @param dag Initial graph. If \code{NULL}, a random graph will be generated
#' @param  ref Reference graph for comparison
#' @param discrete \code{TRUE} if discrete data. If \code{FALSE}, gaussian data
#'             assumed
#' @param nstep Number of steps
#' @param Big Parameter for maximum numbers in exponents
#' @param verbose Level of verbosity
#' @param representation Method options for Gibbs sampling;
#'                      \code{c('gm','parent.set','edge.set')}
#' @param q Number of nodes to sample in each step in Gibbs sampling
#' @param npr Periods for collecting printing statistics
#' @param kappa Maximum in-degree per node assumed
#' @param burn.in Initial periods to throw away before collecting statistics
#' @export
mc.sample <- function(xi, dag=NULL, ref=NULL,
                      discrete=FALSE,
                      nstep=1000, Big=100, verbose=3, method='gibbs',
                      scoring='ml', representation='gm',
                      progress.bar=FALSE, burn.in=100,
                      hyper=NULL, q=2, npr=100, nplot=npr, kappa=3){

  m <- nrow(xi) # no. of samples
  p <- ncol(xi) # no. of nodes
  nodes <- colnames(xi)

  if(is.null(dag)){
    if(discrete)
      dag <- rgraph(nodes=nodes,mean.degree=1.0)
    else
      dag <- rgraph.gauss(p=p, prob=.1)
  }
  A <- as(dag,'matrix')
  A <- apply(A,1:2,as.numeric)
  colnames(A) <- colnames(xi)

  sumd <- 0
  cnt <- 0

  if(discrete)
    s1 <- sc <- multinom.score(xi=xi, dag=dag)
  else{
    if(scoring=='ml') score <- new('GaussL0penObsScore',xi)
    else s1 <- sc <- mvn.score(xi=xi, hyper=hyper, A=A)
  }

  if(method=='gibbs'){
    ac <- parent.sets(nodes=colnames(xi), kappa)
    if(representation=='gm'){
      cache <- local.score(xi=xi, ac=ac, kappa=kappa, discrete=discrete,
                           score=score, progress.bar=progress.bar)
      path <- path.count(dag=dag)$C
    }
  }

  istep <- iburned <- 0
  edge.prob <- matrix(0, nrow=p, ncol=p)   # edge.probability
  rownames(edge.prob) <- colnames(edge.prob) <- nodes
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
        qtable <- select_edges(A=A, xi=xi, q=q, score=score,
                               discrete=discrete)
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

        Pawgh <- partition.pset(A=A, xi=xi, q=q, ac=ac, path=path,
                      kappa=kappa,cache=cache, progress.bar=progress.bar)
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
            sc <- cache[w,pa]
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

    if(verbose==2){
      if(!is.null(ref)){
        d <- distance(m=A, mref=as(ref,'matrix'))
        cnt <- cnt + 1
        sumd <- sumd + d
      }
      if(cnt==npr){
        ddag <- graphAM(adjMat=A,edgemode='directed')
        if(discrete)
          llk <- multinom.score(xi, dag=ddag)
        else
          llk <- score$global.score(as(A,'GaussParDAG'))
        cat('istep = ',istep,', log LK = ',llk,', mean distance = ',
            sumd/cnt,sep='')
        if(is.null(ref) | verbose < 3) cat('\n')
        if(istep %% nplot==0){
          plot(ref, main='True')
          plot(graphAM(adjMat=A,edgemode='directed'), main=paste0('Distance=',d))
        }
        sumd <- cnt <- 0
        if(istep>burn.in){
          iburned <- iburned + 1
          edge.prob <- edge.prob + A
        }
      }
    }
    if(istep>=nstep) break
  }
  colnames(A) <- colnames(xi)
  dag <- graphAM(adjMat=A, edgemode='directed')
  edge.prob <- edge.prob/iburned

  return(list(dag=dag, edge.prob=edge.prob))
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
