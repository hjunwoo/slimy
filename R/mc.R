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
                      nstep=1000, verbose=3, scoring='ml', cache=NULL,
                      progress.bar=FALSE, burn.in=100, map=FALSE,
                      hyper=NULL, q=2, npr=100, nplot=npr, kappa=3,
                      nprime=1,attrs=NULL,
                      g=1e10, ncores=1, useC=FALSE){

  if(burn.in >= nstep) burn.in <- nstep-1
  Big <- 100    # exp(-x)=0 for x > Big
  m <- nrow(xi) # no. of samples
  p <- ncol(xi) # no. of nodes
  nodes <- colnames(xi)

  if(is.null(dag))
      dag <- rgraph(nodes=nodes,discrete=discrete,mean.degree=0.1,
                    max.degree=kappa)
  A <- as(dag,'matrix')
  A <- apply(A,1:2,as.numeric)
  colnames(A) <- colnames(xi)

  sumd <- 0
  cnt <- 0

  ac <- parent.sets(nodes=colnames(xi), kappa)

  if(is.null(cache))
    cache <- local.score(xi=xi, ac=ac, kappa=kappa,
                             discrete=discrete, scoring=scoring,
                             score=score, hyper=hyper, g=g,
                             progress.bar=progress.bar, ncores=ncores)
  path <- path.count(dag=dag)$C

  istep <- iburned <- 0
  edge.prob <- matrix(0, nrow=p, ncol=p)   # edge.probability
  rownames(edge.prob) <- colnames(edge.prob) <- nodes
  Map <- NULL   # MAP graph and its score
  Elm <- NULL   # Mean log lkh

  while(TRUE){

    Pawgh <- partition.pset(A=A, xi=xi, q=q, ac=ac, path=path,
                      kappa=kappa,cache=cache, progress.bar=progress.bar,
                      ncores=ncores, useC=useC)
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

    istep <- istep + 1

    if(verbose>=2){
      if(!is.null(ref)){
        d <- distance(m=A, mref=as(ref,'matrix'))
        cnt <- cnt + 1
        sumd <- sumd + d
      }
      if(cnt==npr){
        ddag <- graphAM(adjMat=A,edgemode='directed')
        if(discrete){
          if(is.null(cache))
            llk <- multinom.score(xi, dag=ddag, scoring=scoring,
                                nprime=nprime)
          else
            llk <- multinom.cache.score(dag=ddag, ac, cache)
        }
        else if(scoring=='ml')
          llk <- score$global.score(as(A,'GaussParDAG'))
        else if(scoring=='bge')
          llk <- mvn.score(xi=xi, hyper=hyper, A=A)
        else if(scoring=='g')
          llk <- g.score.global(xi=xi, A=A, g=g, ac=ac, cache=cache)
        else stop('Unknown scoring')

        cat('istep = ',istep,', log LK = ',llk/m/p,', mean distance = ',
            sumd/cnt,'\n',sep='')
        if(istep %% nplot==0){
          if(!is.null(attrs)){
            plot(ref, main='True',attrs=attrs)
            plot(graphAM(adjMat=A,edgemode='directed'),
                 main=paste0('Distance=',d),attrs=attrs)
          }else{
            plot(ref, main='True')
            plot(graphAM(adjMat=A,edgemode='directed'),
                 main=paste0('Distance=',d))
          }
        }
        sumd <- cnt <- 0
        if(istep > burn.in){
          iburned <- iburned + 1
          edge.prob <- edge.prob + A
          if(is.null(Map))
            Map <- list(dag=ddag, score=llk)
          else if(llk > Map$score)
            Map <- list(dag=ddag, score=llk)
          Elm <- c(Elm, llk)
        }
      }
    }
    if(istep>=nstep) break
  }
  colnames(A) <- colnames(xi)
  dag <- graphAM(adjMat=A, edgemode='directed')
  edge.prob <- edge.prob/iburned
  Elm <- mean(Elm)/m/p
  z <- list(dag=dag, map=Map, edge.prob=edge.prob, elm=Elm)

  return(z)
}

neighbor <- function(A){

  p <- nrow(A)
  nbr <- c(0,0)
  for(i in seq(1,p)) for(j in seq(1,p)){
    if(i==j) next
    if(A[i,j]==0){
      Ap <- A
      Ap[i,j] <- 1
      if(!is.DAG(Ap)) next
    }
    nbr <- rbind(nbr,c(i,j))
  }

  return(nbr)
}

#' Enumerate and store local scores of all parent sets
#' @export

compute.score <- function(xi, kappa=3, discrete=FALSE, scoring='ml',
                          hyper=NULL, g=1e10, progress.bar=TRUE,nprime=1,
                          ncores=1){

  ac <- parent.sets(nodes=colnames(xi), kappa)
  if(!discrete)
    if(scoring=='ml') score <- new('GaussL0penObsScore',xi)

  if(scoring=='bge'){
    if(is.null(hyper)) stop('Hyperparameters for bge score missing')
    B <- hyper$B
    n <- nrow(B)
    W <- precision(v=hyper$v, B=B)
    Sig <- solve(a=W, b=diag(n))
    rownames(Sig) <- colnames(Sig) <- rownames(B)
    hyper <- list(nu=hyper$nu, alpha=hyper$alpha, v=hyper$v,
                mu0=hyper$mu0, Sig=Sig)
  }

  cache <- local.score(xi=xi, ac=ac, kappa=kappa, discrete=discrete,
                       scoring=scoring, score=score, g=g, hyper=hyper,
                       nprime=nprime,
                       progress.bar=progress.bar, ncores=ncores)

  return(cache)
}
