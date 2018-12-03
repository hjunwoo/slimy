#' Main MCMC sampler
#' @param xi Data frame of continuous or discrete data with samples in rows and
#'           nodes in columns
#' @param dag Initial graph. If \code{NULL}, a random graph will be generated
#' @param  ref Reference graph for comparison
#' @param discrete \code{TRUE} if discrete data. If \code{FALSE}, gaussian data
#'             assumed
#' @param nstep Number of steps
#' @param verbose Level of verbosity
#' @param representation Method options for Gibbs sampling;
#'                      \code{c('gm','parent.set','edge.set')}
#' @param q Number of nodes to sample in each step in Gibbs sampling
#' @param npr Periods for collecting printing statistics
#' @param kappa Maximum in-degree per node assumed
#' @param burn.in Initial periods to throw away before collecting statistics
#' @export
mc.sample <- function(object, init.dag=NULL, nstep=1000, verbose=3,
                      kappa=NULL,
                      progress.bar=FALSE, burn.in=100, map=FALSE,
                      q=2, npr=100, nplot=npr, nprime=1,attrs=NULL,
                      init.deg=2, frq.update=1, track.field=FALSE,
                      dmax=3,dy=0.01, init='sample',
                      pois=NULL,
                      ncores=1, update.n=NULL, useC=TRUE){

  type <- object@data.type
  if(burn.in >= nstep) burn.in <- nstep-1
  m <- object@nsample
  p <- object@p
  nodes <- object@nodes
  if(is.null(kappa)) kappa <- object@kappa
  prior <- object@prior
  hyper <- object@hyper

  discrete <- type=='discrete'
  if(is.null(init.dag))   # initial graph
    dag <- rgraph(nodes=nodes, discrete=discrete,
                       mean.degree=init.deg, max.degree=kappa)
  A <- as(dag,'matrix')
  A <- apply(A,1:2,as.numeric)
  colnames(A) <- nodes

  ref <- object@ref.dag
  if(type=='mvln') xi <- object@latent.var

  sumd <- 0
  cnt <- 0

  ac <- object@ac
  cache <- object@cache
  if(sum(dim(ac))==0){
    ac <- parent.sets(nodes=nodes, kappa)
    object@ac <- ac
  }

  if(type=='counts'){
    ci <- object@data
    po <- pois_stat(yi=ci)
    if(track.field) xi0 <- object@latent.var  # reference field
    s0 <- sqrt(diag(po$sigma))
    sg <- po$sigma/outer(s0,s0,'*')
    if(init=='sample')
      xi <- MASS::mvrnorm(n=m, mu=rep(0,p), Sigma=sg)
    else if(init=='input')
      xi <- object@latent.var
    else
      xi <- matrix(rnorm(n=m),nrow=m,ncol=p)
    const <- -sum(lfactorial(ci))
    const <- const + lgamma(hyper$a+0.5*nsample) - 0.5*nsample*log(pi)
    if(hyper$a > 0) const <- const - lgamma(hyper$a)
    if(hyper$b > 0) const <- const + hyper$a*log(2*hyper$b)
    colnames(xi) <- nodes
    object@latent.var <- xi
    po$sigma <- sqrt(diag(po$sigma))
    if(!is.null(pois))
      po <- pois
  }
  else if(sum(dim(cache))==0) # score absent
    cache <- local.score(object, kappa=kappa,
                         progress.bar=progress.bar, ncores=ncores)
  path <- path.count(dag=dag)$C

  istep <- iburned <- 0
  edge.prob <- matrix(0, nrow=p, ncol=p)   # edge.probability
  rownames(edge.prob) <- colnames(edge.prob) <- nodes
  Map <- NULL   # MAP graph and its score
  Elm <- NULL   # Mean log lkh

  while(TRUE){

    W <- sample(nodes,size=q)
    if(type=='counts'){
      if(frq.update > 0) if(istep %% frq.update==0){
        object <- update.field(object, W=W, hyper=hyper,
                               po=po, A=A, dmax=dmax,dy=dy,
                               update.n=update.n, useC=useC)
        xi <- object@latent.var
      }
      object <- local.score(object, kappa=kappa, po=po,
                           progress.bar=progress.bar,
                           ncores=ncores)
      cache <- object@cache
    }

    Pawgh <- partition.pset(A=A, q=q, W=W, ac=ac, path=path,
                      kappa=kappa,cache=cache,
                      progress.bar=progress.bar,
                      ncores=ncores)
    A1 <- sample.subgraph(A=A, Pawgh=Pawgh, ac=ac, cache=cache,
                          path=path)
    path <- path.update(A0=A, A1=A1, path)
    A <- A1

    istep <- istep + 1

    if(verbose>=2){
      if(!sum(dim((ref@adjMat)))==0){
        d <- distance(m=A, mref=as(ref,'matrix'))
        sumd <- sumd + d
      }
      cnt <- cnt + 1

      if(cnt==npr){
        ddag <- graphAM(adjMat=A,edgemode='directed')
        if(discrete){
          if(is.null(cache))
            llk <- multinom.score(xi, dag=ddag, hyper=hyper)
          else
            llk <- multinom.cache.score(dag=ddag, ac, cache)
        }
        else if(type=='counts'){
          llk <- pois.score.global(ci=ci,xi=xi,A=A,hyper=hyper,
                                   po=po, ac=ac)+const
          if(track.field) mse <- mean(c(xi-xi0)^2)
        }
        else if(prior=='g')
          llk <- g.score.global(xi=xi, A=A, hyper=hyper,
                                ac=ac, cache=cache)
        else if(prior=='diag')
          llk <- diag.score.global(xi=xi, A=A, hyper=hyper, ac=ac,
                                   cache=cache)
        else stop('Unknown type/prior')

        cat('istep = ',istep,', log LK = ',llk/m/p,
            ', mean distance = ', sumd/cnt,sep='')
        if(type=='counts' & track.field)
          cat(', MSE(field) = ', mse, sep='')
        cat('\n')
        if(nplot>0){
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
        }
        sumd <- cnt <- 0
        if(istep > burn.in){
          iburned <- iburned + 1
          edge.prob <- edge.prob + A
          if(is.null(Map))
            Map <- list(dag=ddag, score=llk/m/p)
          else if(llk > Map$score)
            Map <- list(dag=ddag, score=llk/m/p)
          Elm <- c(Elm, llk)
        }
      }
    }
    if(istep>=nstep) break
  }
  Elm <- mean(Elm)/m/p

  colnames(A) <- colnames(xi)
  dag <- graphAM(adjMat=A, edgemode='directed')
  object@dag <- dag
  edge.prob <- edge.prob/iburned
  object@edge.prob <- edge.prob
  object@map <- Map
  object@mlk <- llk/m/p
  object@ac <- ac
  object@cache <- cache
  object@emlk <- Elm

  return(object)
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

compute.score <- function(object=NULL, kappa=3, prior='g',
                          hyper=NULL, progress.bar=TRUE,
                          ncores=1){

  object@kappa <- kappa
  object@ac <- parent.sets(nodes=object@nodes, kappa)
  object@prior <- prior
  if(!is.null(hyper)) object@hyper <- hyper

  object <- local.score(object, kappa=kappa,
                        progress.bar=progress.bar,
                        ncores=ncores)

  return(object)
}

# Sample new parent set using factorization result

sample.subgraph <- function(A, Pawgh, ac, cache,path){

  nh <- length(Pawgh$PaH)
  prob <- exp(Pawgh$KH)
  prob <- prob/sum(prob)
  h <- sample(seq_len(nh), size=1, prob=prob)  # Partition index

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

  return(A1)
}
