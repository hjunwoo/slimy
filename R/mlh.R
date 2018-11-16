#' Multiple MC sampling with varying hyperparameters
#' @export
margL.scan <- function(xi, ref=NULL, hyper.range, discrete=FALSE,
                       nstep=1000, verbose=3, burn.in=100,
                       progress.bar=FALSE, attrs=NULL,
                       scoring='K2', cut=0.5,
                       q=3, kappa=3, npr=100, nplot=npr, ncores=1){

  nrun <- length(hyper.range)
  elv <- NULL
  elx <- vector('list',nrun)
  for(irun in seq_len(nrun)){
    if(discrete){
      nprime <- hyper.range[irun]
      cat("Sampling under N' = ",nprime,'\n',sep='')
      g=NULL
    }
    else{
      g <- hyper.range[irun]
      cat('Sampling under g = ',g,'\n',sep='')
      nprime=NULL
    }
    mc <- mc.sample(xi=xi, ref=ref, discrete=discrete,
              nstep=nstep, verbose=verbose, scoring=scoring,  g=g,
              burn.in=burn.in, q=q, npr=npr, nplot=nplot,
              progress.bar=progress.bar,
              kappa=kappa, nprime=nprime, attrs=attrs, ncores=ncores)
    cat('Mean log marginal L = ',mc$elm,'\n\n',sep='')
    map.dist <- distance(mc$map$dag, ref)
    post <- posterior.graph(mc$edge.prob,cut=cut)
    mean.dist <- distance(post, ref)
    if(discrete)
      x <- data.frame(nprime=nprime, lmarg=mc$elm, map.dist=map.dist,
                    mean.dist=mean.dist)
    else
      x <- data.frame(g=g, lmarg=mc$elm, map.dist=map.dist,
                    mean.dist=mean.dist)
    elv <- rbind(elv,x)
    elx[[irun]] <- list(dag=mc$dag, map=mc$map, edge=mc$edge.prob)
  }

  z <- list(elv=elv, elx=elx)
  return(z)
}

#' Run multiple chains with varying initial graphs
#' @export
mc.multichain <- function(nchain=2, xi, ref=NULL, discrete=FALSE,
                      nstep=1000, verbose=3, scoring='ml', cache=NULL,
                      progress.bar=FALSE, burn.in=100, map=FALSE,
                      hyper=NULL, q=2, npr=100, nplot=npr, kappa=3,
                      nprime=1,attrs=NULL, init.deg=0,
                      g=1e10, ncores=1){

  if(is.null(cache)){
    ac <- parent.sets(nodes=colnames(xi), kappa=kappa)
    cache <- local.score(xi=xi, ac=ac, kappa=kappa,
                       discrete=discrete, scoring=scoring,
                       score=score, hyper=hyper, g=g,
                       progress.bar=progress.bar, ncores=ncores)
  }

  for(k in seq_len(nchain)){
    cat('Chain #', k, ': \n')
    mc <- mc.sample(xi, dag=NULL, ref=ref,
             discrete=discrete, nstep=nstep, verbose=verbose,
             scoring=scoring, cache=cache, progress.bar=progress.bar,
             burn.in=burn.in, map=map, hyper=hyper, q=q, npr=npr,
             nplot=nplot, kappa=kappa, nprime=nprime, attrs=attrs,
             init.deg=init.deg, g=g, ncores=ncores)
    cat('Mean log LK = ',mc$elm,'\n\n',sep='')
    if(k==1) max.chain <- mc
    else if(mc$elm > max.chain$elm)
      max.chain <- mc
  }
  cat('Max mean log LK = ',max.chain$elm,'\n\n',sep='')
  return(max.chain)
}
