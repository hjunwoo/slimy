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
