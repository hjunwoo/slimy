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
local.score <- function(object, kappa, po=NULL, progress.bar, ncores){

  nodes <- object@nodes
  p <- length(nodes)
  ac <- object@ac
  cache <- matrix(0, nrow=p, ncol=ncol(ac))
  rownames(cache) <- nodes

  if(object@data.type != 'counts')
    cat('Computing local scores ...\n')

  bundle <- list(object=object, po=po)
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
    Rmpi::mpi.bcast.Robj2slave(object)
    lcache <- Rmpi::mpi.applyLB(seq_len(nac), FUN=fill.cache,
                                bundle)
  }

  maxs <- -Inf
  for(iac in seq_len(nac))
    for(x in lcache[[iac]]) if(!is.na(x)) if(x>maxs) maxs <- x

  for(iac in seq_len(nac)){
    z <- lcache[[iac]]
    for(i in seq_len(length(z))) if(!is.na(z[i])) z[i] <- z[i] - maxs
    cache[,iac] <- z
  }

  object@cache <- cache
  return(object)
}

fill.cache <- function(iac, bundle){

  object <- bundle$object
  po <- bundle$po

  ac <- object@ac
  hyper <- object@hyper
  prior <- object@prior
  nodes <- object@nodes
  p <- length(nodes)
  type <- object@data.type
  if(type %in% c('counts','mvln')){
    ci <- object@data
    xi <- object@latent.var
  }
  else xi <- object@data


  pa <- nodes[which(ac[,iac]==1)]
  lcache <- double(p)
  for(i in seq_len(p)){
    w <- nodes[i]
    if(w %in% pa) sc <- NA
    else{
      wpa <- nodes[nodes %in% c(w,pa)]
      nw <- length(wpa)
      A <-  matrix(0, nrow=nw, ncol=nw)
      rownames(A) <- colnames(A) <- wpa
      A[pa,w] <- 1
      if(!is.DAG(A)) sc <- NA
      else{
        if(type=='discrete') sc <-
            multinom.local.score(xi, w, pa, prior=prior,
                                 hyper=hyper)
        else if(type=='counts')
          sc <- pois.score(ci=ci,xi=xi,node=w, pa=pa, hyper=hyper,
                           po=po)
        else if(prior=='g')
          sc <- g.score(xi=xi, node=w, pa=pa, hyper=hyper)
        else if(prior=='diag')
          sc <- diag.score(xi=xi, node=w, pa=pa, hyper=hyper)
        else stop('Unknown prior')
      }
    }
    lcache[i] <- sc
  }

  return(lcache)
}
