#' Compute sample mean and variance for poisson log-normal model
#' @export
pois_stat <- function(yi){

  nsample <- nrow(yi)
  nodes <- colnames(yi)
  p <- length(nodes)

  alpha <- colMeans(yi)   # sample mean
  v <- cov(yi)            # covariance matrix
# sigma <- 1 + v/outer(alpha,alpha,'*')
#  sigma <- log(sigma)
#  mu <- log(alpha) - diag(sigma)/2
  sigma <- log(1 + (diag(v)-alpha)/alpha^2)
  mu <- log(alpha) - sigma/2
  names(mu) <- names(sigma) <- nodes
  x <- list(mu=mu, sigma=sigma)

  return(x)
}

pois.score <- function(ci, xi, node, pa, hyper, po){

    nsample <- nrow(ci)
    y <- xi[,node]
    v <- 2*hyper$b + crossprod(y)
    npa <- length(pa)
    if(npa > 0){
      xpa <- as.matrix(xi[,pa])
      xtx <- solve(diag(npa)+hyper$v*crossprod(xpa))
      xy <- crossprod(xpa,y)
      v <- v - hyper$v*t(xy) %*% xtx %*% xy
    }
    z <- -(0.5*nsample+hyper$a)*log(v)
    if(npa>0)
      z <- z + 0.5*determinant(xtx,log=TRUE)$modulus

    m <- po$mu[node]
    sg <- po$sigma[node]
    score <- z + sum((sg*y+m)*ci[,node]-exp(sg*y+m))

    return(score)
}

# sample and update latent field xi
update.field <- function(ci, xi, hyper, po, A, dmax=5, dy=0.01){

  nsample <- nrow(ci)
  nodes <- colnames(ci)
  p <- ncol(ci)
  grid.y <- seq(-dmax,dmax,by=dy)

  for(w in nodes){
    pa <- nodes[which(A[,w]==1)]
    npa <- length(pa)
    yx <- xi[,w]
#    for(k in seq_len(nsample)){
    for(k in sample(nsample,size=10)){
      prob <- NULL
      for(y in grid.y){
        yx[k] <- y
        v <- 2*hyper$b + crossprod(yx)
        if(npa > 0){
          xpa <- as.matrix(xi[,pa])
          xtx <- solve(diag(npa)+hyper$v*crossprod(xpa))
          xy <- crossprod(xpa,yx)
          v <- v - hyper$v*t(xy) %*% xtx %*% xy
        }
        if(v<=0) browser()
        z <- -(0.5*nsample+hyper$a)*log(v)
        if(npa>0)
          z <- z + 0.5*determinant(xtx,log=TRUE)$modulus
        m <- po$mu[w]
        s <- po$sigma[w]
        lkh <- sum((s*yx+m)*ci[,w]-exp(s*yx+m))
        prob <- c(prob,z + lkh)
      }
      prob <- prob - max(prob)
      prob <- exp(prob)
      prob <- prob/sum(prob)
      ys <- sample(grid.y, size=1, prob=prob)
      xi[k,w] <- ys
    }
#   print(w)
  }
  return(xi)
}

pois.score.global <- function(ci, xi, A, hyper, po,
                              ac=NULL, cache=NULL){

  nodes <- colnames(A)
  llk <- 0
  for(w in nodes){
    pa <- nodes[which(A[,w]!=0)]
    if(is.null(cache))
      llk <- llk + pois.score(ci=ci,xi=xi, node=w, pa=pa, hyper=hyper,
                              po=po)
    else{
      z <- apply(ac,2,function(x){sum(x!=A[,w])})
      k <- which(z==0)
      llk <- llk + cache[w,k]
    }
  }

  return(llk)
}
