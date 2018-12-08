#' Compute sample mean and variance for poisson log-normal model
#' @export
pois_stat <- function(yi){

  nsample <- nrow(yi)
  nodes <- colnames(yi)
  p <- length(nodes)

  alpha <- colMeans(yi)   # sample mean
  v <- cov(yi)            # covariance matrix
  sigma <- 1 + v/outer(alpha,alpha,'*')
  sigma <- log(sigma)
  sigma <- Matrix::nearPD(x=sigma)$mat # nearest positive definite mat
  sigma <- as.matrix(sigma)
  mu <- log(alpha) - diag(sigma)/2

  names(mu) <- colnames(sigma) <- rownames(sigma) <- nodes
  x <- list(mu=mu, sigma=sigma)

  return(x)
}

pois.score <- function(ci, xi, node, pa, hyper, po){

    nsample <- nrow(ci)
    y <- xi[,node]
    v <- 2*hyper$b + crossprod(y)  # y^t*y
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

update.field <- function(object, W, hyper, po, A, xupdate,
                         update.n,
                         dmax,dy, useC=TRUE){

  if(is.null(update.n)) update.n <- object@nsample
  else update.n <- min(update.n, object@nsample)

  if(useC){
    ci <- object@data
    xi <- object@latent.var
    nodes <- object@nodes
    w <- match(W, nodes)-1   # node ID
    update.n <- c(update.n)
    seed <- c(runif(n=1, max=1000))
    if(xupdate=='gibbs')
      xi2 <- update_field(ci, xi, w, hyper, po, A, dmax, dy,
                          update.n, seed)
    else # Metropolis
      xi2 <- update_metro(ci, xi, w, hyper, po, A, dy,
                          update.n, seed)
    colnames(xi2) <- nodes
    object@latent.var <- xi2
  }
  else
    object <- update_fieldR(object=object, hyper=hyper,
                            W=W, po=po, A=A,
                            dmax=dmax, dy=dy, xupdate=xupdate,
                            update.n=update.n)

  return(object)
}

# sample and update latent field xi
update_fieldR <- function(object, hyper, W, po, A, dmax=3, dy=0.01,
                          xupdate, update.n){

  ci <- object@data
  xi <- object@latent.var

  nsample <- nrow(ci)
  nodes <- colnames(ci)
  p <- ncol(ci)
  grid.y <- seq(-dmax,dmax,by=dy)

  for(w in W){
    pa <- nodes[which(A[,w]==1)]
    npa <- length(pa)
    yx <- xi[,w]
    xpa <- as.matrix(xi[,pa])
    if(npa>0)
      xtx <- solve(diag(npa)+hyper$v*crossprod(xpa))
    else xtx <- NULL
    m=po$mu[w]
    sg=po$sig[w]
    cw=ci[,w]
    for(k in sample(nsample,size=update.n)){
      if(xupdate=='gibbs'){
        prob <- NULL
        for(y in grid.y){
          yx[k] <- y
          x <- zprob(yx=yx,hyper=hyper,xpa=xpa,xtx=xtx,nsample=nsample,
                       cw=cw,m=m,sg=sg)
          prob <- c(prob,x)
        }
        prob <- prob - max(prob)
        prob <- exp(prob)
        prob <- prob/sum(prob)
        ys <- sample(grid.y, size=1, prob=prob)
        xi[k,w] <- ys
      }
      else{    # Metropolis-Hastings
        x0 <- zprob(yx=yx,hyper=hyper,xpa=xpa,xtx=xtx,nsample=nsample,
                    cw=cw,m=m,sg=sg)
        delta <- rnorm(n=1,mean=0, sd=dy)
        yx[k] <- yx[k] + delta
        x1 <- zprob(yx=yx,hyper=hyper,xpa=xpa,xtx=xtx,nsample=nsample,
                    cw=cw,m=m,sg=sg)
        accept <- FALSE
        if(x1>x0) accept <- TRUE
        else{
          prob <- exp(x1-x0)
          if(prob>runif(n=1)) accep <- TRUE
        }
        if(accept)
          xi[k,w] <- yx[k]
        else
          yx[k] <- yx[k] - delta
      }
    }
  }

  object@latent.var <- xi
  return(object)
}

zprob <- function(yx,hyper,xpa,xtx,nsample,cw,m,sg){

  v <- 2*hyper$b + crossprod(yx)
  npa <- ncol(xpa)
  if(npa>0){
    xy <- crossprod(xpa,yx)
    v <- v - hyper$v*t(xy) %*% xtx %*% xy
  }
  z <- -(0.5*nsample+hyper$a)*log(v)
  if(npa>0)
    z <- z + 0.5*determinant(xtx,log=TRUE)$modulus
  lkh <- sum((sg*yx+m)*cw-exp(sg*yx+m))
  x <- as.numeric(z + lkh)

  return(x)
}

pois.score.global <- function(ci, xi, A, hyper, po, ac){

  nodes <- colnames(A)
  llk <- 0
  llk2 <- 0
  for(w in nodes){
    pa <- nodes[which(A[,w]!=0)]
    llk <- llk + pois.score(ci=ci,xi=xi, node=w, pa=pa, hyper=hyper,
                              po=po)
  }
  return(llk)
}
