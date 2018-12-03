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

update.field <- function(object, W, hyper, po, A, update.n,
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
    xi2 <- update_field(ci, xi, w, hyper, po, A, dmax, dy, update.n,
                        seed)
    colnames(xi2) <- nodes
    object@latent.var <- xi2
  }
  else
    object <- update_fieldR(object=object, W=W, po=po, A=A,
                            update.n=update.n)

  return(object)
}

# sample and update latent field xi
update_fieldR <- function(object, W, po, A, dmax=3, dy=0.01,
                         update.n){

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
    for(k in sample(nsample,size=update.n)){
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
        z <- -(0.5*nsample+hyper$a)*log(v)
        if(npa>0)
          z <- z + 0.5*determinant(xtx,log=TRUE)$modulus
        m <- po$mu[w]
        sg <- po$sigma[w]
        lkh <- sum((sg*yx+m)*ci[,w]-exp(sg*yx+m))
        prob <- c(prob,z + lkh)
      }
      prob <- prob - max(prob)
      prob <- exp(prob)
      prob <- prob/sum(prob)
      ys <- sample(grid.y, size=1, prob=prob)
      xi[k,w] <- ys
    }
  }

  object@latent.var <- xi
  return(object)
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
