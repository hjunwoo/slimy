#' compute the g-prior score
#' @export
g.score <- function(xi, node, pa, hyper){

  nsample <- nrow(xi)
  npa <- length(pa)
  g <- hyper$g
  if('a' %in% names(hyper)) a <- hyper$a
  else a <- 0
  if('b' %in% names(hyper)) b <- hyper$b
  else b <- 0

  y <- xi[,node]
  v <- 2*b + crossprod(y)
  if(npa >= 1){
    xpa <- as.matrix(xi[,pa])
    xtx <- solve(crossprod(xpa))
    fit <- xtx %*% crossprod(xpa,y)
    v <- v - g/(1+g)*crossprod(xpa %*% fit)
  }
  z <- -0.5*npa*log(1+g)-(0.5*nsample+a)*log(v)
  z <- z + a*log(2) + lgamma(a+nsample/2) - nsample*log(pi)/2
  if(a*b >0) z <- z + a*log(b) - lgamma(a)

  return(z)
}

g.score.global <- function(xi, A, hyper, ac=NULL, cache=NULL){

  nodes <- colnames(A)
  llk <- 0
  for(w in nodes){
    pa <- nodes[which(A[,w]!=0)]
    if(is.null(cache))
      llk <- llk + g.score(xi=xi, node=w, pa=pa, hyper=hyperb )
    else{
      z <- apply(ac,2,function(x){sum(x!=A[,w])})
      k <- which(z==0)
      llk <- llk + cache[w,k]
    }
  }

  return(llk)
}


#' Compute the diagonal prior score
#' @export
diag.score <- function(xi, node, pa, hyper){

  nsample <- nrow(xi)
  npa <- length(pa)

  y <- xi[,node]
  v <- 2*hyper$b + crossprod(y-mean(y))  # y^t*y
  if(npa >= 1){
    xpa <- as.matrix(xi[,pa])
    xtx <- solve(diag(npa)+hyper$v*crossprod(xpa))
    xy <- crossprod(xpa,y)
    v <- v - hyper$v*t(xy) %*% xtx %*% xy
  }
  z <- -(0.5*nsample+hyper$a)*log(v)
  if(npa>0)
    z <- z + 0.5*determinant(xtx,log=TRUE)$modulus

  return(z)
}

diag.score.global <- function(xi, A, hyper, ac=NULL, cache=NULL){

  nodes <- colnames(A)
  llk <- 0
  for(w in nodes){
    pa <- nodes[which(A[,w]!=0)]
    if(is.null(cache))
      llk <- llk + diag.score(xi=xi, node=w, pa=pa, hyper=hyper)
    else{
      z <- apply(ac,2,function(x){sum(x!=A[,w])})
      k <- which(z==0)
      llk <- llk + cache[w,k]
    }
  }

  return(llk)
}

#' Multivariate gaussian score
#' @param xi data matrix (rows= samples, columns=variables)
#' @param A adjacency matrix (n x n); A_ij =1 iff xi -> xj
#' @param B adjacency of prior network
#' @param mu0 Prior mean; vector of length n
#' @param v   Conditional variance; vector of length n
#' @export
#'
mvn.score <- function(xi, hyper.par, A){

  #compute hyperparameters

  m <- nrow(xi)
  N <- ncol(A)
  nodes <- rownames(A)
  up <- down <- list()
  for(i in nodes){
    x <- nodes[which(A[,i]!=0)]  # parents of i
    down[[i]] <- x
    up[[i]] <- c(x,i)
  }

  sc <- 0
  for(i in seq_len(N)){
    idx <- up[[i]]
    Ti <- as.matrix(hyper$T0[idx,idx])
    Tmi <- as.matrix(hyper$Tm[idx,idx])
    r <- rho(n=length(idx), nu=hyper$nu, alpha=hyper$alpha, m=m,
             T0=Ti, Tm=Tmi)
    sc <- sc + r

    idx <- down[[i]]
    if(length(idx)==0) next
    Ti <- as.matrix(hyper$T0[idx,idx])
    Tmi <- as.matrix(hyper$Tm[idx,idx])
    r <- rho(n=length(idx), nu=hyper$nu, alpha=hyper$alpha, m=m,
             T0=Ti, Tm=Tmi)
    sc <- sc - r
  }
  return(sc)
}

rho <- function(n, nu, alpha, m, T0, Tm){

  N <- (n/2)*log(nu/(nu+m)) - (n*m/2)*log(2*pi)
  c12 <- dcna(n, alpha,m)
  dT0 <- as.numeric(determinant(T0)$modulus)
  dTm <- as.numeric(determinant(Tm)$modulus)
  r <- N + c12 + (alpha/2)*dT0 - ((alpha+m)/2)*dTm

  return(r)
}

# compute precision matrix from B
precision <- function(v, B){

  n <- length(v)
  w <- 1/v[1]
  if(n>1){
    for(i in seq(1,n-1)){
      bip <- B[seq(1,i), i+1]
      wp <- w + outer(bip, bip)/v[i+1]
      wp <- rbind(wp, -t(bip)/v[i+1])
      wp <- cbind(wp, c(-bip/v[i+1],1/v[i+1]))
      w <- wp
    }
  }
  return(w)
}

lcna <- function(n, alpha){

  x <- 0
  for(i in seq_len(n))
    x <- x + lgamma((alpha+1-i)/2)
  x <- x + (alpha*n/2)*log(2) + (n*(n-1)/4)*log(pi)
  x <- -x
  return(x)
}

dcna <- function(n, alpha, m){

  x <- 0
  for(i in seq_len(n))
    x <- x + lgamma((alpha+m+1-i)/2) - lgamma((alpha+1-i)/2)
  x <- x + (m*n/2)*log(2)

  return(x)
}

#' @export
hyper.par <- function(xi, nu, alpha, v, mu0, Sig){

   n <- nrow(Sig)
#  W <- precision(v=v, B=B)
#  Sig <- solve(a=W, b=diag(n))
#  rownames(Sig) <- colnames(Sig) <- rownames(B)
  T0 <- Sig*nu*(alpha-n-1)/(nu+1)
  Xm <- colMeans(xi)
  m <- nrow(xi)
  mul <- (nu*mu0+m*Xm)/(nu+m)
  dxi <- as.matrix(xi) - matrix(Xm, nrow=m, ncol=n, byrow=TRUE)
  Sm <- t(dxi) %*% dxi
  Tm <- T0 + Sm + nu*m/(nu+m)*outer(mu0-Xm, mu0-Xm, '*')

  h <- list(nu=nu,alpha=alpha,v=v,mu0=mu0,T0=T0,Tm=Tm)
  return(h)
}

#' compute the normal-inv-gamma-prior score
#' @export
ninvg.score <- function(xi, node, pa, hyper){

  npa <- length(pa)
  if(npa == 0) return(0)

  nsample <- nrow(xi)
  y <- xi[,node]
  xpa <- as.matrix(xi[,pa])
  dy <- y-hyper$mu*rowSums(xpa)
  z1 <- sum(dy^2)
  z2 <- sum((t(xpa) %*% dy)^2)
  z <- -(hyper$a+0.5)*log(1 + (z1+z2)/2/hyper$b)

  return(z)
}

ninvg.score.global <- function(xi, A, hyper, ac=NULL, cache=NULL){

  nodes <- colnames(A)
  llk <- 0
  for(w in nodes){
    pa <- nodes[which(A[,w]!=0)]
    if(is.null(cache))
      llk <- llk + ninvg.score(xi=xi, node=w, pa=pa, hyper=hyper)
    else{
      z <- apply(ac,2,function(x){sum(x!=A[,w])})
      k <- which(z==0)
      llk <- llk + cache[w,k]
    }
  }

  return(llk)
}
