#' Multivariate gaussian score
#' @param xi data matrix (rows= samples, columns=variables)
#' @param A adjacency matrix (n x n); A_ij =1 iff xi -> xj
#' @param B adjacency of prior network
#' @param mu0 Prior mean; vector of length n
#' @param v   Conditional variance; vector of length n
#' @export
#'
mvn.score <- function(xi, hyper, A){

  #compute hyperparameters

  m <- nrow(xi)
  N <- ncol(A)
  up <- down <- list()
  for(i in seq_len(N)){
    x <- which(A[,i]!=0)  # parents of i
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
# c1 <- lcna(n, alpha)
# c2 <- lcna(n, alpha+m)
  c12 <- dcna(n, alpha,m)
  dT0 <- as.numeric(determinant(T0)$modulus)
  dTm <- as.numeric(determinant(Tm)$modulus)
# r <- N + c1-c2 + (alpha/2)*dT0 - ((alpha+m)/2)*dTm
  r <- N + c12 + (alpha/2)*dT0 - ((alpha+m)/2)*dTm

  return(r)
}

# compute precision matrix from B
precision <- function(v, B){

  n <- length(v)
  w <- 1/v[1]
  for(i in seq(1,n-1)){
    bip <- B[seq(1,i), i+1]
    wp <- w + outer(bip, bip)/v[i+1]
    wp <- rbind(wp, -t(bip)/v[i+1])
    wp <- cbind(wp, c(-bip/v[i+1],1/v[i+1]))
    w <- wp
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
hyper <- function(xi, nu, alpha, v, mu0, B){

  n <- nrow(B)
  W <- precision(v=v, B=B)
  Sig <- solve(a=W, b=diag(n))
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
