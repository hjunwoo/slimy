#' Class \code{slimySet} for storing data and graph inference results
#'
#' @slot nodes Vector of node names
#' @slot p     Number of nodes
#' @slot data.type \code{'continuous'} for Gaussian-type data,
#'                 \code{'discrete'} for multinomial-type data,
#'                 or \code{'counts'} for Poisson-type count data.
#' @slot data  Data frame containing data with samples in rows and
#'             nodes in columns
#' @slot latent.var Latent variable associated with data
#'             for \code{counts}.
#' @slot nsample Number of samples
#' @slot ref.dag Reference (true) graph underlying data, if known.
#' @slot prior Type of prior used.
#'            \code{'g'} for g-prior,
#'            \code{'multinomial'} for multinomial,
#'            \code{'diagonal'} for diagonal prior.
#' @slot hyper List of hyperparameters.
#' @slot kappa Maximum in-degree.
#' @slot ac Matrix of all possible parent sets.
#' @slot cache Local scores. Same dimension as \code{ac}.
#' @slot dag Graph inferred at the moment
#' @slot mlk Log marginal likelihood per node per sample
#' @slot edge.prob Posterior edge probability matrix
#' @slot map Maximum a posterior graph; list of components
#'           \code{dag}, \code{graphAM} object of the MAP graph,
#'           and \code{score}, the corresponding \code{mlk}.
#' @return Object of class \code{slimySet}
#' @export slimySet
setClass('slimySet',
         slots=c(nodes='vector',
                 p='numeric',
                 data.type='character',
                 data='matrix',
                 nsample='numeric',
                 latent.var='matrix',
                 ref.dag='graphAM',
                 prior='character',
                 hyper='list',
                 kappa='numeric',
                 ac='matrix',
                 cache='matrix',
                 dag='graphAM',
                 mlk='numeric',
                 emlk='numeric',
                 edge.prob='matrix',
                 map='list')
)

#' Create \code{slimySet} object
#'
#' @param data.type Data type
#' @param data Input data
#' @return Object of class \code{slimySet}.
#' @export
slimySet <- function(data.type, data, ref.dag=NULL, hyper=NULL){

  if(class(data)!='matrix') data <- as.matrix(data)
  x <- new('slimySet',data.type=data.type,data=data)
  x@nodes <- colnames(data)
  x@p <- ncol(data)
  x@nsample <- nrow(data)
  if(is.null(ref.dag)) x@ref.dag <- graphAM()
  else x@ref.dag <- ref.dag
  if(is.null(hyper)) x@hyper <- list()
  else x@hyper <- hyper

  if(data.type=='mvln'){
    musig <- pois_stat(data)
    mu <- musig$mu
    sig <- diag(musig$sigma)
    xi <- (log(1+data) - mu)/sig
    x@latent.var <- xi
  }

  return(x)
}

#' Dispaly object
#'
#' @param object Object of class \code{slimySet}
#' @return \code{NULL}
#' @export
setMethod('show', signature='slimySet',
          definition=function(object){
            cat('Class:', class(object),'\n')
            cat('Dim:', object@p,'\n')
            cat('Data type:', object@data.type,'\n')
            cat('Sample size:', object@nsample,'\n')
            cat('Nodes: ')
            print(head(object@nodes))
          })
