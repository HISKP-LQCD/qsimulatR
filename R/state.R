check_qstate <- function(object) {
  stopifnot(object@nbits > 0)
  stopifnot(object@nbits <= 16)
  N <- 2^object@nbits
  stopifnot(N == length(object@coefs))
}

genStateString <- function(int, nbits) {
  x <- intToBits(int)
  str <- "|"
  for(i in c(nbits:1)) {
    bitstr <- "0"
    if(x[i] > 0) bitstr <- "1"
    str <- paste0(str, bitstr)
  }
  str <- paste0(str, ">")
  return(str)
}

qstatecoefs <- function(y) {
  if (is(y, 'qstate')) y@coefs
  else y
}

eps <- 1.e-12
#' The qstate class
#'
#' This class represents a quantum state
#'
#' @slot nbits The number of qbits
#' @slot coefs The 2^nbits complex valued vector of coefficients
#' 
#' @name qstate
#' @rdname qstate
#' @aliases qstate-class
#' @exportClass qstate
setClass("qstate",
         representation(nbits="integer",
                        coefs="complex"),
         prototype(nbits=1L, coefs=c(1. + 0*1i, 0*1i)),
         validity=check_qstate)

## "constructor" function
qstate <- function(nbits=1L, coefs=c(1+0*1i, rep(0*1i, times=2^nbits-1))) {
  return(new("qstate", nbits=as.integer(nbits), coefs=as.complex(coefs)))
}

setMethod("show", signature(object = "qstate"),
          function(object) {
            N <- 2^object@nbits
            for(i in as.integer(c(1:2^object@nbits))) {
              if(abs(object@coefs[i]) > eps) {
                if(i != 1) cat("+ ")
                else cat("  ")
                cat(object@coefs[i], genStateString(int=i-1L, nbits=object@nbits), "\n")
              }
            }
          }
          )

setMethod("*", c("matrix", "qstate"),
          function(e1, e2) {
            e2@coefs = drop(e1 %*% e2@coefs)
            validObject(e2)
            return(e2)
          }
          )

setClass("qgate",
         representation(bits="integer"),
         prototype(bits=c(1L)))

#' @export
applygate <- function(bit=1L, qstate) {
  stopifnot(bit > 0 && bit <= qstate@nbits)
  nbits <- qstate@nbits
  N  <- 2^nbits
  mcf <- array(c(1,1,1,-1), dim=c(2,2))/sqrt(2)
  ii <- seq(0, 2^bit-1, by=2^(bit-1))
  res <- c()
  for(i in c(0:(N/2-1))) {
    for(j in c(0,1)) {
      res[j*2^(bit-1) + i*2^(nbits-bit) + 1] <-  sum(mcf[j+1,]*qstate@coefs[ii + i*2^(nbits-bit) + 1])
    }
  }
  cat(res, "\n")
  return(qstate(nbits=nbits, coefs=as.complex(res)))
}

