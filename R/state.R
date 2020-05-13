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
         prototype(nbits=1L, coefs=c(1. + 0i, 0i)),
         validity=check_qstate)

## "constructor" function
qstate <- function(nbits=1L, coefs=c(1+0i, rep(0i, times=2^nbits-1))) {
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

setClass("sqgate",
         representation(bit="integer",
                        M="array"),
         prototype(bit=c(1L), M=array(as.complex(c(1,0,0,1)), dim=c(2,2))))

sqgate <- function(bit=1L, M=array(as.complex(c(1,0,0,1)), dim=c(2,2))) {
  return(new("sqgate", bit=as.integer(bit), M=M))
}

setMethod("*", c("sqgate", "qstate"),
          function(e1, e2) {
            stopifnot(e1@bit > 0 && e1@bit <= e2@nbits)
            bit <- e1@bit
            nbits <- e2@nbits
            ii <- c(0, 2^(bit-1))
            res <- c()
            ## outside
            for(k in c(0:(2^(nbits-bit)-1))) {
              ## inside
              for(i in c(0:(2^(bit-1)-1))) {
                ## the 2x2 matrix
                for(j in c(0,1)) {
                  ll <- j*2^(bit-1) + i + k*2^bit + 1
                  rr <- ii + i + k*2^(bit) + 1
                  res[ll] <-  sum(e1@M[j+1,]*e2@coefs[rr])
                }
              }
            }
            return(qstate(nbits=nbits, coefs=as.complex(res)))
          }
          )


## first bit is control bit
setClass("cnotgate",
         representation(bits="integer"),
         prototype(bits=c(1L,2L)))

cnotgate <- function(bits=c(1, 2)) return(new("cnotgate", bits=as.integer(bits)))
                        
setMethod("*", c("cnotgate", "qstate"),
          function(e1, e2) {
            stopifnot(length(e1@bits) == 2)
            stopifnot(all(e1@bits > 0) && all(e1@bits < e2@nbits))
            stopifnot(e1@bits[1] != e1@bits[2])
            ## control bit == 1
            al <- array(as.integer(c(0:(2^e2@nbits-1))), dim=c(2^e2@nbits,1))
            cb <- apply(al, MARGIN=1, FUN=function(x, y) bitwAnd(x, as.integer(y))>0, y=2^(e1@bits[1]-1))
            ## target bit
            tb <- apply(al, MARGIN=1, FUN=function(x, y) bitwAnd(x, as.integer(y))>0, y=2^(e1@bits[2]-1))
            x <- which(cb & tb)
            y <- which(cb & !tb)
            tmp <- e2@coefs[x]
            e2@coefs[x] <- e2@coefs[y]
            e2@coefs[y] <- tmp
            return(e2)
          }
          )
