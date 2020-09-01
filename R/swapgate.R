#' The SWAP gate
#'
#' This class represents a generic SWAP gate
#'
#' @slot bits Integer vector of length 2. The two bits to swap.
#'
#' @include state.R
#' 
#' @examples
#' x <- H(1) * qstate(nbits=2)
#' z <- SWAP(c(1,2)) * (H(1) * x)
#'
#' @name swapgate
#' @rdname swapgate
#' @aliases swapgate-class
#' @exportClass swapgate
setClass("swapgate",
         representation(bits="integer"),
         prototype(bits=c(1L,2L)))

#' @export
swapgate <- function(bits=c(1, 2)) return(methods::new("swapgate", bits=as.integer(bits)))

#' The SWAP gate
#'
#' @param bits integer vector of length two, containing the bits to swap.
#'
#' @return
#' An S4 class 'swapgate' object is returned
#' @export
SWAP <- function(bits=c(1, 2)) return(methods::new("swapgate", bits=as.integer(bits)))

#' times-swapgate-qstate
#'
#' Applies a SWAP gate to a quantum state.
#'
#' @param e1 object of S4 class 'swapgate'
#' @param e2 object of S4 class 'qstate'
#'
#' @return
#' An object of S4 class 'qstate'
setMethod("*", c("swapgate", "qstate"),
          function(e1, e2) {
            stopifnot(length(e1@bits) == 2)
            stopifnot(all(e1@bits > 0) && all(e1@bits <= e2@nbits))
            stopifnot(e1@bits[1] != e1@bits[2])
            ## control bit == 1
            al <- 0:(2^e2@nbits-1)
            b1 <- is.bitset(al, bit=e1@bits[1])
            b2 <- is.bitset(al, bit=e1@bits[2])
            x <- which(b1 & !b2)
            y <- which(!b1 & b2)
            e2@coefs[c(x,y)]  <- e2@coefs[c(y,x)]
            ## again the circuit needs extension for plotting
            ngates <- length(e2@circuit$gatelist)
            e2@circuit$gatelist[[ngates+1]] <- list(type="SWAP", bits=c(e1@bits, NA))

            return(e2)
          }
          )


#' The CSWAP gate
#'
#' This class represents a generic SWAP gate, also called Fredkin gate
#'
#' @slot bits Integer vector of length 2. First two bits are the control bits,
#'            third the target bit.
#'
#' @examples
#' x <- qstate(nbits=3)
#' z <- CSWAP(c(1,2,3)) * (H(1) * x)
#'
#' @name cswapgate
#'
#' @rdname cswapgate
#' @aliases cswapgate-class
#' @exportClass cswapgate
setClass("cswapgate",
         representation(bits="integer"),
         prototype(bits=c(1L, 2L, 3L)))

#' @export
cswapgate <- function(bits=c(1, 2)) return(methods::new("cswapgate", bits=as.integer(bits)))

#' The CSWAP or Fredkin gate
#'
#' @param bits integer vector of length two, the first bit being the control and the second
#' the target bit.
#'
#' @aliases fredkin
#' @return
#' An S4 class 'cswapgate' object is returned
#' @export
CSWAP <- function(bits=c(1, 2, 3)) return(methods::new("cswapgate", bits=as.integer(bits)))

#' @rdname CSWAP
#' @aliases CSWAP
#' @export
fredkin <- function(bits=c(1, 2, 3)) return(methods::new("cswapgate", bits=as.integer(bits)))

#' times-cswapgate-qstate
#'
#' Applies a CSWAP gate to a quantum state.
#'
#' @param e1 object of S4 class 'cswapgate'
#' @param e2 object of S4 class 'qstate'
#'
#' @return
#' An object of S4 class 'qstate'
setMethod("*", c("cswapgate", "qstate"),
          function(e1, e2) {
            stopifnot(length(e1@bits) == 3)
            stopifnot(length(e1@bits) == length(unique(e1@bits)))
            stopifnot(all(e1@bits > 0) && all(e1@bits <= e2@nbits))
            ## control bit == 1
            al <- 0:(2^e2@nbits-1)
            cb <- is.bitset(al, bit=e1@bits[1]) 
            b1 <- is.bitset(al, bit=e1@bits[2])
            b2 <- is.bitset(al, bit=e1@bits[3])
            x <- which(cb & (b1 & !b2))
            y <- which(cb & (!b1 & b2))
            e2@coefs[c(x,y)]  <- e2@coefs[c(y,x)]
            ## again the circuit needs extension for plotting
            ngates <- length(e2@circuit$gatelist)
            e2@circuit$gatelist[[ngates+1]] <- list(type="CSWAP", bits=c(e1@bits))

            return(e2)
          }
          )
