#' The CNOT gate
#'
#' This class represents a generic CNOT gate
#'
#' @slot bits Integer vector of length 2. First bit is the control bit,
#'            second the target bit.
#'
#' @include state.R
#' 
#' @examples
#' x <- qstate(nbits=2)
#' ## A Bell state
#' z <- CNOT(c(1,2)) * (H(1) * x)
#'
#' @name cnotgate
#' @rdname cnotgate
#' @aliases cnotgate-class
#' @exportClass cnotgate
setClass("cnotgate",
         representation(bits="integer"),
         prototype(bits=c(1L,2L)))

#' @export
cnotgate <- function(bits=c(1, 2)) return(methods::new("cnotgate", bits=as.integer(bits)))

#' The CNOT gate
#'
#' @param bits integer vector of length two, the first bit being the control and the second
#' the target bit.
#'
#' @return
#' An S4 class 'cnotgate' object is returned
#' @export
CNOT <- function(bits=c(1, 2)) return(methods::new("cnotgate", bits=as.integer(bits)))

#' times-cnotgate-qstate
#'
#' Applies a CNOT gate to a quantum state.
#'
#' @param e1 object of S4 class 'cnotgate'
#' @param e2 object of S4 class 'qstate'
#'
#' @return
#' An object of S4 class 'qstate'
setMethod("*", c("cnotgate", "qstate"),
          function(e1, e2) {
            stopifnot(length(e1@bits) == 2)
            stopifnot(all(e1@bits > 0) && all(e1@bits <= e2@nbits))
            stopifnot(e1@bits[1] != e1@bits[2])
            ## control bit == 1
            al <- 0:(2^e2@nbits-1)
            cb <- is.bitset(al, bit=e1@bits[1])
            ## target bit
            tb <- is.bitset(al, bit=e1@bits[2])
            x <- which(cb & tb)
            y <- which(cb & !tb)
            e2@coefs[c(x,y)]  <- e2@coefs[c(y,x)]
            ## again the circuit needs extension for plotting
            ngates <- length(e2@circuit$gatelist)
            e2@circuit$gatelist[[ngates+1]] <- list(type="CNOT", bits=c(e1@bits, NA))

            return(e2)
          }
          )


#' The CCNOT gate
#'
#' This class represents a generic CNOT gate
#'
#' @slot bits Integer vector of length 2. First two bits are the control bits,
#'            third the target bit.
#'
#' @examples
#' x <- qstate(nbits=3)
#' z <- CCNOT(c(1,2,3)) * (H(1) * x)
#'
#' @name ccnotgate
#'
#' @rdname ccnotgate
#' @aliases ccnotgate-class
#' @exportClass ccnotgate
setClass("ccnotgate",
         representation(bits="integer"),
         prototype(bits=c(1L, 2L, 3L)))

#' @export
ccnotgate <- function(bits=c(1, 2)) return(methods::new("ccnotgate", bits=as.integer(bits)))

#' The CCNOT or toffoli gate
#'
#' @param bits integer vector of length two, the first bit being the control and the second
#' the target bit.
#'
#' @aliases toffoli
#' @return
#' An S4 class 'ccnotgate' object is returned
#' @export
CCNOT <- function(bits=c(1, 2, 3)) return(methods::new("ccnotgate", bits=as.integer(bits)))

#' @rdname CCNOT
#' @aliases CCNOT
#' @export
toffoli <- function(bits=c(1, 2, 3)) return(methods::new("ccnotgate", bits=as.integer(bits)))

#' times-ccnotgate-qstate
#'
#' Applies a CCNOT (or toffoli) gate to a quantum state.
#'
#' @param e1 object of S4 class 'ccnotgate'
#' @param e2 object of S4 class 'qstate'
#'
#' @return
#' An object of S4 class 'qstate'
setMethod("*", c("ccnotgate", "qstate"),
          function(e1, e2) {
            stopifnot(length(e1@bits) == 3)
            stopifnot(length(e1@bits) == length(unique(e1@bits)))
            stopifnot(all(e1@bits > 0) && all(e1@bits <= e2@nbits))
            ## control bit == 1
            al <- 0:(2^e2@nbits-1)
            cb <- is.bitset(al, bit=e1@bits[1]) & is.bitset(al, bit=e1@bits[2])
            ## target bit
            tb <- is.bitset(al, bit=e1@bits[3])
            x <- which(cb & tb)
            y <- which(cb & !tb)
            e2@coefs[c(x,y)]  <- e2@coefs[c(y,x)]
            ## again the circuit needs extension for plotting
            ngates <- length(e2@circuit$gatelist)
            e2@circuit$gatelist[[ngates+1]] <- list(type="CCNOT", bits=c(e1@bits))

            return(e2)
          }
          )
