check_cnqgate  <- function(object) {
  stopifnot(length(object@cbits) > 1)
  stopifnot(all(object@cbits > 0))
  stopifnot(length(object@cbits) == length(unique(object@cbits)))
  stopifnot(length(object@tbit) == 1)
  stopifnot(length(object@cbits) == length(object@inverse))
}


#' n-fold controlled single qubit gate
#'
#' This class represents a generic n-fold controlled gate
#'
#' @slot cbits Integer. Integer vector of control bits.
#' @slot tbit Integer. Target bit.
#' @slot gate sqgate. The single qubit gate.
#' @slot inverse Logical. Boolean vector of same length as \code{cbits}. If 
#' \code{TRUE}, the corresponding control bit is negated.
#'
#' @details
#' The qubits are counted from 1 to \code{nbits} starting with the least
#' significant bit.
#'
#' @include state.R
#' @include sqgate.R
#' 
#' @examples
#' x <- H(1) * qstate(nbits=3)
#' ## application of the CCX (CCNOT) gate to bits 1,2 and 3
#' z <- cnqgate(cbits=c(1L, 2L), tbit=3L, gate=X(3L)) * x
#' z
#' ## the same, but differently implemented
#' z <- CCNOT(c(1,2,3)) * x
#' z
#' 
#' @name cnqgate
#' @rdname cnqgate
#' @aliases cnqgate-class
#' @exportClass cnqgate
setClass("cnqgate",
         representation(cbits="integer",
                        tbit="integer",
                        gate="sqgate",
                        inverse="logical"),
         # prototype(cbits=c(1L, 2L),
         #                  tbit=3L,
         #           gate=Id(3L),
         #                  inverse=rep(FALSE, length(cbits))),
         validity=check_cnqgate
)

#' @include state.R
#' @include sqgate.R
#' 
#' @export
cnqgate <- function(cbits=c(1L, 2L), tbit=3L, gate=Id(3L), inverse=rep(FALSE, length(cbits))) {
  return(methods::new("cnqgate", cbits=as.integer(cbits), tbit=as.integer(tbit), gate=gate, inverse=as.logical(inverse)))
}

#' times-cnqgate-qstate
#'
#' Applies n-fold controlled single qubit gate to a quantum state.
#'
#' @param e1 object of S4 class 'cnqgate'
#' @param e2 object of S4 class 'qstate'
#' @return
#' An object of S4 class 'qstate'
#'
setMethod("*", c("cnqgate", "qstate"),
          function(e1, e2) {
            stopifnot(all(e1@cbits <= e2@nbits))
            bits <- e1@cbits
            target <- e1@tbit
            nbits <- e2@nbits
            res <- e2@coefs
            al <- 0:(2^nbits-1)
            cb <- apply(sapply(bits, is.bitset, x=al), 1, function(set) all(xor(set, e1@inverse)))
            ii <- is.bitset(al, bit=target) & cb
            kk <- (!is.bitset(al, bit=target)) & cb
            res[kk] <- e1@gate@M[1,1]*e2@coefs[kk] + e1@gate@M[1,2]*e2@coefs[ii]
            res[ii] <- e1@gate@M[2,1]*e2@coefs[kk] + e1@gate@M[2,2]*e2@coefs[ii]
            ## the gatelist needs to be extended for plotting
            circuit <- e2@circuit
            ngates <- length(circuit$gatelist)
            circuit$gatelist[[ngates+1]] <- list(type=e1@gate@type, bits=c(bits, target, NA), inverse=e1@inverse)
            if(e1@gate@type == "Rx" || e1@gate@type == "Ry" || e1@gate@type == "Rz") circuit$gatelist[[ngates+1]]$angle <- -Re(2*1i*log(e1@M[2,2]))
            circuit$gatelist[[ngates+1]]$controlled <- TRUE

            result <- qstate(nbits=nbits, coefs=as.complex(res), basis=e2@basis, noise=e2@noise, circuit=circuit)
            if(e1@gate@type == "ERR" || ! any(bits %in% e2@noise$bits) || e2@noise$p < runif(1)){
              return(result)
            }else{
              return(noise(bits[bits %in% e2@noise$bits], error=e2@noise$error, args=e2@noise$args) * result)
            }
          }
)
