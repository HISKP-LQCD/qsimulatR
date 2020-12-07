check_ccqgate  <- function(object) {
  stopifnot(length(object@bits) == 3)
  stopifnot(all(object@bits > 0))
  stopifnot(length(object@bits) == length(unique(object@bits)))
}


#' A twice controlled single qubit gate
#'
#' This class represents a generic controlled gate
#'
#' @slot bits Integer. Integer vector of bits. The first two are the control bits, the third the target bit.
#' @slot gate sqgate. The single qubit gate.
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
#' ## application of the CCX (CCNOT) gate to bit 1,2,3
#' z <- ccqgate(bits=c(1L, 2L, 3L), gate=X(3L)) * x
#' z
#' ## the same, but differently implemented
#' z <- CCNOT(c(1,2,3)) * x
#' z
#' 
#' @name ccqgate
#' @rdname ccqgate
#' @aliases ccqgate-class
#' @exportClass ccqgate
setClass("ccqgate",
         representation(bits="integer",
                        gate="sqgate"),
         prototype(bits=c(1L, 2L, 3L),
                   gate=Id(3L)),
         validity=check_ccqgate
         )

#' @include state.R
#' @include sqgate.R
#' 
#' @export
ccqgate <- function(bits=c(1L, 2L, 3L), gate=Id(3L)) {
  return(methods::new("ccqgate", bits=as.integer(bits), gate=gate))
}

#' times-ccqgate-qstate
#'
#' Applies a twice controlled single qubit gate to a quantum state.
#'
#' @param e1 object of S4 class 'ccqgate'
#' @param e2 object of S4 class 'qstate'
#' @return
#' An object of S4 class 'qstate'
#'
setMethod("*", c("ccqgate", "qstate"),
          function(e1, e2) {
            stopifnot(all(e1@bits > 0) && all(e1@bits <= e2@nbits))
            bits <- e1@bits
            nbits <- e2@nbits
            res <- e2@coefs
            al <- 0:(2^e2@nbits-1)
            cb <- is.bitset(al, bit=e1@bits[1]) & is.bitset(al, bit=e1@bits[2])
            ii <- is.bitset(al, bit=bits[3]) & cb
            kk <- (!is.bitset(al, bit=bits[3])) & cb
            res[kk] <- e1@gate@M[1,1]*e2@coefs[kk] + e1@gate@M[1,2]*e2@coefs[ii]
            res[ii] <- e1@gate@M[2,1]*e2@coefs[kk] + e1@gate@M[2,2]*e2@coefs[ii]
            ## the gatelist needs to be extended for plotting
            circuit <- e2@circuit
            ngates <- length(circuit$gatelist)
            circuit$gatelist[[ngates+1]] <- list(type=e1@gate@type, bits=c(e1@bits, NA))
            if(e1@gate@type == "Rx" || e1@gate@type == "Ry" || e1@gate@type == "Rz") circuit$gatelist[[ngates+1]]$angle <- -Re(2*1i*log(e1@M[2,2]))
            circuit$gatelist[[ngates+1]]$controlled <- TRUE
            
            return(qstate(nbits=nbits, coefs=as.complex(res), basis=e2@basis, circuit=circuit))
          }
          )
