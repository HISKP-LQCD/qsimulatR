check_cqgate  <- function(object) {
  stopifnot(length(object@bits) == 2)
  stopifnot(all(object@bits > 0))
  stopifnot(length(object@bits) == length(unique(object@bits)))
}


#' A controlled single qubit gate
#'
#' This class represents a generic controlled gate
#'
#' @slot bits Integer. Integer vector of bits. The first is the control bit, the second the target bit.
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
#' x <- H(1) * qstate(nbits=2)
#' ## application of the CX (CNOT) gate to bit 1,2
#' z <- cqgate(bits=c(1L, 2L), gate=X(2L)) * x
#' z
#' ## the same as, but differently implemented
#' z <- CNOT(c(1,2)) * x
#' z
#' 
#' @name cqgate
#' @rdname cqgate
#' @aliases cqgate-class
#' @exportClass cqgate
setClass("cqgate",
         representation(bits="integer",
                        gate="sqgate"),
         prototype(bits=c(1L, 2L),
                   gate=Id(2L)),
         validity=check_cqgate
         )

#' @include state.R
#' @include sqgate.R
#' 
#' @export
cqgate <- function(bits=c(1L, 2L), gate=Id(2L)) {
  return(methods::new("cqgate", bits=as.integer(bits), gate=gate))
}

#' times-cqgate-qstate
#'
#' Applies a controlled single qubit gate to a quantum state.
#'
#' @param e1 object of S4 class 'cqgate'
#' @param e2 object of S4 class 'qstate'
#' @return
#' An object of S4 class 'qstate'
#'
setMethod("*", c("cqgate", "qstate"),
          function(e1, e2) {
            stopifnot(all(e1@bits > 0) && all(e1@bits <= e2@nbits))
            bits <- e1@bits
            nbits <- e2@nbits
            res <- e2@coefs
            al <- 0:(2^e2@nbits-1)
            cb  <- is.bitset(al, bit=bits[1])
            ii <- is.bitset(al, bit=bits[2]) & cb
            kk <- (!is.bitset(al, bit=bits[2])) & cb
            res[kk] <- e1@gate@M[1,1]*e2@coefs[kk] + e1@gate@M[1,2]*e2@coefs[ii]
            res[ii] <- e1@gate@M[2,1]*e2@coefs[kk] + e1@gate@M[2,2]*e2@coefs[ii]
            ## the gatelist needs to be extended for plotting
            circuit <- e2@circuit
            ngates <- length(circuit$gatelist)
            circuit$gatelist[[ngates+1]] <- list(type=e1@gate@type, bits=c(e1@bits, NA))
            if(e1@gate@type == "Rx" || e1@gate@type == "Ry" || e1@gate@type == "Rz" || grepl("^R[0-9]+", e1@gate@type))
              circuit$gatelist[[ngates+1]]$angle <- Arg(e1@gate@M[2,2])
            circuit$gatelist[[ngates+1]]$controlled <- TRUE
            
            return(qstate(nbits=nbits, coefs=as.complex(res), basis=e2@basis, circuit=circuit))
          }
          )
