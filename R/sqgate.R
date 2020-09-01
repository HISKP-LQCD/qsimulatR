check_sqgate  <- function(object) {
  stopifnot(length(object@bit) == 1)
  stopifnot(object@bit > 0)
  stopifnot(all(dim(object@M) == c(2,2)))
}

#' A single qubit gate
#'
#' This class represents a generic single qubit gate
#'
#' @slot bit Integer. The single bit to act on.
#' @slot M complex valued array. The 2x2 matrix representing the
#' gate
#'
#' @details
#' The qubits are counted from 1 to \code{nbits} starting with the least
#' significant bit.
#'
#' @include state.R
#' 
#' @examples
#' x <- qstate(nbits=2)
#' ## application of the X (NOT) gate to bit 1
#' z <- sqgate(bit=1L, M=array(as.complex(c(0,1,1,0)), dim=c(2,2))) * x
#' z
#' 
#' @name sqgate
#' @rdname sqgate
#' @aliases sqgate-class
#' @exportClass sqgate
setClass("sqgate",
         representation(bit="integer",
                        M="array",
                        type="character"),
         prototype(bit=c(1L),
                   M=array(as.complex(c(1,0,0,1)), dim=c(2,2)),
                   type="Id"),
         validity=check_sqgate)

#' @export
sqgate <- function(bit=1L, M=array(as.complex(c(1,0,0,1)), dim=c(2,2)), type="Id") {
  return(methods::new("sqgate", bit=as.integer(bit), M=M, type=type))
}

#' times-sqgate-qstate
#'
#' Applies a single qubit gate to a quantum state.
#'
#' @param e1 object of S4 class 'sqgate'
#' @param e2 object of S4 class 'qstate'
#' @return
#' An object of S4 class 'qstate'
#'
setMethod("*", c("sqgate", "qstate"),
          function(e1, e2) {
            stopifnot(e1@bit > 0 && e1@bit <= e2@nbits)
            bit <- e1@bit
            nbits <- e2@nbits
            res <- c()
            al <- 0:(2^e2@nbits-1)
            ii <- is.bitset(al, bit=bit)
            kk <- !ii
            res[kk] <- e1@M[1,1]*e2@coefs[kk] + e1@M[1,2]*e2@coefs[ii]
            res[ii] <- e1@M[2,1]*e2@coefs[kk] + e1@M[2,2]*e2@coefs[ii]
            ## the gatelist needs to be extended for plotting
            circuit <- e2@circuit
            ngates <- length(circuit$gatelist)
            circuit$gatelist[[ngates+1]] <- list(type=e1@type, bits=c(e1@bit, NA, NA))
            if(e1@type == "Rz") circuit$gatelist[[ngates+1]]$angle <- -Re(2*1i*log(e1@M[2,2]))
            
            return(qstate(nbits=nbits, coefs=as.complex(res), basis=e2@basis, circuit=circuit))
          }
          )

#' The Hadarmard gate
#'
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- qstate(nbits=2)
#' z <- H(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
H <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1,1,1,-1)), dim=c(2,2))/sqrt(2), type="H"))
}
#' The identity gate
#'
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- qstate(nbits=2)
#' z <- Id(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
Id <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1,0,0,1)), dim=c(2,2)), type="Id"))
}
#' The Rz gate
#' 
#' @param bit integer. The bit to which to apply the gate
#' @param theta numeric. angle
#' 
#' @examples
#' x <- qstate(nbits=2)
#' z <- Rz(1, pi/4) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
Rz <- function(bit, theta=0.) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(exp(-1i*theta/2), 0, 0, exp(1i*theta/2))), dim=c(2,2)), type="Rz"))
}
#' The S gate
#' 
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- X(1) * qstate(nbits=2)
#' z <- S(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
S <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1,0,0,1i)), dim=c(2,2)), type="S"))
}
#' The Tgate gate
#' 
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- X(1)*qstate(nbits=2)
#' z <- Tgate(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
Tgate <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1., 0, 0, exp(1i*pi/4))), dim=c(2,2)), type="T"))
}
#' The X gate
#' 
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- qstate(nbits=2)
#' z <- X(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
X <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(0., 1., 1., 0.)), dim=c(2,2)), type="X"))
}
#' The Y gate
#' 
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- qstate(nbits=2)
#' z <- Y(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
Y <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(0., -1i, 1i, 0.)), dim=c(2,2)), type="Y"))
}
#' The Z gate
#' 
#' @param bit integer. The bit to which to apply the gate
#'
#' @examples
#' x <- X(1) * qstate(nbits=2)
#' z <- Z(1) * x
#' z
#' 
#' @return
#' An S4 class 'sqgate' object is returned
#' @export
Z <- function(bit) {
  return(methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1., 0., 0., -1.)), dim=c(2,2)), type="Z"))
}
