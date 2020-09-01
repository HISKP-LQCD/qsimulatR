#' Method measure
#' @name measure
#' @rdname measure-methods
#'
#' @description
#' performs a masurement on a `qstate` object.
#'
#' @include state.R
#' 
#' @param e1 object to measure
#' @param bit bit to project on
#' @docType methods
#' @exportMethod measure
setGeneric("measure", function(e1, bit) attributes(e1))

#' @rdname measure-methods
#' @aliases measure
#'
#' @details \code{measure(e1)} performs a projection of the total wave function (i.e. all qubits).
#'
#' @return
#' \code{measure(e1)} returns a `qstate` object representing the state projected on.
#'
#' @examples
#' ## project the total wave function
#' x <- H(1) * (H(2) * qstate(nbits=2))
#' measure(x)
setMethod("measure", c("qstate"),
          function(e1) {
            prob <- Re(e1@coefs * Conj(e1@coefs))
            e1@coefs <- as.complex(stats::rmultinom(n=1, size=1, prob=prob))
            return(e1)
          }
          )

#' @rdname measure-methods
#' @aliases measure
#'
#' @details \code{measure(e1, bit)} performs a projection/measurement of the abit `bit`.
#'
#' @return
#' \code{measure(e1, bit)} returns a list with an element `psi` the `qstate` object projected onto
#' and an element `value` with the value of the qubit `bit`.
#'
#' @importFrom stats runif
#' @examples
#' ## measure the separate bits
#' x <- H(1) * (H(2) * qstate(nbits=2))
#' measure(x, 1)
#' measure(x, 2)
setMethod("measure", c("qstate", "numeric"),
          function(e1, bit) {
            stopifnot(bit %in% c(1:e1@nbits))
            prob <- Re(e1@coefs * Conj(e1@coefs))
            N <- 2^e1@nbits
            ii <- which(!is.bitset(0:(N-1), bit))
            is0 <- sum(prob[ii])
            value <- 0
            coefs <- e1@coefs
            if(runif(1) < is0) coefs[-ii] <- 0i
            else {
              coefs[ii] <- 0i
              value <- 1
            }
            ngates <- length(e1@circuit$gatelist)
            cbit <- e1@circuit$ncbits+1
            e1@circuit$gatelist[[ngates+1]] <- list(type="measure", bits=c(bit, cbit, NA))
            e1@circuit$ncbits <- cbit
            return(list(psi=qstate(nbits=e1@nbits, coefs=as.complex(coefs), basis=e1@basis, circuit=e1@circuit), value=value))
          }
          )
