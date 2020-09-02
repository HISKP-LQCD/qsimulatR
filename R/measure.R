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

truth.line <- function(i, gate, nbits) {
  eps <- 1e-14
  s <- rep(0, 2^nbits)
  s[i] <- 1

  x <- qstate(nbits, coefs=s)
  out <- gate * x
  k <- which(abs(out@coefs) > eps)
  out.bits <- rbind(sapply(k-1, genStateNumber, nbits=nbits))
  out.num <- apply(out.bits, 1, function(bits){
                     if(all(as.logical(bits))) return(1)
                     else if(!any(as.logical(bits))) return(0)
                     else return(NA)
                  })

  return(c(genStateNumber(i-1, nbits), out.num))
}

#' Method truth.table
#' @name truth.table
#' @rdname truth-table
#'
#' @details calculates the quantum truth table of the gate `e1`.
#' If a basis state is transformed to a superposition of basis states by the gate, the result is NA.
#'
#' @param e1 gate to measure.
#' @param nbits number of bits the gate acts on.
#' @param bits optional vector of length `nbits` containing the qubit order in the gate.
#'
#' @return
#' returns a data.frame containing the truth table. Each row corresponds
#' to one input-output combination. Each column to one specific bit.
#'
#' @include state.R
#'
#' @examples
#' ## truth table for a single bit gate
#' truth.table(X, 1)
#' ## for a 2-bit gate
#' truth.table(CNOT, 2)
#' ## for a 2-bit gate with reversed controll and target bits
#' truth.table(CNOT, bits=2:1)
#' ## for a general controlled gate
#' truth.table(cqgate, 2, gate=H(2))
#' 
#' @exportMethod truth.table
setGeneric("truth.table", function(e1, nbits, bits, ...) {
             stopifnot(!missing(nbits) || !missing(bits))
             if(missing(bits)) bits <- 1:nbits
             else nbits <- length(bits)

             gate <- e1(bits, ...)

             tab <- sapply(1:(2^nbits), truth.line, gate=gate, nbits=nbits)

             res <- as.data.frame(t(tab))
             names(res)[1:nbits]             <- paste0(rep("In", nbits),  (nbits-1):0)
             names(res)[(nbits+1):(2*nbits)] <- paste0(rep("Out", nbits), (nbits-1):0)
             return(res)
          }
)
