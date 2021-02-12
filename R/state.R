#' is.bitset
#'
#' checks whether or not a bit is set in target
#'
#' @param x target vector
#' @param bit integer. The bit to check
#'
#' @return a boolean vector
is.bitset <- function(x, bit) {
  return((x %/% 2^(bit-1)) %% 2 != 0)
}

#' normalise
#'
#' Normalises a complex vector to 1
#'
#' @param x complex valued vector
#'
#' @return
#' Returns the normalised complex valued vector
normalise <- function(x) {
  n = 1/sqrt(sum(Re(x*Conj(x))))
  return(as.complex(x*n))
}

check_qstate <- function(object) {
  stopifnot(object@nbits > 0)
  stopifnot(object@nbits <= 24)
  N <- 2^object@nbits
  stopifnot(N == length(object@coefs))
  stopifnot(1 == length(object@basis) || N == length(object@basis))
}

#' genStateNumber
#'
#' function to generate the bit representation for a specific basis state
#'
#' @param int integer number representing the basis state
#' @param nbits integer. The number of qubits
#'
#' @return a integer vector of length `nbits`
#' 
#' @examples
#' genStateNumber(5, 4)
#' genStateNumber(2, 2)
#'
#' @export
genStateNumber <- function(int, nbits) {
  x <- intToBits(int)
  i <- nbits:1
  return(ifelse(x[i] > 0, 1, 0))
}

#' genStateString
#'
#' function to generate the string for a specific basis state
#'
#' @param int integer number representing the basis state
#' @param nbits integer. The number of qubits
#' @param collapse character. String to fill in between separate bits
#'
#' @return a character
#' 
#' @examples
#' genStateString(5, 4)
#' genStateString(2, 2, collapse=">|")
#'
#' @export
genStateString <- function(int, nbits, collapse="") {
  str <- paste0(genStateNumber(int, nbits), collapse=collapse)
  str <- paste0("|", str, ">")
  return(str)
}

#' genComputationalBasis
#'
#' function to generate the basis strings for given number of
#' bits
#'
#' @param nbits integer. The number of qubits
#' @param collapse character. String to fill in between separate bits
#'
#' @return a character vector of length 2^nbits
#' 
#' @examples
#' genComputationalBasis(4)
#' genComputationalBasis(2, collapse=">|")
#'
#' @export
genComputationalBasis <- function(nbits, collapse="") {
  N <- 2^nbits
  basis <- sapply(0:(N-1), genStateString, nbits=nbits, collapse=collapse)
  return(basis)
}

#' genNoise
#'
#' function to generate the noise list
#'
#' See function \code{noise} for details.
#'
#' @param nbits integer. The number of qubits
#' @param p probability with which noise is applied after every gate
#' @param bits integer or integer array. The bit to which to apply the gate.
#' @param error String containing the error model.
#' @param ... Additional arguments to be stored in \code{args}.
#'
#' @return a list containing \code{p}, \code{bits}, \code{error} and 
#' \code{args}
#' 
#' @examples
#' genNoise(4)
#' genNoise(2, p=1, error="small", sigma=0.1)
#'
#' @export
genNoise <- function(nbits, p=0, bits=1:nbits, error="any", ...) {
  if(missing(nbits) && missing(bits)){
    stop("nbits or bits have to be provided!")
  }
  args <- list(...)
  noise <- list(p=p, bits=bits, error=error, args=args)
  return(noise)
}

qstatecoefs <- function(y) {
  if (methods::is(y, 'qstate')) y@coefs
  else y
}

#' The qstate class
#'
#' This class represents a quantum state
#'
#' @slot nbits The number of qubits
#' @slot coefs The 2^nbits complex valued vector of coefficients
#' @slot basis String or vector of strings. A single string will be interpreted 
#' as the \code{collapse}-parameter in \code{genComputationalBasis}. A vector 
#' of length 2^nbits yields the basis directly.
#' @slot noise List containing the probability \code{p} some noise is applied 
#' to one of the \code{bits} after a gate application, the model 
#' \code{error} of this noise and further arguments \code{args} to be passed to the 
#' function \code{noise}. See function \code{noise} for details.
#' The list \code{noise} can be generated with \code{genNoise}.
#' @slot circuit List containing the number of non-quantum bits \code{ncbits}
#' and a list of gates \code{gatelist} applied to the original state.
#' Filled automatically as gates are applied, required for plotting.
#'
#' @details
#' The qubits are counted from 1 to \code{nbits} starting with the least
#' significant bit.
#' 
#' @examples
#' x <- qstate(nbits=2)
#' x
#'
#' x <- qstate(nbits=2, coefs=as.complex(sqrt(rep(0.25, 4))), basis=",")
#' x
#'
#' x <- qstate(nbits=1, coefs=as.complex(sqrt(rep(0.5, 2))), basis=c("|dead>", "|alive>"))
#' x
#'
#' x <- qstate(nbits=2, noise=genNoise(nbits=2, p=1))
#' Id(2) * x
#' 
#' x <- qstate(nbits=3, noise=genNoise(p=1, bits=1:2, error="small", sigma=0.1))
#' Id(2) * x
#' 
#' @name qstate
#' @rdname qstate
#' @aliases qstate-class
#' @exportClass qstate
setClass("qstate",
         representation(nbits="integer",
                        coefs="complex",
                        basis="character",
                        noise="list",
                        circuit="list"),
         prototype(nbits=1L,
                   coefs=c(1. + 0i, 0i),
                   basis=genComputationalBasis(1L),
                   noise=genNoise(1L),
                   circuit=list(ncbits=0, gatelist=list())),
         validity=check_qstate)

## "constructor" function
#' @export
qstate <- function(nbits=1L,
                   coefs=c(1+0i, rep(0i, times=2^nbits-1)),
                   basis=genComputationalBasis(nbits=nbits),
                   noise=genNoise(nbits=nbits),
                   circuit=list(ncbits=0, gatelist=list())) {
  return(methods::new("qstate", nbits=as.integer(nbits),
                      coefs=normalise(coefs),
                      basis=as.character(basis),
                      noise=noise,
                      circuit=circuit))
}

## Convert to numeric vector if and only if all imaginary parts equal zero
SquashIm <- function(x) {
  if (all(Im(z <- zapsmall(x))==0)) as.numeric(z)
  else x
}

setMethod("show", signature(object = "qstate"),
          function(object) {
            N <- 2^object@nbits
            first <- TRUE
            coefs <- SquashIm(object@coefs)
            eps <- 1.e-12
            for(i in 1:N) {
              if(abs(object@coefs[i]) > eps) {
                if((i != 1) && !first) cat(" + ")
                else cat("   ")
                if(length(object@basis) == 1){
                  cat("(", coefs[i], ")\t*", genStateString(i-1, object@nbits, collapse=object@basis), "\n")
                }else{
                  cat("(", coefs[i], ")\t*", object@basis[i], "\n")
                }
                first <- FALSE
              }
            }
          }
          )




#' times-matrix-qstate
#'
#' Applies a single qubit gate to a quantum state.
#'
#' @param e1 object of S4 class 'matrix'
#' @param e2 object of S4 class 'qstate'
#' @return
#' An object of S4 class 'qstate'
#'
setMethod("*", c("matrix", "qstate"),
          function(e1, e2) {
            e2@coefs = drop(e1 %*% e2@coefs)
            methods::validObject(e2)
            return(e2)
          }
          )

#' times-number-qstate
#'
#' Multiplies a quantum gate by a global (phase) factor.
#'
#' @param e1 object of S4 class 'complex'
#' @param e2 object of S4 class 'qstate'
#' @return
#' An object of S4 class 'qstate'
#'
setMethod("*", c("complex", "qstate"),
          function(e1, e2) {
            stopifnot(abs(abs(e1)-1) < 4*.Machine$double.eps)
            e2@coefs = e1 * e2@coefs
            methods::validObject(e2)
            return(e2)
          }
          )
