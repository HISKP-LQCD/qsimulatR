is.bitset <- function(x, bit) {
  return((x %/% 2^(bit-1)) %% 2 != 0)
}

normalise <- function(x) {
  n = 1/sqrt(sum(Re(x*Conj(x))))
  return(as.complex(x*n))
}

check_qstate <- function(object) {
  stopifnot(object@nbits > 0)
  stopifnot(object@nbits <= 16)
  N <- 2^object@nbits
  stopifnot(N == length(object@coefs))
  stopifnot(1 == length(object@basis) || N == length(object@basis))
}

genStateString <- function(int, nbits, collapse="") {
  x <- intToBits(int)
  i <- nbits:1
  str <- paste0(ifelse(x[i] > 0, 1, 0), collapse=collapse)
  str <- paste0("|", str, ">")
  return(str)
}

#' genComputationalBasis
#'
#' function to generate the basis strings for given number of
#' bits
#'
#' @param nbits integer. The number of qbits
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
  basis <- c()
  N <- 2^nbits
  basis <- sapply(0:(N-1), genStateString, nbits=nbits, collapse=collapse)
  return(basis)
}

qstatecoefs <- function(y) {
  if (methods::is(y, 'qstate')) y@coefs
  else y
}

eps <- 1.e-12
#' The qstate class
#'
#' This class represents a quantum state
#'
#' @slot nbits The number of qbits
#' @slot coefs The 2^nbits complex valued vector of coefficients
#' @slot basis String or vector of strings. A single string will be interpreted 
#' as the \code{collapse}-parameter in \code{genComputationalBasis}. A vector 
#' of length 2^nbits yields the basis directly.
#'
#' @details
#' The qbits are counted from 1 to \code{nbits} starting with the least
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
#' @name qstate
#' @rdname qstate
#' @aliases qstate-class
#' @exportClass qstate
setClass("qstate",
         representation(nbits="integer",
                        coefs="complex",
                        basis="character",
                        circuit="list"),
         prototype(nbits=1L,
                   coefs=c(1. + 0i, 0i),
                   basis=genComputationalBasis(1L),
                   circuit=list(ncbits=0, gatelist=list())),
         validity=check_qstate)

## "constructor" function
#' @export
qstate <- function(nbits=1L,
                   coefs=c(1+0i, rep(0i, times=2^nbits-1)),
                   basis=genComputationalBasis(nbits=nbits),
                   circuit=list(ncbits=0, gatelist=list())) {
  return(methods::new("qstate", nbits=as.integer(nbits),
                      coefs=normalise(coefs),
                      basis=as.character(basis),
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

setMethod("plot", signature(x = "qstate"),
          function(x, ...) {
            nbits <- x@nbits
            ncbits <- x@circuit$ncbits
            n <- nbits + ncbits
            ngates <- length(x@circuit$gatelist)
            plot(NA, ann=FALSE, xlim=c(0,ngates+1), ylim=c(0,n+1), axes=FALSE, frame.plot=FALSE)
            for(i in c(n:(ncbits+1))) {
              lines(x=c(0, ngates+1), y=c(i, i))
            }
            if(ncbits > 0) {
              for(i in c(ncbits:1)) {
                lines(x=c(0, ngates+1), y=c(i-0.05, i-0.05))
                lines(x=c(0, ngates+1), y=c(i+0.05, i+0.05))
              }
            }
            gatelist <- x@circuit$gatelist
            for(i in c(1:ngates)) {
              if(is.na(gatelist[[i]]$bit2)) {
                legend(x=i, y=gatelist[[i]]$bit1,
                       gatelist[[i]]$type, xjust=0.5, yjust=0.5,
                       x.intersp=-0.5, y.intersp=0.1,
                       bg="white")
              }
              else {
                if(gatelist[[i]]$type == "CNOT") {
                  points(x=i, y=gatelist[[i]]$bit1, pch=19, cex=1.5)
                  points(x=i, y=gatelist[[i]]$bit2, pch=10, cex=2.5)
                  lines(x=c(i,i), y=c(gatelist[[i]]$bit1, gatelist[[i]]$bit2))
                }
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
#' The qbits are counted from 1 to \code{nbits} starting with the least
#' significant bit.
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
                   type="1"),
         validity=check_sqgate)

#' @export
sqgate <- function(bit=1L, M=array(as.complex(c(1,0,0,1)), dim=c(2,2)), type="1") {
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
            ii <- c(0, 2^(bit-1))
            res <- c()
            al <- 0:(2^e2@nbits-1)
            ii <- is.bitset(al, bit=bit)
            kk <- !ii
            res[kk] <- e1@M[1,1]*e2@coefs[kk] + e1@M[1,2]*e2@coefs[ii]
            res[ii] <- e1@M[2,1]*e2@coefs[kk] + e1@M[2,2]*e2@coefs[ii]
            ## the gatelist needs to be extended for plotting
            circuit <- e2@circuit
            ngates <- length(circuit$gatelist)
            circuit$gatelist[[ngates+1]] <- list(type=e1@type, bit1=e1@bit, bit2=NA)
            
            return(qstate(nbits=nbits, coefs=as.complex(res), basis=e2@basis, circuit=circuit))
          }
          )

## original version of the "*" above, significantly slower:
#### outside
##for(k in c(0:(2^(nbits-bit)-1))) {
##  ## inside
##  for(i in c(0:(2^(bit-1)-1))) {
##    ## the 2x2 matrix
##    for(j in c(0,1)) {
##      ll <- j*2^(bit-1) + i + k*2^bit + 1
##      rr <- ii + i + k*2^(bit) + 1
##      res[ll] <-  sum(e1@M[j+1,]*e2@coefs[rr])
##    }
##  }
##}


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

#' The CNOT gate
#'
#' This class represents a generic CNOT gate
#'
#' @slot bits Integer vector of length 2. First bit is the control bit,
#'            second the target bit.
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
#' @aliases "*"
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
            tmp <- e2@coefs[x]
            e2@coefs[x] <- e2@coefs[y]
            e2@coefs[y] <- tmp
            ## again the circuit needs extension for plotting
            ngates <- length(e2@circuit$gatelist)
            e2@circuit$gatelist[[ngates+1]] <- list(type="CNOT", bit1=e1@bits[1], bit2=e1@bits[2])

            return(e2)
          }
          )

#' Method measure
#' @name measure
#' @rdname measure-methods
#'
#' @description
#' performs a masurement on a `qstate` object.
#' 
#' @param e1 object to measure
#' @param bit bit to project on
#' @exportMethod measure
setGeneric("measure", function(e1, bit) attributes(e1))

#' @rdname measure-methods
#' @aliases measure
#'
#' @details \code{measure(e1)} performs a projection of the total wave function (i.e. all qbits).
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
#' and an element `value` with the value of the qbit `bit`.
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
            circuit <- e1@circuit
            return(list(psi=qstate(nbits=e1@nbits, coefs=as.complex(coefs), basis=e1@basis, circuit=circuit), value=value))
          }
          )

