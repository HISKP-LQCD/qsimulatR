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
#' @param repetitions number of measurements
#' @docType methods
#' @exportMethod measure
setGeneric("measure", function(e1, bit=NA, repetitions=NA) attributes(e1))

#' @rdname measure-methods
#' @aliases measure
#'
#' @details \code{measure(e1, bit, repetitions)} performs `repetitions` many 
#' projections/measurements of the qubit `bit`. If `bit` is not given 
#' explicitly, all qubits are projected.
#'
#' @return
#' \code{measure(e1, bit, repetitions)} returns a list with the measured `bit`, 
#' the number of `repetitions`, the probability distribution of all states 
#' `prob` and the results vector `value`. If all bits are measured, the 
#' basis is added to the list as `basis`. The collapsed state is stored as 
#' `psi` if exactly one measurement is performed.
#' In the case of a single qubit measurement `value` is of length `repetitions` 
#' and contains all the results of this projection. Otherwise `value` is of 
#' length 2^nbits and it contains the counts how often each state has been 
#' obtained.
#'
#' @importFrom stats runif
#' @examples
#' ## measure the separate bits
#' x <- H(1) * (H(2) * qstate(nbits=2))
#' summary(measure(x, bit=1))
#' hist(measure(x, rep=100))
setMethod("measure", c("qstate"),
          function(e1, bit=NA, repetitions=1) {
            stopifnot(is.na(bit) || bit %in% c(1:e1@nbits))
            prob <- Re(e1@coefs * Conj(e1@coefs))

            res <- list(bit=bit, repetitions=repetitions, prob=prob)
            if(is.na(bit)){
              res$nbits <- e1@nbits
              res$basis <- e1@basis
              value <- stats::rmultinom(n=1, size=repetitions, prob=prob)[,1]
              if(repetitions == 1) coefs <- value
            }else{
              N <- 2^e1@nbits
              ii <- which(is.bitset(0:(N-1), bit))
              is1 <- sum(prob[ii])
              value <- ifelse(runif(repetitions) < is1, 1, 0)
              if(repetitions == 1){
                coefs <- e1@coefs
                if(value == 1) coefs[-ii] <- 0i
                else coefs[ii] <- 0i
              }
            }
            res$value <- value

            if(repetitions == 1){
              ngates <- length(e1@circuit$gatelist)
              if(is.na(bit)){
                e1@circuit$gatelist[(ngates+1):(ngates+e1@nbits)] <- lapply(1:e1@nbits, function(bit) list(type="measure", bits=c(bit, e1@circuit$ncbits+bit, NA)))
                e1@circuit$ncbits <- e1@circuit$ncbits + e1@nbits
              }else{
                cbit <- e1@circuit$ncbits+1
                e1@circuit$gatelist[[ngates+1]] <- list(type="measure", bits=c(bit, cbit, NA))
                e1@circuit$ncbits <- cbit
              }
              res$psi <- qstate(nbits=e1@nbits, coefs=as.complex(coefs), basis=e1@basis, circuit=e1@circuit)
            }

            attr(res, "class") <- c("measurement", "list")
            return(invisible(res))
          }
          )

#' Summarize a quantum measurement
#'
#' @param object as returned by \code{measure}
#' @param ... Generic parameters to pass on, not used here.
#'
#' @importFrom methods show
#' @return
#' No return value.
#' 
#' @export
summary.measurement <- function(object, ...) {
  if(is.na(object$bit)){
    cat("All bits have been measured", object$repetitions, "times with the outcome:\n")
    tmp.state <- qstate(object$nbits, basis=object$basis)
    tmp.state@coefs <- as.complex(object$value)
    show(tmp.state)
  }else{
    cat("Bit", object$bit, "has been measured", object$repetitions, "times with the outcome:\n")
    ones <- sum(object$value)
    zeros <- object$repetitions - ones
    cat("0: ", zeros, "\n1: ", ones, "\n")
  }
}

#' Plot the histogram of a quantum measurement
#'
#' @param x object as returned by \code{measure}
#' @param only.nonzero are the states with zero measurements to be plotted?
#' @param by.name shall the xlabel contain the basis names? If `FALSE`, the 
#' index number is used.
#' @param freq shall the total counts be plotted? If not, the values are 
#' normalised to 1.
#' @param ... Generic parameters to pass on to \code{barplot()}
#'
#' @importFrom graphics barplot
#' 
#' @return
#' No return value.
#' 
#' @export
hist.measurement <- function(x, only.nonzero=TRUE, by.name=only.nonzero, freq=TRUE, ...) {
  if(is.na(x$bit)){
    if(only.nonzero) mask <- which(x$value > 0)
    else mask <- 1:(length(x$value))
    counts <- x$value[mask]
    if(by.name){
      if(length(x$basis) > 1) names.arg <- x$basis[mask]
      else{
        names.arg <- sapply(mask-1, genStateString, nbits=x$nbits, collapse=x$basis)
      }
    }else{
      names.arg <- mask
    }
  }else{
    names.arg <- 0:1
    ones <- sum(x$value)
    zeros <- x$repetitions - ones
    counts <- c(zeros, ones)
  }

  if(freq){
    ylab <- "Counts"
  }else{
    counts <- counts/x$repetitions
    ylab <- "Probability"
  }

  barplot(counts, ylab=ylab, names.arg=names.arg, ...)
}

truth.line <- function(i, gate, nbits) {
  eps <- 1e-14
  s <- rep(0, 2^nbits)
  s[i] <- 1
  x <- qstate(nbits, coefs=s)

  ## gate could be something more complicated than a default gate with overloaded '*'
  if(is.function(gate)){
    out <- gate(x)
  }else{
    out <- gate * x
  }
  k <- which(abs(out@coefs) > eps)
  out.bits <- rbind(sapply(k-1, genStateNumber, nbits=nbits))
  out.num <- apply(X=out.bits, MARGIN=1,
                   FUN=function(bits){
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
#' If a basis state is transformed to a superposition of basis states by
#' the gate, the result is 'NA'.
#'
#' @param e1 gate to measure.
#' @param nbits number of bits the gate acts on.
#' @param bits optional vector of length `nbits` containing the qubit order in the gate.
#' @param ... additional parameters to passed be on to 'e1'
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
#' ## for an arbitrary circuit (here a swap implementation using only CNOTs)
#' myswap <- function(bits){ function(x){ CNOT(bits) * (CNOT(rev(bits)) * (CNOT(bits) * x))}}
#' truth.table(myswap, 2)
#' 
#' @exportMethod truth.table
setGeneric("truth.table",
           function(e1, nbits, bits, ...) {
             stopifnot(!missing(nbits) || !missing(bits))
             if(missing(bits)) bits <- 1:nbits
             else nbits <- length(bits)
             
             gate <- e1(bits, ...)
             
             tab <- sapply(X=1:(2^nbits), FUN=truth.line, gate=gate, nbits=nbits)
             
             res <- as.data.frame(t(tab))
             names(res)[1:nbits]             <- paste0(rep("In", nbits),  nbits:1)
             names(res)[(nbits+1):(2*nbits)] <- paste0(rep("Out", nbits), nbits:1)
             return(res)
           }
           )
