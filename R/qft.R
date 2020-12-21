#' qft
#'
#' Quantum Fourier Trafo
#'
#' @description
#' performs the quantum Fourier Trafo on the qstate x and the
#'   specified list of qubits.
#'
#' @param x qstate
#' @param inverse boolean. If 'TRUE', perform inverse transform
#' @param bits integer. list of qubits to include in the trafo. if
#'   missing, `bits=c(1:n)` is assumed, with `n` the number of qubits
#'   in `x`.
#'
#' @include cqgate.R
#' 
#' @return
#' a `qstate` object with the quantum Fourier trafo of input `x`.
#'
#' @details
#' The Fourier Trafo is defined as
# \eqn{|j\rangle \to 1/\sqrt{N}\sum_{k=0}^{N-1} e^{2\pi i j k/N}|k\rangle}
#' \deqn{|j> -> 1/sqrt(N) sum_k=0^N_1 exp(2 pi i j k/N) |k>}
#' the inverse with the oposite sign in the exponential.
#'
#' @examples
#' x <- qstate(3)
#' y <- qft(x)
#' z <- qft(y, inverse=TRUE)
#' 
#' @export
qft <- function(x, inverse=FALSE, bits) {
  check_qstate(x)
  n <- x@nbits
  if(missing(bits)) bits <- c(1:n)
  else stopifnot((length(bits) <= n) &&
                 all(bits <= n) && all(bits > 0) &&
                 (length(bits) == length(unique(bits))))
  
  y <- x
  sign <- +1
  if(is.numeric(inverse)) sign <- inverse
  else if(inverse) sign <- -1
  for(bit in rev(seq_along(bits))) {
    y <- H(bits[bit]) * y
    if(bit > 1) {
      for(i in c((bit-1):1)) {
        y <- cqgate(bits=c(bits[i], bits[bit]), gate=Ri(bits[i], bit-(i-1), sign=sign)) * y
      }
    }
  }
  ## reverse order
  n <- length(bits)
  if(n > 1) {
    for(k in c(1:floor(n/2))) {
      y <- SWAP(c(bits[k], bits[n-(k-1)])) * y
    }
  }
  return(invisible(y))
}


#' cqft
#'
#' Controlled Quantum Fourier Trafo
#'
#' @description
#' performs the controlled quantum Fourier Trafo on the qstate x and the
#'   specified list of qubits.
#'
#' @param c integer. a single control qubit.
#' @param x qstate. state the qft will applied to
#' @param inverse boolean. If 'TRUE', perform inverse transform
#' @param bits integer. list of qubits to include in the trafo. if
#'   missing, `bits=c(1:n)[-c]` is assumed, with `n` the number of qubits
#'   in `x`.
#'
#' @include cqgate.R
#' 
#' @return
#' a `qstate` object with the quantum Fourier trafo of input `x`.
#'
#' @details
#' The Fourier Trafo is defined as
# \eqn{|j\rangle \to 1/\sqrt{N}\sum_{k=0}^{N-1} e^{2\pi i j k/N}|k\rangle}
#' \deqn{|j> -> 1/sqrt(N) sum_k=0^N_1 exp(2 pi i j k/N) |k>}
#' the inverse with the oposite sign in the exponential.
#'
#' @examples
#' x <- qstate(3)
#' y <- cqft(1, x)
#' z <- cqft(1, y, inverse=TRUE)
#' 
#' @export
cqft <- function(c, x, inverse=FALSE, bits) {
  if(missing(c)) stop("control bit must be specified!")
  stopifnot(length(c) == 1)
  check_qstate(x)
  n <- x@nbits
  if(missing(bits)) bits <- c(1:n)[-c]
  else stopifnot((length(bits) <= n) &&
                 all(bits <= n) && all(bits > 0) &&
                 (length(bits) == length(unique(bits))))
  stopifnot(all(bits != c))

  y <- x
  sign <- +1
  if(is.numeric(inverse)) sign <- inverse
  else if(inverse) sign <- -1
  for(bit in rev(seq_along(bits))) {
    y <- cqgate(bits=c(c, bits[bit]), gate=H(bits[bit])) * y
    if(bit > 1) {
      for(i in c((bit-1):1)) {
        y <- ccqgate(bits=c(c, bits[i], bits[bit]), gate=Ri(bits[i], bit-(i-1), sign=sign)) * y
      }
    }
  }
  ## reverse order
  n <- length(bits)
  for(k in c(1:floor(n/2))) {
    y <- CSWAP(c(c, bits[k], bits[n-(k-1)])) * y
  }
  return(invisible(y))
}
