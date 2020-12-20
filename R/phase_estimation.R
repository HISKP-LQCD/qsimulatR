#' phase_estimation
#'
#' phase estimation algorithm
#'
#' @param bitmask integer. Vector of qubits for the t qubit wide
#' register needed for the phase estimation
#' @param FUN a function implementing the controlled application of
#' a unitary operator U to the power 2^(j-1) to the state x. It's first
#' argument must be the control qubit 'c', the second the integer 'j'
#' and the third the state 'x'. Additional parameters can be passed
#' via '...'.
#' @param x   a 'qstate' object
#' @param ... additional parameter to be passed on to 'FUN'
#'
#' @examples
#' ## NOT^k = Id if k even
#' cnotwrapper <- function(c, j, x, t) {
#'   if(j == 1) return(CNOT(c(c, t)) * x)
#'   return(Id(t) * x)
#' }
#' x <- X(1) * qstate(3)
#' ## X has eigenvalues lambda=1 and lambda=-1
#' ## thus phases 0 and 1/2
#' x <- phase_estimation(bitmas=c(2:3), FUN=cnotwrapper, x=x, t=1)
#' x
#' 
#' @export
phase_estimation <- function(bitmask, FUN, x, ...) {
  ## Hadarmard gates for all t phase estimation registers
  for(i in bitmask) {
    x <- H(i) * x
  }
  j <- 1
  for(i in bitmask) {
    x <- FUN(c=i, j=j, x=x, ...)
    j <- j+1
  }
  return(qft(x, bits=bitmask, inverse=TRUE))
}
