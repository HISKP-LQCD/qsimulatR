#' The qsimulatR Package
#'
#' A simulator for a quantum computer
#'
#' A quantum computer simulator framework. General single qubit gates
#' and general controlled single qubit gates can be easily defined.
#' For convenience, it currently directly provides
#' most common gates (X, Y, Z, H, Z, S, T, Rx, Ry, Rz, CNOT, SWAP,
#' toffoli or CCNOT, CSWAP). 'qsimulatR' supports plotting
#' of circuits and is able to export circuits into IBM's 'Qiskit' python
#' package, which can be run on IBM's real quantum hardware.
#' 'qsimulatR' currently works for up to 24 qubits (a virtual
#' restriction, which can be lifted).
#' 
#' @name qsimulatR
#' @docType package
#' @author Johann Ostemeyer, Carsten Urbach, \email{urbach@hiskp.uni-bonn.de}
#'
#' @keywords internal
"_PACKAGE"
NULL
