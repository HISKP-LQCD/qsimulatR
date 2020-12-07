annotate_bitnames <- function(i, y, cbit=FALSE, qubitnames=NULL, cbitnames=NULL) {
  if(!cbit && !is.null(qubitnames) && (i <= length(qubitnames))) {
    text(x=0.25, y=y, labels=qubitnames[i], pos=2)
  }
  if(cbit && !is.null(cbitnames) && (i <= length(cbitnames))) {
    text(x=0.25, y=y, labels=cbitnames[i], pos=2)
  }
}

#' plot-qstate
#'
#' @description
#' Plots a circuit corresponding to a qstate object
#'
#' @aliases plot
#' 
#' @param x qstate object
#' @param y not used here
#' @param ... additional parameters to be passed on
#'
#' @importFrom graphics plot lines points arrows legend text
#' @return nothing is returned, but a plot created
#'
#' @include state.R
#' 
#' @examples
#' x <- qstate(2)
#' y <- H(1) * x
#' z <- CNOT(c(1,2)) * y
#' plot(z)
#' 
#' @exportMethod plot
setMethod("plot", signature(x = "qstate", y = "missing"),
          function(x, y, ...) {
            nbits <- x@nbits
            ncbits <- x@circuit$ncbits
            n <- nbits + ncbits
            ngates <- length(x@circuit$gatelist)
            ## compute xlim first
            ipos <- rep(1, times=n)
            if(ngates > 0) {
              gatelist <- x@circuit$gatelist
              for(i in c(1:ngates)) {
                if(is.na(gatelist[[i]]$bits[2])) {
                  ipos[gatelist[[i]]$bits[1]] <- ipos[gatelist[[i]]$bits[1]] + 1
                }
                else {
                  ipos[1:n] <- max(ipos) + 1
                }
              }
            }
            xmax <- max(ipos)
            ## prepare empty plot
            plot(NA, ann=FALSE, xlim=c(0,xmax), ylim=c(0,n+1), axes=FALSE, frame.plot=FALSE)
            ## plot qubit lines
            for(i in c(n:(ncbits+1))) {
              lines(x=c(0.3, xmax), y=c(i, i))
              annotate_bitnames(i=n-i+1, y=i, ...)
            }
            ## plot classical bit lines
            if(ncbits > 0) {
              for(i in c(ncbits:1)) {
                lines(x=c(0.3, xmax), y=c(i-0.025, i-0.025))
                lines(x=c(0.3, xmax), y=c(i+0.025, i+0.025))
                annotate_bitnames(i=ncbits-i+1, y=i, cbit=TRUE, ...)
              }
            }
            if(ngates > 0) {
              ipos <- rep(1, times=n)
              ## plot gates
              gatelist <- x@circuit$gatelist
              for(i in c(1:ngates)) {
                ## single qubit gates
                if(is.na(gatelist[[i]]$bits[2])) {
                  type <- gatelist[[i]]$type
                  if(gatelist[[i]]$type == "Rx" || gatelist[[i]]$type == "Ry" || gatelist[[i]]$type == "Rz") {
                    type <- paste0(gatelist[[i]]$type, "(", format(gatelist[[i]]$angle, digits=3), ")") 
                  }
                  legend(x=ipos[gatelist[[i]]$bits[1]], y=n+1-gatelist[[i]]$bits[1],
                         type, xjust=0.5, yjust=0.5,
                         x.intersp=-0.5, y.intersp=0.1,
                         bg="white")
                  ipos[gatelist[[i]]$bits[1]] <- ipos[gatelist[[i]]$bits[1]] + 1
                }
                ## multi qubit gates
                else {
                  xp <- max(ipos)
                  ipos[1:n] <- xp + 1 
                  if(!is.null(gatelist[[i]]$controlled)) {
                    if(gatelist[[i]]$controlled) {
                      type <- gatelist[[i]]$type
                      if(gatelist[[i]]$type == "Rx" || gatelist[[i]]$type == "Ry" || gatelist[[i]]$type == "Rz") {
                        type <- paste0(gatelist[[i]]$type, "(", format(gatelist[[i]]$angle, digits=3), ")") 
                      }
                      points(x=xp, y=n+1-gatelist[[i]]$bits[1], pch=19, cex=1.5)
                      lines(x=c(xp, xp), y=n+1-c(gatelist[[i]]$bits[1], gatelist[[i]]$bits[2]))
                      legend(x=xp, y=n+1-gatelist[[i]]$bits[2],
                             type, xjust=0.5, yjust=0.5,
                             x.intersp=-0.5, y.intersp=0.1,
                             bg="white")
                    }
                  }
                  if(gatelist[[i]]$type == "CNOT") {
                    points(x=xp, y=n+1-gatelist[[i]]$bits[1], pch=19, cex=1.5)
                    points(x=xp, y=n+1-gatelist[[i]]$bits[2], pch=10, cex=2.5)
                    lines(x=c(xp, xp), y=n+1-c(gatelist[[i]]$bits[1], gatelist[[i]]$bits[2]))
                  }
                  if(gatelist[[i]]$type == "CCNOT") {
                    points(x=xp, y=n+1-gatelist[[i]]$bits[1], pch=19, cex=1.5)
                    points(x=xp, y=n+1-gatelist[[i]]$bits[2], pch=19, cex=1.5)
                    points(x=xp, y=n+1-gatelist[[i]]$bits[3], pch=10, cex=2.5)
                    lines(x=c(xp,xp), y=n+1-range(gatelist[[i]]$bits))
                  }
                  if(gatelist[[i]]$type == "SWAP") {
                    points(x=xp, y=n+1-gatelist[[i]]$bits[1], pch=4, cex=2.5)
                    points(x=xp, y=n+1-gatelist[[i]]$bits[2], pch=4, cex=2.5)
                    lines(x=c(xp, xp), y=n+1-c(gatelist[[i]]$bits[1], gatelist[[i]]$bits[2]))
                  }                
                  if(gatelist[[i]]$type == "CSWAP") {
                    points(x=xp, y=n+1-gatelist[[i]]$bits[1], pch=19, cex=1.5)
                    points(x=xp, y=n+1-gatelist[[i]]$bits[2], pch=4, cex=2.5)
                    points(x=xp, y=n+1-gatelist[[i]]$bits[3], pch=4, cex=2.5)
                    lines(x=c(xp, xp), y=n+1-range(gatelist[[i]]$bits))
                  }
                  if(gatelist[[i]]$type == "measure") {
                    lines(x=c(xp, xp), y=c(n+1-gatelist[[i]]$bits[1], ncbits+1-gatelist[[i]]$bits[2]))
                    points(x=xp, y=n+1-gatelist[[i]]$bits[1], pch=19, cex=1.5)
                    legend(x=xp, y=ncbits+1-gatelist[[i]]$bits[2],
                           "M", xjust=0.5, yjust=0.5,
                           x.intersp=-0.5, y.intersp=0.1,
                           bg="white")
                    arrows(x0=xp-0.2, x1=xp+0.2,
                           y0=ncbits+1-gatelist[[i]]$bits[2]-0.2,
                           y1=ncbits+1-gatelist[[i]]$bits[2]+0.2,
                           length=0.1)
                  }
                }
              }
            }
          }
          )
