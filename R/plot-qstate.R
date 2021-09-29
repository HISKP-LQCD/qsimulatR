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
#' @param gate_x_plot if the number of gates exeed this number multiple plots will
#' be produced, default value 10
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
          function(x, y, gate_x_plot=12 ,...) {
            if(gate_x_plot<=0) stop("gate_x_plot must be >0")
            nbits <- x@nbits
            ncbits <- x@circuit$ncbits
            n <- nbits + ncbits
            ngates <- length(x@circuit$gatelist)
            
            if(ngates > 0) {
              ipos <- rep(1, times=n)
              ## plot gates
              gatelist <- x@circuit$gatelist
              for(i in c(1:ngates)) {
                if(max(ipos)%%gate_x_plot==1){
                  ipos <- rep(1, times=n)
                  xmax <- gate_x_plot
                  ## prepare empty plot
                  plot(NA, ann=FALSE, xlim=c(0,xmax), ylim=c(0,n+1), axes=FALSE, frame.plot=FALSE)
                  ## plot qubit lines
                  for(ii in c(n:(ncbits+1))) {
                    lines(x=c(0.3, xmax), y=c(ii, ii))
                    annotate_bitnames(i=n-ii+1, y=ii, ...)
                  }
                  ## plot classical bit lines
                  if(ncbits > 0) {
                    for(ii in c(ncbits:1)) {
                      lines(x=c(0.3, xmax), y=c(ii-0.025, ii-0.025))
                      lines(x=c(0.3, xmax), y=c(ii+0.025, ii+0.025))
                      annotate_bitnames(i=ncbits-ii+1, y=ii, cbit=TRUE, ...)
                    }
                  }
                }
                
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
                
                # draw an arrow at the bottom left corner if the plot continue
                if(max(ipos)%%gate_x_plot==1 && i!=ngates ){
                  text(x=xmax-1.2, y=0.1, labels="continue")
                  arrows(x0=xmax-0.4, x1=xmax,
                         y0=0.1,
                         y1=0.1,
                         length=0.1)
                  #the newline command is required for pdf output
                  cat("\n\n")
                }
                
              }
            }
            else{
              ##emplty plot
              ipos <- rep(1, times=n)
              xmax <- gate_x_plot
              ## prepare empty plot
              plot(NA, ann=FALSE, xlim=c(0,xmax), ylim=c(0,n+1), axes=FALSE, frame.plot=FALSE)
              ## plot qubit lines
              for(ii in c(n:(ncbits+1))) {
                lines(x=c(0.3, xmax), y=c(ii, ii))
                annotate_bitnames(i=n-ii+1, y=ii, ...)
              }
              ## plot classical bit lines
              if(ncbits > 0) {
                for(ii in c(ncbits:1)) {
                  lines(x=c(0.3, xmax), y=c(ii-0.025, ii-0.025))
                  lines(x=c(0.3, xmax), y=c(ii+0.025, ii+0.025))
                  annotate_bitnames(i=ncbits-ii+1, y=ii, cbit=TRUE, ...)
                }
              }
            }
          
          }
          )
