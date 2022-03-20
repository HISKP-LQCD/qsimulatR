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
            gpos <- c(1:ngates)
            if(ngates > 0) {
              gatelist <- x@circuit$gatelist
              for(i in c(1:ngates)) {
                if(is.na(gatelist[[i]]$bits[2])) {
                  gpos[i] <- ipos[gatelist[[i]]$bits[1]]
                  ipos[gatelist[[i]]$bits[1]] <- ipos[gatelist[[i]]$bits[1]] + 1
                }
                else {
                  gpos[i] <- max(ipos)
                  ipos[1:n] <- max(ipos) + 1
                }
              }
              for(i in c(1:ngates)) {
                ii <- which(gpos == i)
                ni <- length(ii)
                moved <- c()
                if(ni > 1) {
                  for(j in c(ni:2)) {
                    for(k in c((j-1):1)) {
                      if(!(j %in% moved || k %in% moved)) {
                        ## one is a multi-qubit gate
                        if(!is.na(gatelist[[ii[j]]]$bits[2]) || is.na(gatelist[[ii[k]]]$bits[2])) {
                          rj <- range(gatelist[[ii[j]]]$bits, na.rm=TRUE)
                          if(rj[1] < gatelist[[ii[k]]]$bits[1] && gatelist[[ii[k]]]$bits[1] < rj[2]) {
                            gpos[gpos > i] <- gpos[gpos > i] + 1
                            gpos[ii[j]] <- gpos[ii[j]] + 1
                            moved <- c(moved, j)
                          }
                        }
                        ## the other is a multi-qubit gate
                        else if(is.na(gatelist[[ii[j]]]$bits[2]) || !is.na(gatelist[[ii[k]]]$bits[2])) {
                          rk <- range(gatelist[[ii[k]]]$bits, na.rm=TRUE)
                          if(rk[1] < gatelist[[ii[j]]]$bits[1] && gatelist[[ii[j]]]$bits[1] < rk[2]) {
                            gpos[gpos > i] <- gpos[gpos > i] + 1
                            gpos[ii[k]] <- gpos[ii[k]] + 1
                            moved <- c(moved, k)
                          }
                        }
                        ## both are multi-qubit gates
                        else if(!is.na(gatelist[[ii[j]]]$bits[2]) || !is.na(gatelist[[ii[k]]]$bits[2])) {
                          rj <- range(gatelist[[ii[j]]]$bits, na.rm=TRUE)
                          rk <- range(gatelist[[ii[k]]]$bits, na.rm=TRUE)
                          if(any(duplicated(c(c(rj[1]:rj[2]),
                                              c(rk[1]:rk[2]))))) {
                            gpos[gpos > i] <- gpos[gpos > i] + 1
                            gpos[ii[j]] <- gpos[ii[j]] + 1
                            moved <- c(moved, j)
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            xmax <- max(gpos) + 1
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
              ## plot gates
              gatelist <- x@circuit$gatelist
              for(i in c(1:ngates)) {
                ## single qubit gates
                if(is.na(gatelist[[i]]$bits[2])) {
                  type <- gatelist[[i]]$type
                  if(gatelist[[i]]$type == "Rx" || gatelist[[i]]$type == "Ry" || gatelist[[i]]$type == "Rz") {
                    type <- paste0(gatelist[[i]]$type, "(", format(gatelist[[i]]$angle, digits=3), ")") 
                  }
                  legend(x=gpos[i], y=n+1-gatelist[[i]]$bits[1],
                         type, xjust=0.5, yjust=0.5,
                         x.intersp=-0.5, y.intersp=0.1,
                         bg="white")
                }
                ## multi qubit gates
                else {
                  xp <- gpos[i]
                  if(!is.null(gatelist[[i]]$controlled)) {
                    if(gatelist[[i]]$controlled) {
                      type <- gatelist[[i]]$type
                      if(gatelist[[i]]$type == "Rx" || gatelist[[i]]$type == "Ry" || gatelist[[i]]$type == "Rz") {
                        type <- paste0(gatelist[[i]]$type, "(", format(gatelist[[i]]$angle, digits=3), ")") 
                      }
                      ncbits <- length(gatelist[[i]]$bits) - 2
                      for(k in 1:ncbits){
                        pch <- 19
                        if(!is.null(gatelist[[i]]$inverse)) if(gatelist[[i]]$inverse[k]) pch <- 1
                        points(x=xp, y=n+1-gatelist[[i]]$bits[k], pch=pch, cex=1.5)
                        lines(x=c(xp, xp), y=n+1-c(gatelist[[i]]$bits[k], gatelist[[i]]$bits[k+1]))
                      }
                      legend(x=xp, y=n+1-gatelist[[i]]$bits[ncbits+1],
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
