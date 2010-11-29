plot.ctlcurves <-
function(x, main, ylim, ylab, min.weights = FALSE, col, lty = 1, ...)
{
  if (min.weights)
  {
    dat <- x$min.weights
    dat.range <- range (dat[(x$par$k != 1),,])
    if (missing (ylab))
      ylab <- "Minimum Weigths"
  }
  else
  {
    dat.range <- range (dat <- x$obj)
    if (missing (ylab))
      ylab <- "Objective Function Value"
  }

    if (!missing (ylim))
    setylim <- ylim
  else #if (link.ylim )
    setylim <- ylim <- dat.range

  if (!missing (main))
    setmain <- main
    
  lty <- rep (lty, length (x$par$k))
  
  if (missing (col))
    col <- x$par$k + 1
  else
    col <- rep (col, length (x$par$k))

  for (i in 1: length (x$par$restr.fact))
  {
    if (missing (main))
      setmain <- paste ("Restriction Factor =", x$par$restr.fact[i])

    plot (0, type = "n", ylim = setylim, xlim = range (x$par$alpha),
          main = setmain, xlab = "", ylab = ylab)
    mtext ("a", side = 1, line = 3, font = 5)

    for (j in 1:length (x$par$k))
    {
            # k == 1 -> min.weights == 1 -> we're not interested in that. 
      if (min.weights && x$par$k[j] == 1)
        next    

      lines (x$par$alpha, dat[j,,i], type="b", col = col[j], lty = lty[j],
            pch = as.character (x$par$k[j]))
    }
  }
}

