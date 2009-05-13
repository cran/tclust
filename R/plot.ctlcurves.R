plot.ctlcurves <-
function(x, main, ylim, ylab, min.weights = FALSE, ...)
{
#	if (missing (plot.idx))
#		par (mfrow = c (1, length (x$par$restr.fact)))

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
#	if (missing (mfrow))
#		mfrow = c(1,length (x$par$restr.fact))
#	if (!is.null (mfrow))
#		par (mfrow = mfrow)

	if (!missing (ylim))
		setylim = ylim
	else #if (link.ylim )
		setylim <- ylim <- dat.range

	if (!missing (main))
		setmain = main

	for (i in 1: length (x$par$restr.fact))
	{
#		if (!missing (plot.idx) && !any (i == plot.idx))
#			next

		if (missing (main))
			setmain = paste ("Restriction Factor =", x$par$restr.fact[i])

#		if (missing (ylim))
#			setylim = range (dat[,,i]) 	

		plot (0, type = "n", ylim = setylim, xlim = range (x$par$alpha), main = setmain, xlab = "alpha", ylab = ylab)

		for (j in 1:length (x$par$k))
		{
			if (min.weights && x$par$k[j] == 1)
				next		# k == 1 -> min.weights == 1 -> we're not interested in that. 

			lines (x$par$alpha, dat[j,,i], type="b", col = x$par$k[j] + 1, lty = 1, pch = as.character (x$par$k[j]))
		}
	}
}

