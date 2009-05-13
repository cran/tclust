plot.DiscrFact.p3 <-
function (x, main = "Doubtful Assignments", pch, col, col.nodoubt = grey (0.8), by.cluster = FALSE, ...)
{
	idxplot = x$assignfact > x$threshold
	n = length (x$x$assign)
	
	if (by.cluster)
	{
		maxassig <- max (x$x$assign)
		
		if (missing (col))
			col <- 1:(x$x$k+1)
		else
			col <- rep (col,  len = maxassig + 1 )	

		col [idxplot] <- col [x$ind[idxplot]+1]
		col [!idxplot] = col.nodoubt

		if (missing (pch))
			pch <- 1:(x$x$k+1)
		else
			pch <- rep (pch,  len = maxassig + 1 )
		pch <- pch [x$x$assign + 1]		
	}
	else
	{
		if (missing (col))
		{
			col <- rep (col.nodoubt, n)
			col [idxplot] <- x$ind[idxplot]+1 
		}
	
		if (missing (pch))
		{
			pch <- rep (1, n)
			pch [idxplot] = 19
		}
	}
	
	if (nrow (x$x$center) == 2)
		plot (x$x, pch = pch, col = col, tol.col = 1, main = main, sub = "", ...)
	else
		plot (x$x, pch = pch, col = col, main = main, sub = "", ...)
}

