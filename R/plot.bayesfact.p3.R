plot.bayesfact.p3 <-
function (x, main = "Doubtful Assignments", ...)
{
	idxplot = x$assignfact > x$threshold
	n = length (x$x$assign)
	col <- rep (8, n)
	col [idxplot] <- x$ind[idxplot]+1 ## was x$ind2 [.....

	pch <- rep (1, n)
	pch [idxplot] = 19

	if (nrow (x$x$center) == 2)
		plot (x$x, pch = pch, col = col, tol.col = 1, main = main, sub = "", ...)
	else
		plot (x$x, pch = pch, col = col, main = main, sub = "", ...)
}

