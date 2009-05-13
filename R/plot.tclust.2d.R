plot.tclust.2d <-
function (x, tol = 0.95, col, labels = c("cluster", "observation"), main, sub, tol.col = 1, tol.lty = 3, text, xlab, ylab, ...)
{

	if (nrow (x$center) != 2)
		stop ("tclust object of dimension 2 expected.")

	if (missing (col))
		col = x$assign + 1

	if (is.null (x$par$x))
		stop ("dataset not included in tclust object - cannot plot object.")
	
	dn <- dimnames (x$par$x)
	if (is.list (dn) && length (dn[[2]]) == 2)
	{
		if (missing (xlab))
			xlab = dn[[2]][1]
		if (missing (ylab))
			ylab = dn[[2]][2]
	}
	else
	{
		if (missing (xlab))
			xlab = "x1"
		if (missing (ylab))
			ylab = "x2"
	}
	
	if (missing (sub))
		sub = paste ("k = ", x$k, #" (", x$par$k, "), 
					", alpha = ", x$par$alpha, #", obj = ", round (x$obj, 2), 
					sep = "")

	if (missing (main))
		#main = "Cluster Assignment"
		main = "Classification"

	x1 <- x$par$x[,1]
	x2 <- x$par$x[,2]

	if (!missing (text))
	{
		if(length (text) != length (x1))
			warning (paste ("parameter text: text array of length", length (x1), "expected"))
	}
	else if (!missing (labels))
	{
		labels <- match.arg(labels)
		if (labels == "cluster")
			text = paste (x$assign)
		else if (labels == "observation")
			text = paste (1:length (x1))
	}

	if (missing (text))
		plot (x1, x2, col = col, main = main, xlab = xlab, ylab = ylab, ...)
	else
	{
		plot (x1, x2, col = col, main = main, type = "n", xlab = xlab, ylab = ylab, ...)
		text (x1, x2, labels = text, col = col)
	}
	mtext(sub, cex = 0.8, line= 0.25)

	if (is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
	{
		if (length (tol.col) == 1)
			tol.col <- rep (tol.col, x$k)
		else if (length (tol.col) != x$k)
		{
			warning (paste ("color vector of length", x$k, ", or single value expected."))
			tol.col <- 1 + (1:x$k)
		}

		tol.fact = sqrt(qchisq(tol, 2))	
		for (k in 1:x$k)
				.doEllipses (eigen = eigen (x$cov[,,k]), center = x$center [,k], lty = tol.lty, col = tol.col [k], size = tol.fact)
	}
}

