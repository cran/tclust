plot.tclust.1d <-
function (x, tol = 0.95, jitter.y = FALSE, col, main, sub, labels = c("cluster", "observation"), text, xlab, ylab = "", ...)
{
	if (x$dim[2] != 1)
		stop ("tclust object of dimension 1 expected.")
	
	if (missing (col))
		col = x$assign + 1
	
	if (is.null (x$par$x))
		stop ("dataset not included in tclust object - cannot plot object.")
	dn <- dimnames (x$par$x)
	
	if (missing (xlab))
		if (is.list (dn) && length (dn[[2]]) == 1)
			xlab = dn[[2]][1]
		else
			xlab = "x1"
	
	if (missing (sub))
		sub = paste ("k = ", x$k, #" (", x$par$k, "), 
					", alpha = ", x$par$alpha, #", obj = ", round (x$obj, 2),
					sep = "")

	if (missing (main))
		#main = "Cluster Assignment"
		main = "Classification"

	if (jitter.y)
		y <- runif (x$dim[1], min = -1) / 2.5
	else
		y <- rep (0, x$dim[1])

	if (!missing (text))
	{
		if(length (text) != length (x$par$x))
		{
			labels = NULL
			warning (paste ("parameter text: text array of length", length (x$par$x), "expected"))
		}
	}
	else if (!missing (labels))
	{
		labels <- match.arg(labels)
		if (labels == "cluster")
			text = paste (x$assign)
		else if (labels == "observation")
			text = paste (1:length (x$par$x))
	}
	
	if (missing (text))
		plot (x$par$x, y, xlab = xlab, ylab = ylab, main = main, col = col, ylim = c(-1,1), axes = FALSE, ...)
	else
	{
		plot (x$par$x, y, xlab = xlab, ylab = ylab, main = main, type = "n", ylim = c(-1,1), axes = FALSE, ...)
		text (x$par$x, y, labels = text, col = col)
	}
	
	mtext(sub, cex = 0.8, line= 0.25)
	
	box ()
	axis (side =1)

	x.c = as.numeric(x$center)
	x.sd = sqrt (as.numeric (x$cov))
	
	.vline (x.c, 3, lty = 2, col = 1 + (1:x$k))
	
	if (is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
	{
		tol.fact = sqrt(qchisq(tol, 1))
		.vline (x.c + x.sd * tol.fact, 2, col = (1:x$k) + 1)
		.vline (x.c - x.sd * tol.fact, 2, col = (1:x$k) + 1)
	}
}

