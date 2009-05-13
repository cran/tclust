plot.tclust.1d <-
function (x, labels = c("cluster", "observation"), text, main, sub, xlab, ylab, pch, col, by.cluster = FALSE, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col, jitter.y = FALSE, ...)
{
	if (x$dim[2] != 1)
		stop ("tclust object of dimension 1 expected.")
	
	if (is.null (x$par$x))
		stop ("dataset not included in tclust object - cannot plot object.")

	if (by.cluster)
	{
		maxassig <- max (x$assign)
		
		if (missing (col))
			col <- 1:(x$k+1)
		else
			col <- rep (col,  len = maxassig + 1 )	
		
		if (missing (pch))
			pch <- 1:(x$k+1)
		else
			pch <- rep (pch,  len = maxassig + 1 )
		col <- col [x$assign + 1]		
		pch <- pch [x$assign + 1]
	}
	else
	{
		if (missing (col))
			col <- x$assign + 1
		if (missing (pch))
			pch <- 1
	}

	dn <- dimnames (x$par$x)
	if (missing (xlab))
		if (is.list (dn) && length (dn[[2]]) == 1)
			xlab = dn[[2]][1]
		else
			xlab = "x1"
	if (missing (ylab))
		ylab <- ""
	
	if (missing (sub))
		sub = paste ("k = ", x$k, #" (", x$par$k, "), 
					", alpha = ", x$par$alpha, #", obj = ", round (x$obj, 2),
					sep = "")

	if (missing (main))
		main = "Classification"

	if (jitter.y)
		y <- runif (x$dim[1], min = -1) / 2.5
	else
		y <- rep (0, x$dim[1])

	if (!missing (text))
	{
		text <- rep (text, length (x$par$x))
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
		plot (x$par$x, y, xlab = xlab, ylab = ylab, main = main, ylim = c(-1,1), axes = FALSE, pch = pch, col = col, ...)
	else
	{
		plot (x$par$x, y, xlab = xlab, ylab = ylab, main = main, type = "n", ylim = c(-1,1), axes = FALSE, pch = pch, ...)
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
		if (missing (tol.col))
			tol.col <- (1:x$k) + 1
		else
			tol.col <- rep (tol.col, x$k)
			
		tol.lty <- rep (tol.lty , x$k)
		tol.lwd <- rep (tol.lwd , x$k)
			
			
		.vline (x.c + x.sd * tol.fact, 2, col = tol.col, lwd = tol.lwd, lty = tol.lty)
		.vline (x.c - x.sd * tol.fact, 2, col = tol.col, lwd = tol.lwd, lty = tol.lty)
	}
}

