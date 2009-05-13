plot.bayesfact.p2 <-
function (x, xlab = "Bayes Factor", ylab = "Clusters", main, xlim, print.bayes = TRUE, ...)
#, draw.legend
{
	n = x$x$dim[1]

	if (missing (main))
		#main = paste ("Mean Overall Bayes Factor =", format (mean (x$assignfact), digits = 3))
		main = "Silhouette Plot"	

	if (missing (xlim))
			xlim = c(x$ylimmin,0)
			
#	plot.new ()
#	par (usr = c (.stretch (xlim, 0.05), .stretch (c(0, n), 0.05) ))
	plot (0, 0, xlim = xlim,ylim = c(1,n), type="n", xlab = xlab, ylab = ylab, main = main, axes = FALSE, ...)
	

	axis (side = 1)

	cs = c (0, cumsum (c (x$x$dim[1] - sum (x$x$clustsize), x$x$clustsize)))

	{
		ylines <- cs[-1]
		ylines <- ylines [-length (ylines)]
		abline (h = ylines, lty = 2)
	}

	cs = (cs [-1] + cs[-length(cs)]) / 2
	axis (side = 2, at = cs, labels = c("O", 1:x$x$k))
	box ()

	cury = 0
	for (k in 0:x$x$k)
	{
		grupo.k <- sort(x$assignfact[x$ind==k])
		gs = length (grupo.k)
		if (gs > 0)
			polygon ( c(0, grupo.k, 0), c (1, 1:gs, gs) + cury, border = 0, col = k + 1)
		
		{	
			ll <- cury
			ul <- cury + gs

			if (k == 0)
				ll <- par("usr") [3]
			if (k == x$x$k)
				ul <- par("usr") [4]

#			lines (rep (x$mean.bayesfact[k + 1], 2), c(ll, ul), lty = 3)
		}
		cury = cury + gs
	}

	xpos = sum (par ("usr")[1:2] * c(1, 4)) / 5
	
	if (print.bayes)
	{
#		text (xpos, cs, paste ("mean BF:", format (x$mean.bayesfact, digits = 4)), adj = 1)
		legend ("topleft", legend = format (x$mean.bayesfact[(x$x$k + 1):1], digits = 4), inset = 0.04, col = 1 + (x$x$k:0), pch = 15, title = "Mean Bayes Factors", box.lwd = 0, bty = "n")	
	}
	abline(v = x$threshold + 1, lty = 2)

#	if (draw.legend)
#		legend (xlim[1] * 0.95, x$x$dim[1], c("Threshold", "Mean Bayes Factor per Cluster"), lty = 2:3)
}

