plot.bayesfact <-
function (x, ...)
{
#	if (nrow (x$center) <= 2)
#	{
		par (mfrow = c(1, 3))
		plot (x$x)	
#	}
#	else
#		par (mfrow = c(1, 2))

	plot.bayesfact.p2 (x, ...)
	plot.bayesfact.p3 (x, ...)	
}

