plot.tclust <-
function (x,  ...)
{
	if (x$dim[2] == 1)
		plot.tclust.1d (x, ...)
	else if (x$dim[2] == 2)
		plot.tclust.2d (x, ...)
	else
		plot.tclust.Nd (x, ...)
}

