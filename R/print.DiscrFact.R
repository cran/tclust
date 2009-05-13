print.DiscrFact<-
function (x, ...)
{
	cat ("Mean overall discriminant factor:", mean (x$assignfact), "\r\n")
#	cat ("\r\nFurther information on this DiscrFact - object goes here\r\n")

	cat ("Mean discriminant factor per cluster:\r\n")
	print (x$mean.DiscrFact)

	idx = x$assignfact > x$threshold

	if (!sum (idx))
	{
		cat ("No decision is considered as doubtful\r\n")
		return
	}

	cat (sum (idx), "decisions are considered as doubtful\r\n")
}

