print.bayesfact <-
function (x, ...)
{
	cat ("Mean overall bayes factor:", mean (x$assignfact), "\r\n")
#	cat ("\r\nFurther information on this bayesfact - object goes here\r\n")

	cat ("Mean bayes factor per cluster:\r\n")
	print (x$mean.bayesfact)

	idx = x$assignfact > x$threshold

	if (!sum (idx))
	{
		cat ("No decision is considered as doubtful\r\n")
		return
	}

	cat (sum (idx), "decisions are considered as doubtful\r\n")
}

