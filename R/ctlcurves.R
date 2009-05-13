ctlcurves <-
function(x, k = 1:4, alpha = seq (0, 0.2, len = 6), restr.fact = 50, trace = 1, ...)
{
	if (length (restr.fact) > 1)
		restr.fact = restr.fact [1]
		
	if (trace == 1)
	{
		cat ("Depending on the dataset and parameters k.max, alpha.count and restr.fact, this\r\nfunction needs some time to compute. Increase parameter \x22trace\x22 for progress information. (Remove this message by setting trace = 0)\r\n")
		flush.console() 
	}

	obj <- min.weights <- array (NA, c (length (k), length (alpha), length (restr.fact)))

	for (i in 1:length(k))
	{
		for (j in 1:length (alpha))
		{
			cur.alpha = alpha[j]
			for (h in 1:length (restr.fact))
			{
				if (trace >= 2)
					cat ("k =", k[i], "; alpha =", alpha[j], "; restr.fact =", restr.fact[h])

				clus <- tclust (x, k = k[i], alpha = alpha[j], restr.fact = restr.fact[h], trace = -1, ...)
				obj[i,j, h] <- clus$obj							##	the objective function criterion
				min.weights[i,j, h] <- min(clus$weights)		##	weights of the smallest group

				if (trace >= 2)
					cat ("; obj = ", clus$obj, ";min (weights) =", min(clus$weights), "\r\n")
			}
		}
	}

	par <- list (x = x, k = k, alpha = alpha, restr.fact = restr.fact)
	ret <- list (obj = obj, min.weights = min.weights, par = par)
	class (ret) = "ctlcurves"
	ret
}

