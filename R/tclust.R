tclust <-
function (x, k = 3, alpha = 0.05, niter = 50, Ksteps = 10, restr = c("eigen", "deter", "sigma"), restr.fact = 2, equal.weights = FALSE, center, scale, store.x = TRUE, drop.empty.clust = TRUE, trace = 0, zero.tol = 1e-16)
{
	if (is.data.frame (x))
		x = as.matrix (x)
	if (!is.numeric (x))
		stop ("parameter x: numeric data frame / matrix / vector expected")

	if (!is.numeric (restr))
	{
		restr <- match.arg(restr)
		if (restr == "eigen")
			restr <- 0
		else if (restr == "deter")
			restr <- 1
		else if (restr == "sigma")
			restr <- 2
	}

	if( !is.matrix (x))
		x <- matrix (x, ncol = 1)

		
	p <- ncol (x)
	n <- nrow (x)

	if (missing (center))
		center = 0
	if (missing (scale))
		scale = 1

	x <- t( (t(x) - center) / scale)

#	scaled = pcaPP:::ScaleAdv (x, center, scale)

	ret <- .C ( "tclust",
				as.integer (c (dim (x), k, niter, Ksteps, equal.weights, restr, trace)),		# integer parameters (last params: m_dwDebugPrint, m_dwConvCount, m_dwIterSuccess, m_dwAllEVZero)
				parN = integer (3),
				as.double (c (alpha, restr.fact, zero.tol)),													# double parameters
				parD = double (1),
				as.double (x),
				center = double (p * k),
				cov = double (p * p * k),
				assign = integer (nrow (x)),
				clustsize = integer (k),
				PACKAGE = "tclust"
	)

	if (drop.empty.clust)
		idxuse = ret$clustsize > 0
	else
		idxuse = rep (TRUE, k)

	if (trace >= 1)
	{
		if (ret$parN[3] == 2)		## not all iterations could be executed
			if (restr == 1)
				warning ("points in the data set are concentrated in k subspaces after trimming") 
			else
				warning ("points in the data set are concentrated in k points after trimming") 

		if (ret$parN[2] && ret$parN[1] / ret$parN[2] < 0.5)
			warning (paste ("less than half of the iterations (", round (ret$parN[2] && ret$parN[1] / ret$parN[2] * 100, 1), "%) converged - increase Ksteps!", sep = ""))
	}
	if (trace >= 0)
		if (any (ret$clustsize[ret$clustsize != 0] < n / 50))
			warning ("clusters with size < n * 0.02 found - try reducing K")

	parlist <- list (k = k, alpha = alpha, restr.fact = restr.fact, niter = niter, Ksteps = Ksteps, equal.weights = equal.weights)

	if (store.x)
		parlist$x <- x

	ret <- list (
		center = array (ret$center, c(p, k)) [,idxuse, drop = FALSE], 
		cov = array (ret$cov, c (p, p, k)) [, , idxuse, drop = FALSE],
		assign = c(0, cumsum (idxuse)) [ret$assign + 2],
		par = parlist,
		k = sum (idxuse),
		obj = ret$parD[1],
		clustsize = ret$clustsize[idxuse],
		weights = ret$clustsize [idxuse] / floor (nrow (x) * (1 - alpha)),
		iter.successful = ret$parN[2],
		iter.converged = ret$parN[1],
		dim = dim (x)
	)

	cmm = matrix (scale, nrow = p, ncol = 1) %*% matrix (scale, nrow = 1, ncol = p)
	for (i in 1:ret$k)
	{
		ret$center[,i] <- ret$center[,i] * scale + center
		ret$cov[,,i] <- ret$cov[,,i] * cmm
	}

	class (ret) <- "tclust"

	return (ret)
}

