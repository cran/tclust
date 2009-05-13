tkmeans <-
function (x, k = 3, alpha = 0.05, niter = 50, Ksteps = 10, center, scale, store.x = TRUE, drop.empty.clust = TRUE, trace = 0, zero.tol = 1e-16)
{
	ret <- tclust (x = x, k = k, restr.fact = 1, restr = "eigen", alpha = alpha, niter = niter, Ksteps = Ksteps, equal.weights = TRUE, store.x = store.x, center = center, scale = scale, drop.empty.clust = drop.empty.clust, trace = trace, zero.tol = zero.tol)
	ret$obj = ret$cov[1,1,1]
	ret
}

