
tclust <-
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, 
          restr = c ("eigen", "deter", "sigma"), restr.fact = 12, 
          equal.weights = FALSE, center, scale, store.x = TRUE, 
          drop.empty.clust = TRUE, trace = 0, warnings = 3, zero.tol = 1e-16
         )
{
#  Disabled arguments:
#  restr = c ("eigen", "deter", "sigma", "dir.eigen", "dir.deter", "prop"), 
#  iter.tune, 
  ovv <- 0      ##  optimization parameter for "optVectors"
  fuzzy <- FALSE
  m <- 2

  usetrace <- .tr (trace, 0)

  x <- .Conv2Matrix (x, substitute(x))

  if (!is.numeric (x))
    stop ("parameter x: numeric matrix/vector expected")

  rd <- .match.restr (restr)
  restr.C <- rd$restr
  deter.C <- rd$deter

#  if (missing (iter.tune))
    iter.tune <- 10
  iter.tune <- rep (iter.tune, len = 3)

  p <- ncol (x)
  n <- nrow (x)

  if (missing (center))
    center = 0
  if (missing (scale))
    scale = 1

#  x <- t( (t(x) - center) / scale)

  scaled <- .ScaleAdv (x, center, scale)
  x <- scaled$x
  scale <- scaled$scale
  center <- scaled$center

  z.size <- ifelse (fuzzy, n * k, 0)

  ret.C <- .C ( "tclust",
        as.integer (c (dim (x), k, fuzzy, nstart, iter.max, equal.weights, 
                                   restr.C, deter.C, usetrace, iter.tune, ovv)),
        parN = integer (5),
        as.double (c (alpha, restr.fact, m, zero.tol)),
        parD = double (2),
        as.double (x),
        center = double (p * k),
        cov = double (p * p * k),
        cluster = integer (nrow (x)),
        size = double (k),
        weights = double (k),
        z = double (z.size),
        er.obj = double (nstart),
        er.conv = integer (nstart),
        PACKAGE = "tclust", DUP = FALSE
  )

  err.exec <- ret.C$parN[4]	## error execution

  if (err.exec)
    stop ()


  parlist <- list (k = k, alpha = alpha, restr.fact = restr.fact, nstart = 
                   nstart, iter.max = iter.max, equal.weights = equal.weights,
                   restr = restr[1], restr.C = restr.C, deter.C = deter.C)

  if (store.x)
    parlist$x <- x

  if (drop.empty.clust)
    idxuse <- which (ret.C$size > 0)
  else
    idxuse <- 1:k

  idxuse <- idxuse [order (ret.C$size[idxuse], decreasing = TRUE)]

  ClusterIDs <- rep (0, k + 1)
  ClusterIDs[idxuse + 1] <- 1:length (idxuse)

  k.real <- length (idxuse)

  int <- list (
      iter.converged = ret.C$parN[1]
    , iter.successful = ret.C$parN[2]
    , dim = dim (x)
    , code = ret.C$parN[3]
    , frac.restr.ok = ret.C$parN[5] / nstart
    , er.obj = ret.C$er.obj
    , er.conv = ret.C$er.conv
  )

  ret <- list (
    centers = array (ret.C$center, c(p, k)) [,idxuse, drop = FALSE]
    , cov = array (ret.C$cov, c (p, p, k)) [, , idxuse, drop = FALSE]
    , cluster = ClusterIDs [ret.C$cluster + 2]
    , par = parlist
    , k = sum (ret.C$size > 0) #k.real
    , obj = ret.C$parD[1]
    , unrestr.fact = ceiling (ret.C$parD[2])
    , size = ret.C$size[idxuse]
    , weights = ret.C$weights[idxuse]
    , int = int
  )

  dn.x <- dimnames (x)
  dn.var <- if (is.null (dn.x[[2]])) paste ("X", 1:ncol (x)) else dn.x[[2]]
  dn.clus <- paste ("C", 1:k.real)

  dimnames (ret$centers) <- list (dn.var, dn.clus)
  dimnames (ret$cov) <- list (dn.var, dn.var, dn.clus)

    ## generate a list of warnings
  ret$warnings <- list ( singular = ret$int$code == 2
                , iter = ret$int$iter.successful && ret$int$iter.converged / 
                         ret$int$iter.successful < 0.5
                , drop = parlist$k > ret$k
                , size = any (ret$size < n / 50)
                , sizep = min (ret$size) <= p
                , smallobj = ret$obj < (-1e+20)
                , restr.lo = ret$unrestr.fact > restr.fact
                , restr.hi = ret$unrestr.fact * 2 < restr.fact
                )
    ##  no error message that the restriction factor is too low if we 
	##    have too small groups.
  if (ret$warnings$sizep || ret$warnings$size)    
    ret$warnings$restr.lo <- FALSE

  .tclust.warn (warnings, ret)

  if (fuzzy)
    ret$z <- matrix (ret.C$z, n, k)[,idxuse]

  cmm <- matrix (scale, nrow = p, ncol = 1) %*%
         matrix (scale, nrow = 1, ncol = p)

  ret$mah <- array (dim = n)
  ret$mah[!ret$cluster] <- NA

  for (i in 1:ret$k)
  {
    idx <- ret$cluster == i
    ret$mah[idx] <- mahalanobis (x[idx, , drop = FALSE], center = ret$centers[, i], cov = ret$cov[,, i])

    ret$centers[,i] <- ret$centers[,i] * scale + center
    ret$cov[,,i] <- ret$cov[,,i] * cmm
  }

  class (ret) <- "tclust"

  return (ret)
}

.tclust.warn <- function (warnings, ret)
{
  if (warnings >= 1)
  {
    if (ret$warnings$iter)
      warning (paste ("Less than 50% of the iterations (", 
        round (ret$int$iter.converged / ret$int$iter.successful * 100, 1), 
               "%) converged - please increase iter.max.", sep = ""))
    if (ret$warnings$smallobj)
      warning ("Due to a very small objective function's value, the final solution does not seem to be reliable.\n  More iterations are probably needed (-> increase \"nstart\").")
  }

  if (warnings >= 2)
  {
    if (ret$warnings$singular)    ## not all iterations could be executed
      if (ret$par$deter.C)
        warning ("All observations are concentrated in k subspaces after trimming.")
      else
        warning ("All observations are concentrated in k points after trimming.")
#        warning ("Points in the data set are concentrated in k subspaces after trimming.")
#      else
#        warning ("points in the data set are concentrated in k points after trimming")

    if (ret$warnings$drop)
    {
      n.drop <- ret$par$k - ret$k
      if (n.drop == 1)
          warning (paste (n.drop, "empty cluster has been detected - try reducing k."))
      else
          warning (paste (n.drop, "empty clusters have been detected - try reducing k."))
    }
    else if (ret$warnings$size)
      warning ("Clusters with size < n * 0.02 found - try reducing k.")
    else if (ret$warnings$size)
     warning ("Clusters with size <= p found - try reducing k.")

  }
  if (warnings >= 3)
  {
    if (ret$par$restr != "sigma")
    {
      if (ret$warnings$restr.lo)
         warning (paste ("The result is artificially constrained due to ", 
         "restr.fact = ",ret$par$restr.fact,".", sep = ""))
#       warning (paste ("The chosen restriction factor (", ret$par$restr.fact,
#       ") artificially restricts the solution.\n",
#       "  This solution implies a restriction factor of ",
#       ceiling (ret$unrestr.fact), ".", sep = ""))

#      if (ret$warnings$restr.hi)    ##  warning currently disabled..
#        warning (paste ("The restriction factor (", ret$par$restr.fact,
#        ") has been chosen too large for this solution.\n",
#        "  This solution implies a restriction factor of ",
#        ceiling (ret$unrestr.fact), ".", sep = ""))
    }
  }
}

.match.restr <- function (restr)
{
  stopifnot (is.character (restr))

  restr <- restr[1]
  restr <- match.arg (restr, c ("eigen", "deter", "sigma", "dir.eigen",
                               "dir.deter", "prop", "none"))
  deter <- (restr == "deter" || restr == "dir.deter")

  if (restr == "eigen" || restr == "deter")
    restr <- 0
  else if (restr == "dir.eigen" || restr == "dir.deter")
    restr <- 1
  else if (restr == "sigma")
    restr <- 2
  else if (restr == "prop")
    restr <- 3
  else if (restr == "none")
    restr <- 4
  else
    stop ("unknown value for parameter restr")

  list (restr = restr, deter = deter)
}

#convplot <- function (x)
#{  ##  function for evaluating the convergence of a run of tclust
#   ord <- order (x$int$er.obj)
#   plot (x$int$er.obj[ord], col = 2-x$er.conv[ord],
#   xlab = "Iteration (ordered)",
#   ylab = "Objective Function", main = "Convergence Plot")
#}
