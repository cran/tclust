tclust <-
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 10, restr = c("eigen",
          "deter", "sigma"), restr.fact = 12, equal.weights = FALSE, center,
          scale, store.x = TRUE, drop.empty.clust = TRUE, trace = 0,
          zero.tol = 1e-16)
{

  usetrace = max (0, trace)
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

#  scaled = pcaPP:::ScaleAdv (x, center, scale)

  ret <- .C ( "tclust",
        as.integer (c (dim (x), k, nstart, iter.max, equal.weights, restr, 
                    usetrace)),
        parN = integer (3),
        as.double (c (alpha, restr.fact, zero.tol)),
        parD = double (1),
        as.double (x),
        center = double (p * k),
        cov = double (p * p * k),
        cluster = integer (nrow (x)),
        size = integer (k),
        PACKAGE = "tclust", DUP = FALSE
  )

  if (drop.empty.clust)
    idxuse <- which (ret$size > 0)
  else
    idxuse <- 1:k

  idxuse <- idxuse [order (ret$size[idxuse], decreasing = TRUE)]

  if (trace >= 1)
  {
    if (ret$parN[3] == 2)    ## not all iterations could be executed
      if (restr == 1)
        warning ("points in the data set are concentrated in k subspaces after trimming") 
      else
        warning ("points in the data set are concentrated in k points after trimming") 

    if (ret$parN[2] && ret$parN[1] / ret$parN[2] < 0.5)
      warning (paste ("less than half of the iterations (", round (ret$parN[2]
               && ret$parN[1] / ret$parN[2] * 100, 1), 
               "%) converged - increase iter.max!", sep = ""))
  }
  if (trace >= 0)
    if (any (ret$size[ret$size != 0] < n / 50))
      warning ("clusters with size < n * 0.02 found - try reducing K")

  parlist <- list (k = k, alpha = alpha, restr.fact = restr.fact, nstart = 
                    nstart, iter.max = iter.max, equal.weights = equal.weights)

  if (store.x)
    parlist$x <- x

  ClusterIDs <- rep (0, k + 1)
  ClusterIDs[idxuse + 1] <- 1:length (idxuse)
  ret <- list (
    centers = array (ret$center, c(p, k)) [,idxuse, drop = FALSE], 
    cov = array (ret$cov, c (p, p, k)) [, , idxuse, drop = FALSE],
    #cluster = c(0, cumsum (idxuse)) [ret$cluster + 2],
    cluster = ClusterIDs [ret$cluster + 2], # + 2, as in C++ trimmed = -1
    par = parlist,
    k = length (idxuse),
    obj = ret$parD[1],
    size = ret$size[idxuse],
    weights = ret$size [idxuse] / floor (nrow (x) * (1 - alpha)),
    iter.successful = ret$parN[2],
    iter.converged = ret$parN[1],
    dim = dim (x)
  )

  cmm = matrix (scale, nrow = p, ncol = 1) %*% matrix (scale, nrow = 1,
                ncol = p)
  for (i in 1:ret$k)
  {
    ret$centers[,i] <- ret$centers[,i] * scale + center
    ret$cov[,,i] <- ret$cov[,,i] * cmm
  }

  class (ret) <- "tclust"

  return (ret)
}
