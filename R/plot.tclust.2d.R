plot.tclust.2d <-
function (x, labels = c ("cluster", "observation"), text, main, sub, xlab, ylab,
          pch, col, by.cluster = TRUE, tol = 0.95, tol.lwd = 1, tol.lty = 3,
          tol.col = 1, main.pre, ...)
{
  if (nrow (x$centers) != 2)
    stop ("tclust object of dimension 2 expected.")

  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  if (by.cluster)
  {
    maxassig <- max (x$cluster)

    if (missing (col))
      col <- 1:(x$k+1)
    else
      col <- rep (col,  len = maxassig + 1 )  

    if (missing (pch))
      pch <- 1:(x$k+1)
    else
      pch <- rep (pch,  len = maxassig + 1 )

    col <- col[x$cluster + 1]    
    pch <- pch[x$cluster + 1]
  }    
  else
  {
    if (missing (col))
      col <- x$cluster + 1
    if (missing (pch))
      pch <- 1
  }

  dn <- dimnames (x$par$x)
  if (is.list (dn) && length (dn[[2]]) == 2)
  {
    if (missing (xlab))
      xlab = dn[[2]][1]
    if (missing (ylab))
      ylab = dn[[2]][2]
  }
  else
  {
    if (missing (xlab))
      xlab = "x1"
    if (missing (ylab))
      ylab = "x2"
  }
  
  if (missing (sub))
    sub <- bquote(paste (k == .(x$par$k), ", ", alpha == .(x$par$alpha)))
#     sub = paste ("k = ", x$k, #" (", x$par$k, "), 
#           ", alpha = ", x$par$alpha, #", obj = ", round (x$obj, 2), 
#           sep = "")

  if (missing (main))
    #main = "Cluster Assignment"
    main <- "Classification"
  if (!missing (main.pre) && !is.null (main.pre))
    main <- paste (main.pre, main)

  x1 <- x$par$x[,1]
  x2 <- x$par$x[,2]

  if (!missing (text))
  {
    text <- rep (text, length (x1))
  }
  else if (!missing (labels))
  {
    labels <- match.arg(labels)
    if (labels == "cluster")
      text = paste (x$cluster)
    else if (labels == "observation")
      text = paste (1:length (x1))
  }

  if (missing (text))
    plot (x1, x2, main = main, xlab = xlab, ylab = ylab, pch = pch, col = col,
          ...)
  else
  {
    plot (x1, x2, col = col, main = main, type = "n", xlab = xlab, ylab = ylab,
          pch = pch, ...)
    text (x1, x2, labels = text, col = col)
  }
  mtext(sub, cex = 0.8, line= 0.25)

  if (is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
  {
    tol.col <- rep (tol.col, x$k)
    tol.lty <- rep (tol.lty, x$k)
    tol.lwd <- rep (tol.lwd, x$k)
      

    tol.fact = sqrt(qchisq(tol, 2))  
    for (k in 1:x$k)
        .doEllipses (eigen = eigen (x$cov[,,k]), center = x$centers[,k],
        lwd = tol.lwd, lty = tol.lty[k], col = tol.col[k], size = tol.fact)
  }
}

