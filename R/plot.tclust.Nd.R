plot.tclust.Nd <-
function (x, labels = c ("cluster", "observation"), text, main, sub, xlab, ylab,
          pch, col, by.cluster = TRUE, main.pre, ...)
{
  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")
  
  if (missing (sub))
    sub <- paste ("k = ", x$k, #" (", x$par$k, "), 
            ", alpha = ", x$par$alpha, #", obj = ", round (x$obj, 2), 
            sep = "")

  if (by.cluster)
  {
    maxassig <- max (x$cluster)
    
    if (missing (col))
      col <- 1:(x$k+1)
    else
      col <- rep (col,  len = maxassig + 1 )  

    if (missing (pch))
      pch <- 1:(x$k+1)#rep (1, x$k+1)  #
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
            
  if (missing (main))
    #main = "Cluster Assignment"
    main <- "Classification"

  if (!missing (main.pre) && !is.null (main.pre))
    main <- paste (main.pre, main)    

  X <- discr_coords (x, x$par$equal.weights)

  if (!missing (text))
  {
    text <- rep (text, x$dim[1])
  }
  else if (!missing (labels))
  {
      labels <- match.arg(labels)
    if (labels == "cluster")
      text = paste (x$cluster)
    else if (labels == "observation")
      text = paste (1:nrow (X))
  }

  if (missing (xlab))
    xlab <- "First discriminant coord."
  if (missing (ylab))
    ylab <- "Second discriminant coord."

  if (missing (text))
    plot (X[,1],X[,2], main = main, axes = FALSE, xlab = xlab, ylab =
    ylab, pch = pch, col = col, ...)
  else
  {
    plot (X[,1], X[,2], main = main, type = "n", axes = FALSE, xlab = xlab,
          ylab= ylab, pch = pch, ...)
    text (X[,1], X[,2], labels = text, col = col)
  }

  mtext(sub, cex = 0.8, line= 0.25)
  
  box ()
}

