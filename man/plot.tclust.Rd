\name{plot.tclust}
\encoding{latin1}
\alias{plot.tclust}
\alias{plot.tclust.1d}
\alias{plot.tclust.2d}
\alias{plot.tclust.Nd}
\title{ Plot Method for tclust Objects }
\description{
  Plots a \code{tclust} object created by \code{\link[tclust]{tclust}}. 
}
\usage{
plot.tclust(x, ...)
plot.tclust.1d(x, labels = c("cluster", "observation"), text, main, sub, xlab, ylab, pch, col, by.cluster = FALSE, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col, jitter.y = FALSE, ...)
plot.tclust.2d(x, labels = c("cluster", "observation"), text, main, sub, xlab, ylab, pch, col, by.cluster = FALSE, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col = 1, ...)
plot.tclust.Nd(x, labels = c("cluster", "observation"), text, main, sub, xlab, ylab, pch, col, by.cluster = FALSE, ...)
}
\arguments{
  \item{x}{ The \code{tclust} object to be printed. }
  \item{labels}{ A string specifying the type of labels to be drawn. Either "cluster" or "observation" can be specified. If specified, parameter \code{pch} is ignored. If omitted, the symbol, specified by \code{pch} is used. }  
  \item{text}{ A vector of length \code{n} (the number of observations) containing strings which are used as labels for each observation. If specified, the parameters \code{text} and \code{pch} are ignored. }  
  \item{main, sub, xlab, ylab, pch, col}{ Arguments passed to \code{\link[graphics]{plot}}. }
  \item{by.cluster}{ Logical value indicating whether parameters \code{pch} and \code{col} refer to observations (\code{FALSE}) or clusters (\code{TRUE}). }
  \item{tol}{ The tolerance interval. 95\% tolerance ellipsoids (assuming normality) are plotted by default. }
  \item{jitter.y}{ Logical value, specifying whether the drawn values shall be jittered in y-direction for better visibility of structures in 1 dimensional data. }
  \item{tol.lwd}{ A vector of length \code{k} or \code{1} containing the line widths of the tolerance ellipsis/lines.}
  \item{tol.lty}{ A vector of length \code{k} or \code{1} containing the line type of the tolerance ellipsis/lines.}
  \item{tol.col}{ A vector of length \code{k} or \code{1} containing the color of the tolerance ellipsis/lines.}
  \item{\dots}{ Parameters passed to \code{\link[graphics]{plot}}-function. }
}
\details{ \code{plot.tclust} selects the appropriate function according to the dimensionality of the \code{tclust} object. One and two dimensional structures are treated separately (e.g. tolerance intervals/ellipses are displayed). Higher dimensional structures are displayed by plotting the two first Fisher's canonical coordinates (evaluated by \code{\link[tclust]{discr_coords}}) and derived from the final cluster assignments (trimmed observations are not taken into account). \code{plot.tclust.Nd} can be called with one or two-dimensional \code{tclust}-objects too. The function fails, if "\code{store.x = FALSE}" is specified in \code{\link[tclust]{tclust}}, because the original data matrix is required here. }
\value{ No return value is provided. }
%\references{ ~put references to the literature/web site here ~ }
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~ }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
#--- EXAMPLE 1------------------------------
sig <- diag (2)
cen <- rep (1,2)
x <- rbind (
	rmvnorm (360, cen * 0,   sig),
	rmvnorm (540, cen * 5,   sig * 6 - 2),
	rmvnorm (100, cen * 2.5, sig * 50)
)
# Two groups and 10\% trimming level
a <- tclust (x, k = 2, alpha = 0.1, restr.fact = 12)
plot (a)
plot (a, labels = "observation")
plot (a, labels = "cluster")
plot (a, by.cluster = TRUE)

#--- EXAMPLE 2------------------------------
x <- c(	rnorm(72,	0,   1),
		rnorm(108, 10,   1),
		rnorm(20,  2.5, 10))
%a <- tclust (x[,1], k = 2, alpha = 0.1, restr.fact = 12)
a <- tkmeans (x, k = 2, alpha = 0.1)
plot (a, jitter.y = TRUE)

}
\keyword{ hplot }
\keyword{ multivariate }
\keyword{ robust }
\keyword{ cluster }