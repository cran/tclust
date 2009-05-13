\name{plot.bayesfact}
\encoding{latin1}
\alias{plot.bayesfact}
\alias{plot.bayesfact.p2}
\alias{plot.bayesfact.p3}
\title{ Plot Method for bayesfact Objects }
\description{
  Next to a plot of the \code{tclust} object which has been used for calculating the \code{bayesfact} object, a silhouette plot indicates the amount of doubtfully assigned observations. A third plot similar to the standard \code{tclust} plot highlights the identified doubtful observations.
}
\usage{
plot.bayesfact(x, ...)
plot.bayesfact.p2(x, xlab = "Bayes Factor", ylab = "Clusters", main, xlim, print.bayes = TRUE, ...)
plot.bayesfact.p3(x, main = "Doubtful Assignments", ...)
}
\arguments{
  \item{x}{ An object of class "\code{bayesfact}" as from \code{bayesfact ()}. }
  \item{xlab, ylab, xlim}{ Arguments passed to funcion \code{\link[tclust]{plot.tclust}}. }
  \item{main}{ Argument passed to funcion \code{\link[graphics]{plot}}. }
%  \item{draw.legend}{ A logical value indicating whether the legend shall be displayed. }
  \item{print.bayes}{ A logical value indicating whether each clusters mean Bayes factor shall be plotted }
  \item{\dots}{ Further arguments are passed to function \code{\link[graphics]{plot}} }
}
\details{
\code{plot.bayesfact.p2} displays a silhouette plot based on the Bayes factors of the observations. A solution with many large Bayes factors is not reliable. Such clusters can be identified with this silhouette plot. %Thus, when comparing clustering solutions, small Bayes factors are preferred. 
Thus \code{plot.bayesfact.p3} displays the dataset, highlighting observations with Bayes factors greater than the given threshold. Function \code{plot.bayesfact} combines the standard plot of a \code{tclust} object, and the two plots introduced here.
}
\value{ No return value is provided.}
\references{ \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering". Preprint available at www.eio.uva.es/infor/personas/langel.html. }
\author{ Agustin Mayo Iscar, Luis Angel \enc{García}{Garcia} Escudero, Heinrich Fritz }
%\note{ ~~further notes~~ ~Make other sections like Warning with \section{Warning }{....} ~}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
sig <- diag (2)
cen <- rep (1,2)
x <- rbind (
	rmvnorm (360, cen * 0,   sig),
	rmvnorm (540, cen * 5,   sig * 6 - 2),
	rmvnorm (100, cen * 2.5, sig * 50)
)

clus.1 <- tclust (x, k=2, alpha=0.1, restr.fact=12)
clus.2 <- tclust (x, k=3, alpha=0.1, restr.fact=1) # "Improper" choice of parameters

bay.1 <- bayesfact(clus.1)
plot(bay.1)

bay.2 <- bayesfact(clus.2)
plot(bay.2)

dev.off ()
plot.bayesfact.p2(bay.1)
plot.bayesfact.p3(bay.2)
}
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
