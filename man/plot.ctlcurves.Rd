\name{plot.ctlcurves}
\alias{plot.ctlcurves}
\encoding{latin1}
\title{ Plot Method for ctlcurves Objects }
\description{
Plotting a \code{ctlcurves} object comparing the target functions values with different values of paramter \code{restr.fact}.
%These functions plot a ctlcurves object in two different styles, either classic or simple. The classical style displays three rows of diagnostic plots, analyzing different clustering results on the same dataset. The first row shows the target functions value. The second row displays the increase of the target functions value gained by increasing k by one, whereas the third row shows the minimum cluster weights for each estimation. The simple plots only display the target functions values, comparing the results of different values of parameter \code{restr.fact}.
}
\usage{
%plot.ctlcurves(x, plot.idx, link.ylim = TRUE, ylim, main, min.weights = FALSE, ...)
plot.ctlcurves(x, main, ylim, ylab, min.weights = FALSE, ...)
%plot.ctlcurves(x, simple, ...)
%plot.ctlcurves.simple(x, plot.idx, link.ylim = TRUE, ylim, main, ...)
%plot.ctlcurves.classic(x, plot.idx, set.idx = 1, ...)
}
\arguments{
  \item{x}{ The ctlcurves object to be printed. }
%  \item{simple}{ Logical value indicating whether the simple plot (\code{plot.ctlcurves.simple}) or the classical plot (\code{plot.ctlcurves.classic}) shall be displayed. If the simulation has been done for more than one factor value, by default the simple plot is displayed for each factor. If only one factor has been chosen, the classical plot is displayed by default. }
%  \item{plot.idx}{ An index vector % (\code{plot.ctlcurves.simple}) or matrix (\code{plot.ctlcurves.classic})
%  giving the indices of the plots which shall be displayed. If provided, no further plot arrangement (\code{\link[graphics]{par} (mfrow = ...)}) is done. }
%  \item{link.ylim}{ A logical value indicating whether the y-ranges of the plots shall be identical. Ignored, if parameter ylim is set. }
%  \item{ylim}{ The y-range of all plots displayed. If provided, parameter \code{link.ylim} is ignored. }
  \item{main}{ A character-string containing the title of the plot. }
  \item{ylim, ylab}{ Argument passed to \code{\link[graphics]{plot}}. }
%  \item{set.idx}{ If the \code{ctlcurves} - object contains more than one factor, only data related to \code{factor[set.idx]} is displayed. }
  \item{min.weights}{Locigal value indicating whether the minimum cluster weights shall be printed instead of the target functions values. By default this option is disabled}
   \item{\dots}{ Further arguments, passed to \code{\link[graphics]{plot}}. }
}
\details{ These curves show the values of the trimmed classification (log-)likelihoods when altering the trimming proportion \code{alpha} and the number of clusters \code{k}. The careful examination of these curves provides valuable information for choosing these parameters in a clustering problem. For instance, an appropriate \code{k} to be chosen is one that we do not observe a clear increase in the trimmed classification likelihood curve for k with respect to the k+1 curve for almost all the range of alpha values. Moreover, an appropriate choice of parameter alpha may be derived by determining where an initial fast increase of the trimmed classification likelihood curve stops for the final chosen k. A more detailed explanation can be found in García-Escudero et al (2009). }
%\value{ no value is returned }
\references{ García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering". Preprint available at www.eio.uva.es/infor/personas/langel.html. }
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
%\note{ ~~further notes~~   ~Make other sections like Warning with  ~}
%\section{Warning }{....}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
sig <- diag (2)
cen <- rep (1,2)
x <- rbind (
	rmvnorm (108, cen * 0,   sig),
	rmvnorm (162, cen * 5,   sig * 6 - 2),
	rmvnorm (30, cen * 2.5, sig * 50)
)

ctl.a <- ctlcurves (x, k = 1:4, restr.fact = 12)
plot(ctl.a)
}
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
