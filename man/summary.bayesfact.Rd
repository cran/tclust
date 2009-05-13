\name{summary.bayesfact}
\encoding{latin1}
\alias{summary.bayesfact}
\title{ Summary Method for class bayesfact}
\description{ The summary method for class \code{bayesfact}.
}
\usage{
summary.bayesfact (object, hide.emtpy = TRUE, show.clust, show.alt, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class "\code{bayesfact}" as from \code{bayesfact ()}. }
  \item{ hide.emtpy }{ A logical value specifying whether clusters without doubtful assignment shall be hidden. }
  \item{ show.clust }{ A logical value specifying whether the number of doubtful assignments per cluster shall be displayed. }
  \item{ show.alt }{ A logical value specifying whether the alternative cluster assignment shall be displayed. }
  \item{\dots}{ Arguments passed to or from other methods. }
}
%\details{}
\value{ No return value is provided.}
\references{ García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering". Preprint available at www.eio.uva.es/infor/personas/langel.html. }
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{ \code{\link{plot.bayesfact}}}
\examples{
sig <- diag (2)
cen <- rep (1,2)
x <- rbind (
	rmvnorm (360, cen * 0,   sig),
	rmvnorm (540, cen * 5,   sig * 6 - 2),
	rmvnorm (100, cen * 2.5, sig * 50)
)
clus.1 <- tclust (x, k=2, alpha=0.1, restr.fact = 12)
clus.2 <- tclust (x, k=3, alpha=0.1, restr.fact = 1)  # "Improper" choice of parameters

bay.1 <- bayesfact(clus.1)
summary(bay.1)

bay.2 <- bayesfact(clus.2)
summary(bay.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
