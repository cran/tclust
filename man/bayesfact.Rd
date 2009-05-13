\name{bayesfact}
\encoding{latin1}
\alias{bayesfact}
\alias{print.bayesfact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayes Factor Analysis for tclust Objects}
\description{ Analyzes a \code{tclust}-object by calculating Bayes factors and comparing the quality of the actual cluster assignments and the second best possible assignment for each observation. Bayes factors, measuring the strength of the "trimming" decision may be defined (!? what does this sentence mean?). Cluster assignments of observations with large Bayes factors are considered as "doubtful" decisions. Silhouette plots give a graphical overview of the Bayes factors distribution (see \code{\link{plot.bayesfact}}). More details can be found in García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009).}
\usage{
bayesfact(x, threshold = 1/10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{tclust} object. }
  \item{threshold}{ %Observations with Bayes factors larger than \code{log (threshold)} are considered as to be "doubtfully" assigned to a cluster.
  A cluster assignment or a trimming decision for an observation with a Bayes factor larger than \code{log(threshold)} is considered as a "doubtful" decision. 
   }
}
\details{
This function compares the actual (best) assignment of each observation to its second best possible assignment. This comparison is based on the bayes factors of each observation, which are calculated here. If the Bayes factor of an observation is larger than a given level (\code{log (threshold)}), the observation is considered as "doubtfully" assigned to a cluster. More detailed information is shown when \link[=plot.bayesfact]{plotting} the returned \code{baysefact} object.
%  The function plots two graphs: The first displays a silhouette plot based on the Bayes factors of the observations. The cluster solution is not reliable when groups with many observations with large Bayes factors appear in the silhouette plot. Thus, when comparing clustering solutions, those having smaller values for the associated mean of Bayes factor are preferable. The second plot displays the dataset, highlighting observations with Bayes factor greater than the given threshold.
}
\value{
The function returns an S3 object of type \code{baysefact} containing the following components:
	\item{x}{ A \code{tclust} object. }
	\item{ylimmin}{ A minimum y-limit calculated for plotting purposes. }
	\item{ind}{ The actual cluster assignment. }
	\item{ind2}{ The second most likely cluster assignment for each observation. }
	\item{disc}{ The (weighted) likelihood of the actual cluster assignment of each observation. }
	\item{disc2}{ The (weighted) likelihood of the second best cluster assignment of each observation. }
	\item{assignfact}{ The factor \code{log (disc/disc2)}. }
	\item{threshold}{ The threshold used for deciding whether \code{assignfact} indicates a "doubtful" assignment. }
	\item{mean.bayesfact}{ A vector of length \code{k + 1} containing the mean Bayes factors for each cluster (including the outliers). }
}
\references{ García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering". Preprint available at www.eio.uva.es/infor/personas/langel.html. }
\author{ Agustin Mayo-Iscar, Luis Angel García-Escudero, Heinrich Fritz }
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
plot(bay.1)

bay.2 <- bayesfact(clus.2)
plot(bay.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
