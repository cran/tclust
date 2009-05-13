\name{geyser2}
\encoding{latin1}
\alias{geyser2}
\title{ Old Faithful Geyser Data }
\description{ 
A bivariate data set obtained from the Old Faithful Geyser, containing the eruption length and the length of the previous eruption for 271 eruptions of this geyser in minutes.
}
\usage{
geyser2
}
\format{
\describe{
\item{\code{Eruption length}}{ The eruption length in minutes. }
\item{\code{Previous eruption length}}{ The length of the previous eruption in minutes. }
}
}
\source{W Härdle

This particular data structure can be obtained using the following code:\cr
\code{f1 <- faithful[,1]}\cr
\code{geyser2 <- cbind (f1[-length(f1)], f1[-1])}\cr
\code{colnames (geyser2) <- c("Eruption length", "Previous eruption length")}
}
\references{ García-Escudero, L.A.; Gordaliza, A. (1999). "Robustness properties of k-means and trimmed k-means". Journal of the Americal Statistical Assoc., Vol. 94, No.447, 956-969.
Härdle, W. (1991) Smoothing Techniques with Implementation in S, New York: Springer.
}
