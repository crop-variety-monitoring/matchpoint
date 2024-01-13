\name{self_dist}

\alias{self_dist}

\title{
Within group distance
}

\description{
Compute the mean (or other statistic) distance between members of a group in a distance matrix. }


\usage{
self_dist(x, fun=mean)
}
 
\arguments{
\item{x}{square matrix or \code{dist}}
\item{fun}{summarizing function that accepts \code{na.rm=TRUE}}
}

\value{
data.frame with the group name and score(s)
}

\seealso{\code{\link{nngb_dist}}}

\examples{
x <- iris[, -5]
d <- dist(x)
attr(d, "Labels") <- as.character(iris[,5])

self_dist(d, mean)
self_dist(d, min)
self_dist(d, quantile)
}

