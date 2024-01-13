\name{nngb_dist}

\alias{nngb_dist}

\title{
Within group distance
}

\description{
Compute the mean (or other statistic) distance between members of a group and their nearest neighbor from another group in a distance matrix. }


\usage{
nngb_dist(x, fun=min)
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

nngb_dist(d, min)
nngb_dist(d, quantile)
}

