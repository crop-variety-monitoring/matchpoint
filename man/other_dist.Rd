\name{other_dist}

\alias{other_dist}

\title{
Distance to other groups
}

\description{
Compute the mean (or other statistic) distance of the members of a group to the mean in a distance matrix. }


\usage{
other_dist(x, fun=mean)
}
 
\arguments{
\item{x}{square matrix or \code{dist}}
\item{fun}{summarizing function that accepts \code{na.rm=TRUE}}
}

\value{
data.frame with the group name and score(s)
}

\seealso{\code{\link{self_dist}}, \code{\link{nngb_dist}}}

\examples{
x <- iris[, -5]
d <- dist(x)
attr(d, "Labels") <- as.character(iris[,5])

other_dist(d, mean)
}

