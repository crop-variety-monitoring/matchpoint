\name{split_groups}

\alias{split_groups}

\title{
Split dissimilar groups
}

\description{
Split dissimilar groups in a distance matrix with prior (known) group labels. Splits are made if their are sub-groups that are more than \code{threshold} apart, unless, if \code{keep_nngb=TRUE} the nearest neighbors. 
}

\usage{
split_groups(x, threshold, keep_nngb=TRUE, verbose=FALSE)
}
 
\arguments{
\item{x}{square matrix or \code{dist}}
\item{threshold}{numeric. sub-groups that are further apart are split}
\item{keep_nngb}{logical. If \code{TRUE} subgroups that are nearest neighbor to each other are not split}
\item{verbose}{logical. If \code{TRUE} the names of the split varieties are printed}
}

\value{
same as \code{x} with possibly updated column/row names (if \code{x} is a matrix) or labels (if \code{x} is a \code{dist})
}

\seealso{\code{\link{nngb_dist}}, \code{\link{group_dend}}}

\examples{
x <- iris[, -5]
dst <- dist(x)
attr(dst, "Labels") <- as.character(iris[,5])

sdst <- split_groups(dst, 0.65, keep_nngb=FALSE, verbose=TRUE)
unique(labels(sdst))

# see group_dend for visualization
}

