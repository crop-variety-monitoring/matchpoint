\name{lump_similar}

\alias{lump_similar}

\title{
Lump cases
}

\description{
Lump cases
}

\usage{
lump_similar(x, threshold, verbose=FALSE)
}
 
\arguments{
\item{x}{square matrix or \code{dist}}
\item{threshold}{numeric. items that less nearer than this may be lumped}
\item{verbose}{logical. If \code{TRUE} the names of the split varieties are printed}
}

\value{
same as \code{x} with possibly updated column/row names (if \code{x} is a matrix) or labels (if \code{x} is a \code{dist})
}

\seealso{\code{\link{split_groups}}}

\examples{
x <- iris[, -5]
dst <- dist(x)
attr(dst, "Labels") <- as.character(iris[,5])

sdst <- lump_similar(dst, 0.65, keep_nngb=FALSE, verbose=TRUE)
unique(labels(sdst))

}

