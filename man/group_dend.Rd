\name{group_dend}

\alias{group_dend}

\title{
Dendrogram for known groups 
}

\description{
Make a dendrogram from a distance matrix with prior (known) group labels. 
}

\usage{
group_dend(x, add=NULL, col=viridis::turbo, cex=0.5, ...)
}
 
\arguments{
\item{x}{square matrix or \code{dist}}
\item{add}{additons to the individual labels that are not used in setting the color}
\item{col}{colors or color generating function for the labels}
\item{cex}{numeric. Size of the labels}
\item{...}{arguments passed to \code{\link{hclust}}}
}

\value{
dendrogram
}


\examples{
i <- c(15,16,17,33,34,56,58,61,62,67,68,69,71,72,73,83,84,85,88,89,91,93,94,95,96,97,99,
		100,102,103,106,107,108,110,112,114,115,118,119,120,122,123,124,126,127,128,
		130,131,132,134,135,136,139,143,147,150)
x <- iris[i, ] 

dst <- dist(x[,-5])
attr(dst, "Labels") <- x[,5]

d1 <- group_dend(dst, col=c("red", "blue", "gray", "green", "orange"))
plot(d1, horiz=TRUE)

sdst <- split_groups(dst, 0.9, keep_nngb=FALSE, verbose=FALSE)
add <- ifelse(labels(sdst) != labels(dst), "  - changed!", "")
d2 <- group_dend(sdst, add, col=c("red", "blue", "gray", "green", "orange"))

par(mar=c(2,0,0,4), mfrow=c(1,2))
plot(d1, horiz=TRUE)
plot(d2, horiz=TRUE)
}

