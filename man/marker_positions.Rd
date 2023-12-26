\name{marker_positions}

\alias{marker_positions}

\title{
Get chromosome position data for SNP markers
}

\description{
Get chromosome position data for SNP markers
}


\usage{
marker_positions(crop="")
}
 
\arguments{
\item{crop}{character. One of \code{"bean", "cassava", "cowpea", "maize", "rice", "teff", "wheat"} or \code{""} for all these crops}
}

\value{
data.frame
}

\examples{
m <- marker_positions("cassava")
head(m)
}

