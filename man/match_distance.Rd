\name{match_distance}

\alias{match_distance}

\title{
Use a distance metric to match field samples to references
}

\description{
Use a distance metric with SNP data to match field samples to references. The Hamming distance is used with a two_row DaRT file. The Nei distance is used with a counts file
}

\usage{
match_distance(x, genotypes, missing_rate=0.25, comp_all=FALSE, filename="", verbose=TRUE)
}
 
 
\arguments{
\item{x}{\code{darter} object. See \code{\link{read_dart}}}

\item{genotypes}{\code{data.frame} with data linking the (dart) IDs (the column names in \code{x$snp}) to variety or field IDS. Column names must include \code{"sample", "variety", "reference"}. "sample" should match the ID for a genotype (column name) in \code{x$snp}. \code{"variety"} should have a name for the reference samples. \code{reference} has logical values that should be \code{TRUE} for the reference samples and \code{FALSE} for the field samples. Columns in \code{x$snp} with a name that is not matched by a value in \code{"sample"} are removed from the analysis}

\item{missing_rate}{numeric between 0 and 1. samples with a higher missing rate are removed}

\item{comp_all}{logical}

\item{filename}{character. Optional Excel output file name}

\item{verbose}{logical. If \code{TRUE} messages may be printed}
}

\value{
list
}

\examples{
fsnp <- system.file("ex/DCas00-0000_SNP_2row.csv", package="matchpoint")
fvar <- system.file("ex/DCas00-0000_variety-info.csv", package="matchpoint")

snps <- matchpoint::read_dart(fsnp)
genotypes <- data.table::fread(fvar)

out <- matchpoint::match_distance(snps, genotypes)

names(out)
head(out$best_match)
}
