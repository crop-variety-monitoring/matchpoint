\name{match_IBS}

\alias{match_IBS}

\title{
Use Identity-By-State with SNP data to match field samples to references
}

\description{
Use Identity-By-State (IBS) with SNP data to match field samples to references
}

\usage{
match_IBS(x, genotypes, markers, MAF_cutoff=0.05, SNP_Missing_Rate=0.2, 
		Ref_Missing_Rate=0.2, Sample_Missing_Rate=0.2, 
		Ref_Heterozygosity_Rate=1, Sample_Heterozygosity_Rate=1,
		IBS_cutoff=0.5, Inb_method="mom.visscher", 
		threads=4, verbose=FALSE, filename="")
}
 
 
\arguments{
\item{x}{\code{darter} object. See \code{\link{read_dart}}}

\item{genotypes}{\code{data.frame} with data linking the (dart) IDs (the column names in \code{x$snp}) to variety or field IDS. Column names must include \code{"sample", "variety", "reference"}. "sample" should match the ID for a genotype (column name) in \code{x$snp}. \code{"variety"} should have a name for the reference samples. \code{reference} has logical values that should be \code{TRUE} for the reference samples and \code{FALSE} for the field samples. Columns in \code{x$snp} with a name that is not matched by a value in \code{"sample"} are removed from the analysis}

\item{markers}{\code{data.frame} with at least these three columns: \code{"MarkerName", "Chr", "Pos"} that store chromosome number and position for each SNP. You can get these with \code{\link{marker_positions}} for a number of crops}

\item{MAF_cutoff}{numeric. Minor Allele Frequency cutoff. Only use makers with a MAF of at least this value. A SNP minor allele frequency of 0.05 is often used in the literature}

\item{SNP_Missing_Rate}{numeric. Markers with a higher missing rate are not used. The default value of 0.2 is a stringent threshold}

\item{Ref_Missing_Rate}{numeric. Reference missing rate threshold}

\item{Ref_Heterozygosity_Rate}{numeric. Currently ignored! Reference heterozygosity rate threshold}

\item{Sample_Heterozygosity_Rate}{numeric. Currently ignored! Reference heterozygosity rate threshold}

\item{Sample_Missing_Rate}{numeric. Sample missing rate threshold }

\item{IBS_cutoff}{numeric. Cutoff}

\item{Inb_method}{character. Inbreeding method. One of "mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"}

\item{threads}{positive integer. Number of threads to use to speed up the computation (four threads may be optimal)}

\item{verbose}{logical. If \code{TRUE} messages may be printed}

\item{filename}{character. If not empty, the the filename is used to write a gds and xlsx file}
}

\value{
list
}

\examples{
fsnp <- system.file("ex/DCas00-0000_SNP.csv", package="matchpoint")
fvar <- system.file("ex/DCas00-0000_variety-info.csv", package="matchpoint")

snps <- matchpoint::read_dart(fsnp)
genotypes <- data.table::fread(fvar)
markers <- matchpoint::marker_positions("")

out <- matchpoint::match_IBS(snps, genotypes, markers)

names(out)
head(out$best_match)
}
