\name{match_IBS}

\alias{match_IBS}

\title{
Match field samples to reference with a dart 1-row SNP file and IBS
}

\description{
Use Identity-By-State (IBS) to match SNP data from field samples to references
}

\usage{
match_IBS(SNPs, genotypes, markers, MAF_cutoff=0.05, SNP_Missing_Rate=0.2, 
		Ref_Missing_Rate=0.2, Sample_Missing_Rate=0.2, 
		Ref_Heterozygosity_Rate=0.1, Sample_Heterozygosity_Rate=0.1,
		IBS_cutoff=0.5, Inb_method="mom.visscher", 
		threads=1, verbose=FALSE, filename="")
}
 
 
\arguments{
\item{SNPs}{\code{data.frame} with the DArT SNP marker data. Rows are genotypes. The first column has should be called \code{"MarkerName"} and have the SNP marker names. The other columns have the marker values that DArT uses: \code{0, 1, 2, NA} . 
}

\item{genotypes}{\code{data.frame} with data linking the (dart) IDs (the column names in \code{SNPs}) to variety or field IDS. Column names must include \code{"plate.id", "field.id", "genotype", "reference"}. "plate.id" should match the ID for a genotype ("TargetID") in the \code{SNPs} data. \code{"genotype"} should have a value (name) for the reference samples. \code{reference} has logical values that should be \code{TRUE} for the reference samples and \code{FALSE} for the field samples. Columns in \code{SNPs} with a name that is not matched by \code{"plate.id"} are removed from the analysis}

\item{markers}{\code{data.frame} with at least these three columns: \code{"MarkerName", "Chr", "Pos"} that store chromosome number and position for each SNP. You can get these with \code{\link{marker_positions}} for a number of crops}

\item{MAF_cutoff}{numeric. Minor Allele Frequency cutoff. Only use makers with a MAF of at least this value. A SNP minor allele frequency of 0.05 is often used in the literature}

\item{SNP_Missing_Rate}{numeric. Markers with a higher missing rate are not used. A SNP missing rate 0.9 is a low threshold, 0.2 is a high threshold}

\item{Ref_Missing_Rate}{numeric. Reference missing rate threshold}

\item{Ref_Heterozygosity_Rate}{numeric. Reference heterozygosity rate threshold}

\item{Sample_Heterozygosity_Rate}{numeric. Reference heterozygosity rate threshold}

\item{Sample_Missing_Rate}{numeric. Sample missing rate threshold }

\item{IBS_cutoff}{numeric. Cutoff}

\item{Inb_method}{character. Inbreeding method. One of "mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"}

\item{threads}{positive integer. Number of threads to use}

\item{verbose}{logical}

\item{filename}{character. If not empty, the the filename is used to write a gds and xlsx file}
}

\value{
list
}

\examples{
fsnp <- system.file("ex/DCas00-0000_SNP.csv", package="matchpoint")
fvar <- system.file("ex/DCas00-0000_genotype-info.csv", package="matchpoint")
snps <- matchpoint::read_dart(fsnp)
genotypes <- data.table::fread(fvar)
markers <- matchpoint::marker_positions("")
out <- match_IBS(snps$snp, genotypes, markers)

names(out)
head(out$best_match)
}
