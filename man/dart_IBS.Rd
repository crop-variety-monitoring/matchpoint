\name{dart_IBS}

\alias{dart_IBS}

\title{
Match field samples to reference with a dart 1-row SNP file and IBS
}

\description{
Match field samples to reference with a dart 1-row SNP file and IBS
}

\usage{
dart_IBS(dart, info, MAF_cutoff=0.05, SNP_Missing_Rate=0.2, 
		Ref_Missing_Rate=0.2, Sample_Missing_Rate=0.2, 
		Ref_Heterozygosity_Rate=0.1, Sample_Heterozygosity_Rate=0.1,
		IBS_cutoff=0.5, Inb_method="mom.visscher", 
		threads=1, verbose=FALSE, filename="")
}
 
\arguments{
\item{dart}{data.frame with the genetic reference data. Rows are genotypes. Columns have the markers. First three columns are "MarkerName", "Chr", "Pos"}

\item{info}{data.frame with data linking the (dart) IDs in the \code{dart} to variety or field IDS. Columns must include "dart.id", "field.id", "variety", "reference"}

\item{MAF_cutoff}{numeric. Minor Allele Frequency cutoff. Only use makers with a MAF of at least this value. A SNP minor allele frequency of 0.05 is often used in the literature}

\item{SNP_Missing_Rate}{numeric. Markers with a higher missing rate are not used. A SNP missing rate 0.9 is a low threshold, 0.2 is a high threshold}

\item{Ref_Missing_Rate}{numeric. Reference missing rate threshold}

Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate,

\item{Sample_Missing_Rate}{numeric. Sample missing rate threshold }

\item{IBS_cutoff}{numeric. Cutoff}

\item{inb_method}{character. Inbreeding method. One of "mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"}

\item{threads}{positive integer. Number of threads to use}

\item{verbose}{logical}

\item{filename}{character. If not empty, the the filename is used to write a gds and xlsx file}
}

\value{
list
}

\examples{
fdart <- system.file("ex/dart.csv", package="matchpoint")
finfo <- system.file("ex/info.csv", package="matchpoint")

dart <- data.frame(data.table::fread(fdart), check.names=FALSE)
info <- data.table::fread(finfo)

out <- matchpoint:::dart_IBS(dart, info)
}

