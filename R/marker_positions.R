
marker_positions <- function(crop="") {
	crop <- tolower(crop)
	if (crop != "") {
		stopifnot(crop %in% c('bean', 'cassava', 'cowpea', 'maize', 'rice', 'teff', 'wheat'))
	}
	x <- readRDS(system.file("ex/markers.rds", package="matchpoint"))
	if (crop != "") {
		x[x$crop == crop, ]
	} else {
		x
	}
}	



.prepare_markers <- function(filename, fcassava) {
	crops <- readxl::excel_sheets(filename)
	crops <- crops[crops != "cassava"]
	x <- lapply(crops, \(crop) {
		vars <- c("DART.ID", "Chr", "Pos")
		if (crop == "cowpea") {
			vars[1] <- "Alt .ID (optional)"
		}
		d <- data.frame(crop=crop, readxl::read_excel(filename, crop)[, vars])
		colnames(d)[2] <- "MarkerName"
		d
	})
	x <- do.call(rbind, x)
	cas <- matchpoint::read_dart(fcassava)
	pos <- data.frame(do.call(rbind, strsplit(cas$snp[,1], "_")))
	pos <- data.frame(crop="cassava", cas$snp[,1], pos[,-1])
	colnames(pos) <- c("crop", "MarkerName", "Chr", "Pos")
	x <- rbind(x, pos)
	x$MarkerName <- toupper(x$MarkerName)
	saveRDS(x, "matchpoint/inst/ex/markers.rds")
}
#filename <- "data-NGA/raw/dart//DarTAG Mid density panel info.xlsx"
#fcassava <- "data-NGA/raw/dart//cassava/Report_DCas23-7954_SNP.csv"

