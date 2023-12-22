
get_info <- function(path) {
	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	fld <- as.data.frame(readxl::read_excel(ffld, 2))[, c("crop", "dart.id", "planting.barcode")]
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	colnames(fld)[3] <- "field.id"
	fld$variety <- ""
	fld$type <- "survey"
	
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	ref <- as.data.frame(readxl::read_excel(fref, 2))[, c("crop", "dart.id", "var.name.full")]
	colnames(ref)[3] <- "variety"
	ref$field.id <- ""
	ref$type <- "reference"
	
	fldref <- rbind(fld, ref)
	fldref$crop <- tolower(fldref$crop)
	fldref
	
}


prepare_dart <- function(path) {

	outpath <- gsub("raw/dart", "clean", path)
	intpath <- gsub("raw/dart", "intermediate", path)

	dir.create(outpath, FALSE, FALSE)
	info <- get_info(path)

	crops <- list.dirs(path, full.names=FALSE)
	crops <- crops[crops != ""]
	crops <- crops[!grepl("/", crops)]

	floc <- list.files(intpath, pattern="location", full.names=TRUE)
	loc <- NULL
	if (length(floc) == 1) {
		loc <- as.data.frame(readxl::read_excel(floc))
	}
	fpan <- file.path(path, "DarTAG Mid density panel info.xlsx")


	for (crop in crops) {
		print(crop); flush.console()
		croppath <- file.path(path, crop)
		cropf <- list.files(croppath, pattern="SNP.csv$", full.names=TRUE, ignore.case=TRUE)
		stopifnot (length(cropf) == 1)
		oname <- gsub("Report_", "", basename(cropf))
		oname <- gsub("_SNP.csv", "", oname)
		
		x <- read_dart(cropf) 
		x <- make_dart_1row(x)
		x$type <- NULL

		if (crop == "cassava") {
			position <- data.frame(do.call(rbind, strsplit(x$snp[,1], "_")))
			position <- cbind(x$snp[,1], position[,-1])
			colnames(position) <- c("MarkerName", "Chr", "Pos")
		} else {
			vars <- c("DART.ID", "Chr", "Pos")
			if (crop == "cowpea") {
				vars[1] <- "Alt .ID (optional)"
			}
			position <- as.data.frame(readxl::read_excel(fpan, crop))[, vars]
			colnames(position)[1] <- "MarkerName"
			position$MarkerName <- toupper(position$MarkerName)
			x$marker$MarkerName <- toupper(x$marker$MarkerName)
			i <- position$MarkerName %in% x$marker$MarkerName
			position <- position[i, ]
			stopifnot(nrow(position) == nrow(x$marker))
			m <- match(position$MarkerName, x$marker$MarkerName)
			if (any(is.na(m))) stop("unknown markers")
			position <- position[m,]
		}
		x$marker$order <- 1:nrow(x$marker)
		m <- merge(position, x$marker, by="MarkerName", all.y=TRUE)
		m <- m[m$order, ]
		m$order <- NULL
		x$marker <- m

		inf <- info[info$crop == crop, ]
		if (!is.null(loc)) {
			n <- nrow(inf)
			i <- match(inf$field.id, loc$field.id)
			inf <- merge(inf, loc, by="field.id", all.x=TRUE)
			stopifnot(nrow(inf) == n)
			inf <- inf[, c(names(info), c("longitude", "latitude"))]
		}
		x$info <- inf

		ibs <- x$snp[,-1]
		cinfo <- info[info$crop == crop, ]
		i <- match(colnames(ibs), cinfo$dart.id)
		unk <- colnames(ibs)[which(is.na(i))]
		if (length(unk) > 0) {
			message(paste0("removing unknown samples:\n", paste(unk, collapse=", ")))
			i <- i[!is.na(i)]		
		}

		j <- info$type[i] == "reference"
		ibs_ref <- cbind(position, ibs[, j])
		ibs_fld <- cbind(position, ibs[, !j])

		write.csv(ibs_ref, file.path(outpath, paste(crop, oname, "ibs_ref.csv", sep="_")), row.names=FALSE)
		write.csv(ibs_fld, file.path(outpath, paste(crop, oname, "ibs_fld.csv", sep="_")), row.names=FALSE)

		write.csv(inf, file.path(outpath, paste(crop, oname, "info.csv", sep="_")), row.names=FALSE)

		writexl::write_xlsx(x, file.path(outpath, paste0(crop, "_", oname, ".xlsx")), format_headers=FALSE)		
	}
}


