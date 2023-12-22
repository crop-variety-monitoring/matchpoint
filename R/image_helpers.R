
get_info <- function(path) {

	ffld <- list.files(path, pattern="identifier", full.names=TRUE)[1]
	stopifnot(file.exists(ffld))
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
		colnames(x$marker)[colnames(x$marker) == "AlleleID"] <- "MarkerName"
		colnames(x$snp)[colnames(x$snp) == "AlleleID"] <- "MarkerName"
		cn <- colnames(x$snp)
		tab <- table(cn)
		dups <- tab[tab > 1]
		if (length(dups) > 0) {
			for (i in 1:length(dups)) {
				n <- dups[i]
				stopifnot(n == 2)
				nm <- names(n)
				cn[cn == nm] <- paste0(nm, c("_R1", "_R2"))
			}
			colnames(x$snp) <- cn
		}
		
		x$marker$MarkerName <- toupper(x$marker$MarkerName)
		x$snp$MarkerName <- toupper(x$snp$MarkerName)

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
			i <- position$MarkerName %in% x$marker$MarkerName
			#j <- (!x$marker$MarkerName %in% position$MarkerName )
			position <- position[i, ]
			if (nrow(position) != nrow(x$marker)) {
				message(paste("SNP data has", nrow(x$marker), "markers. Panel has", nrow(position), "matches")) 
			}
			m <- match(x$marker$MarkerName, position$MarkerName)
			if (any(is.na(m))) message("unknown markers found")
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

		ibs <- merge(x$marker[, 1:3], x$snp, all=TRUE)
		cinfo <- info[info$crop == crop, ]
		mnames <- cinfo$dart.id
#		prefix <- paste0(oname, "_")
#		prefix <- gsub("_moreOrders", "", prefix)
#		nopref <- which(substr(mnames, 1, 3) != substr(oname, 1, 3)) 
#		if (length(nopref) > 0) {
#			mnames[nopref] <- paste0(prefix, mnames[nopref])
#		}
		cns <- colnames(ibs)[-c(1:3)]
		cns <- gsub("_R.$", "", cns)
		i <- match(cns, mnames)
		unk <- cns[which(is.na(i))]
		if (length(unk) > 0) {
			n <- length(unk)
			if (n > 8) {
				unk <- c(unk[1:6], "...")
			}
			message(paste0("removing ", n, " unknown samples:\n", paste(unk, collapse=", ")))
			i <- i[!is.na(i)]		
		}
		
		j <- which(info$type[i] == "reference") + 3

		ibs_ref <- ibs[, j]
		ibs_fld <- ibs[, -j]

		write.csv(ibs_ref, file.path(outpath, paste0(oname, "_ibs-ref.csv")), row.names=FALSE)
		write.csv(ibs_fld, file.path(outpath, paste0(oname, "_ibs-fld.csv")), row.names=FALSE)

		write.csv(inf, file.path(outpath, paste0(oname, "_info.csv")), row.names=FALSE)

		writexl::write_xlsx(x, file.path(outpath, paste0(oname, ".xlsx")), format_headers=FALSE)		
	}
}


