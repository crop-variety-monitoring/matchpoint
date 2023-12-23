
get_info <- function(path) {

	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	ffld <- ffld[substr(basename(ffld), 1, 2) != "~$"]
	stopifnot(file.exists(ffld))
	fld <- as.data.frame(readxl::read_excel(ffld, 2))[, c("crop", "dart.id", "planting.barcode")]
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	colnames(fld)[3] <- "field.id"
	fld$variety <- ""
	fld$type <- "survey"
	
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	fref <- fref[substr(basename(fref), 1, 2) != "~$"]
	ref <- as.data.frame(readxl::read_excel(fref, 2))[, c("crop", "dart.id", "var.name.full")]
	colnames(ref)[3] <- "variety"
	ref$field.id <- ""
	ref$type <- "reference"
	
	fldref <- rbind(fld, ref)
	fldref$crop <- tolower(fldref$crop)
	unique(fldref)	
}


prepare_dart <- function(path) {

	outpath <- gsub("raw/dart", "clean", path)
	intpath <- gsub("raw/dart", "intermediate", path)

	dir.create(outpath, FALSE, FALSE)
	varinfo <- get_info(path)

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
			message("adding _D1 _D2 to duplicates")
			for (i in 1:length(dups)) {
				n <- dups[i]
				stopifnot(n == 2)
				nm <- names(n)
				cn[cn == nm] <- paste0(nm, c("_D1", "_D2"))
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
			#j <- !(x$marker$MarkerName %in% position$MarkerName )
			position <- position[i, ]
			if (nrow(position) != nrow(x$marker)) {
				message(paste("SNP data has", nrow(x$marker), "markers. Panel has", nrow(position), "matches")) 
			}
			m <- match(x$marker$MarkerName, position$MarkerName)
			if (any(is.na(m))) message(paste0(sum(is.na(m)), " unknown markers found"))
		}
		x$marker$order <- 1:nrow(x$marker)
		m <- merge(position, x$marker, by="MarkerName", all.y=TRUE)
		m <- m[m$order, ]
		m$order <- NULL
		x$marker <- m

		inf <- varinfo[varinfo$crop == crop, ]
		if (!is.null(loc)) {
			n <- nrow(inf)
			i <- match(inf$field.id, loc$field.id)
			inf <- merge(inf, loc, by="field.id", all.x=TRUE)
			stopifnot(nrow(inf) == n)
			inf <- inf[, c(names(inf), c("longitude", "latitude"))]
		}
		x$info <- inf

		ibs <- merge(x$marker[, 1:3], x$snp, all=TRUE)
		
		refnms <- as.character(inf[inf$crop == crop, "dart.id"])
		dartnms <- gsub("_D.$", "", colnames(ibs)[-c(1:3)])
		i <- match(dartnms, refnms)
		
		unk <- dartnms[which(is.na(i))]
		if (length(unk) > 0) {
			if (oname == "DMz23-7957_7228") {
				## undeclared 
				unk <- unk[!(unk %in% c('327', '328', '329', '330', '331', '332', '333', '334', '335', '336', '337', '338', '339', '340', '341', '342', '343', '344'))]
			} else if (oname ==  "DCob23-7823") {
				unk <- unk[!(unk %in% c('A5145', 'A5146', 'A5147', 'A5148', 'A5149', 'A5150', 'A5151', 'A5152', 'A5153', 'A5154', 'A5155', 'A5156', 'A5157', 'A5158', 'A5159', 'A5160', 'A5161', 'A5162', 'A5163', 'A5164', 'A5165', 'A5166', 'A5167', 'A5168', 'A5169', 'A5170', 'A5171', 'A5172', 'A5173', 'A5174', 'A5175', 'A5176', 'A5177', 'A5178', 'A5179', 'A5180', 'A5181', 'A5182', 'A5183', 'A5184', 'A5185', 'A5186', 'A5187', 'A5188', 'A5189', 'A5190', 'A5191', 'A5192', 'A5193', 'A5195', 'A5196', 'A5197', 'A5198', 'A5199', 'A5200', 'A5201', 'A5202', 'A5203', 'A5204', 'A5205', 'A5206', 'A5207', 'A5208', 'A5209', 'A5210', 'A5211', 'A5212', 'A5213', 'A5214', 'A5215', 'A5216', 'A5217', 'A5218', 'A5219', 'A5220', 'A5221', 'A5222', 'A5223', 'A5224', 'A5225', 'A5226'))]
			}
			n <- length(unk)
			if (n > 0) {
				if (n > 8) {
					unk <- c(unk[1:6], "...")
				}
				message(paste0("removing ", n, " unknown samples:\n", paste(unk, collapse=", ")))
			}
			i <- i[!is.na(i)]		
		}
		
		j <- which(inf$type[i] == "reference") + 3
		ibs_ref <- ibs[, j]
		ibs_fld <- ibs[, -j]

		write.csv(ibs_ref, file.path(outpath, paste0(oname, "_ibs-ref.csv")), row.names=FALSE)
		write.csv(ibs_fld, file.path(outpath, paste0(oname, "_ibs-fld.csv")), row.names=FALSE)
		write.csv(inf, file.path(outpath, paste0(oname, "_info.csv")), row.names=FALSE)
		writexl::write_xlsx(x, file.path(outpath, paste0(oname, ".xlsx")), format_headers=FALSE)		
	}
}


