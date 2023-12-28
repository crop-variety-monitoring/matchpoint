
assign_info <- function(filename) {
	info <- read.csv(filename)
	out <- info["TargetID"]
	out$Genotype <- info$sample
	out$RefType <- info$variety
	out$SampleType <- ifelse(info$reference, NA, "field")
	out
}	


fake_it <- function(dart_file, outdir=".") {

#  dart_file <- "input/NGA/DCas23-7954_SNP.csv"
#  outdir <- "C:/github/cropvarmon/matchpoint/inst/ex"
	
	dart_file2 <- gsub(".csv$", "_2row.csv", dart_file)
	r1 <- data.frame(data.table::fread(dart_file, header=FALSE), check.names=FALSE)
	r2 <- data.frame(data.table::fread(dart_file2, header=FALSE), check.names=FALSE)

	srow <- sum(r1[1:30, 1] == "*") + 1
	scol <- sum(r1[1, 1:50] == "*") + 1

	take_anon_sample <- function(x, sid, nc=16) {
		y <- x[(srow+1):nrow(x), 1:nc]
		y <- y[y[,1] %in% sid, ]
		x <- x[1:(nrow(y)+srow), ]
		x[(srow+1):nrow(x), 1:nc] <- y
		x
	}
	
	s <- sort(sample(srow:nrow(r1), (nrow(r1)-srow-1)/2))
	sid <- r1[s,1]
	r1 <- take_anon_sample(r1, sid, scol-1)
	r2 <- take_anon_sample(r2, sid, scol)

	geno_file <- gsub("SNP.csv$", "variety-info.csv", dart_file)
	geno <- data.frame(data.table::fread(geno_file), check.names=FALSE)
	geno <- geno[sample(nrow(geno)), ]
	k <- match(geno$sample, r1[7,])
	geno <- geno[!is.na(k), ]

	i <- geno$variety == ""
	g <- geno[!i, ]
	g$variety <- paste0("var_", 1:nrow(g))
	g <- g[gtools::mixedorder(g$variety), ]

	gg <- geno[i, ]
	gg$longitude <- round(gg$longitude + stats::runif(nrow(gg), -0.5, 0.5), 2) 
	gg$latitude <- round(gg$latitude + stats::runif(nrow(gg), -0.5, 0.5), 2) 

	s <- sample(nrow(gg)/2)
	ss <- grep("undetermined", gg$inventory)
	s <- unique(c(s, ss))

	remid <- gg$sample[s]
	smpls <- unlist(r1[srow,])
	j <- smpls %in% remid

	r1s <- r1[, !j]
	r2s <- r2[, !j[-1]]
	gs <- gg[-s,]
	
	gs$inventory <- as.integer(as.factor(gs$inventory))

	gen <- rbind(g, gs)

	utils::write.csv(gen, file.path(outdir, "DCas00-0000_variety-info.csv"), row.names=FALSE, na="")

	utils::write.table(r1s, file.path(outdir, "DCas00-0000_SNP.csv"), row.names=FALSE, col.names=FALSE, na="", sep=",")
	utils::write.table(r2s, file.path(outdir, "DCas00-0000_SNP_2row.csv"), row.names=FALSE, col.names=FALSE, na="", sep=",")
}




get_varinfo <- function(path) {

	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	ffld <- ffld[substr(basename(ffld), 1, 2) != "~$"]
	if (length(ffld) == 0) {
		stop("no identifier file found. Check path?")
	}
	fld <- as.data.frame(readxl::read_excel(ffld, 2))[, c("crop", "dart.id", "planting.barcode")]
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	colnames(fld)[2:3] <- c("sample", "inventory")
	fld$variety <- ""
	fld$reference <- FALSE
	
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	fref <- fref[substr(basename(fref), 1, 2) != "~$"]
	ref <- as.data.frame(readxl::read_excel(fref, 2))[, c("crop", "dart.id", "var.name.full")]
	colnames(ref)[2:3] <- c("sample", "variety")
	ref$inventory <- ""
	ref$reference <- TRUE
	
	fldref <- rbind(fld, ref)
	fldref$crop <- tolower(fldref$crop)
	unique(fldref)	
}


prepare_dart <- function(path, outpath) {

	intpath <- gsub("raw/dart", "intermediate", path)
	dir.create(outpath, FALSE, TRUE)
	
	varinfo <- matchpoint:::get_varinfo(path)

	crops <- list.dirs(path, full.names=FALSE)
	crops <- crops[crops != ""]
	crops <- crops[!grepl("/", crops)]

	floc <- list.files(intpath, pattern="location", full.names=TRUE)
	loc <- NULL
	if (length(floc) == 1) {
		loc <- as.data.frame(readxl::read_excel(floc))
	}

	for (crop in crops) {
		print(crop); utils::flush.console()
		croppath <- file.path(path, crop)
		cropf <- list.files(croppath, pattern="SNP.csv$", full.names=TRUE, ignore.case=TRUE)
		stopifnot(length(cropf) == 1)
		countf <- list.files(croppath, pattern="Counts.csv$", full.names=TRUE, ignore.case=TRUE)
		stopifnot(length(countf) == 1)

		x <- matchpoint::read_dart(cropf) 
		cnts <- matchpoint::read_dart(countf) 
		i <- match(cnts$geno$genotype, x$geno$genotype)
		x$geno$TargetID <- cnts$geno$TargetID[i]
		
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

		x <- matchpoint::make_dart_1row(x)
		x$type <- NULL

		pos <- matchpoint::marker_positions(crop)
		i <- pos$MarkerName %in% x$marker$MarkerName
		#j <- !(x$marker$MarkerName %in% pos$MarkerName )
		pos <- pos[i, ]
		if (nrow(pos) != nrow(x$marker)) {
			message(paste("SNP data has", nrow(x$marker), "markers. Panel has", nrow(pos), "matches")) 
		}
		m <- match(x$marker$MarkerName, pos$MarkerName)
		if (any(is.na(m))) message(paste0(sum(is.na(m)), " unknown markers found"))
		x$marker$order <- 1:nrow(x$marker)
		m <- merge(pos, x$marker, by="MarkerName", all.y=TRUE)
		m <- m[m$order, ]
		m$order <- NULL
		x$marker <- m

		inf <- varinfo[varinfo$crop == crop, ]
		if (!is.null(loc)) {
			n <- nrow(inf)
			inf <- merge(inf, loc, by="inventory", all.x=TRUE)
			stopifnot(nrow(inf) == n)
			inf <- inf[, c(names(varinfo), c("longitude", "latitude"))]
		}

		#i <- match(inf$sample, cnts$geno$genotype)
		#inf$TargetID <- cnts$geno$TargetID[i]
		x$info <- merge(inf, x$geno, by.x="sample", by.y="genotype", all.x=TRUE)

#		ibs <- merge(x$marker[, 1:3], x$snp, all=TRUE)
		dartnms <- gsub("_D.$", "", colnames(x$snp)[-1])
		refnms <- as.character(inf[inf$crop == crop, "sample"])
		i <- match(dartnms, refnms)

#		write.csv(ibs[, c(1:3, which(!is.na(i))+3)], file.path(outpath, paste0(oname, "_IBS.csv")), na="", row.names=FALSE)

		unk <- dartnms[which(is.na(i))]

		if (length(unk) > 0) {
			if (x$order == "DMz23-7957_7228") {
				## undeclared 
				unk <- unk[!(unk %in% c('327', '328', '329', '330', '331', '332', '333', '334', '335', '336', '337', '338', '339', '340', '341', '342', '343', '344'))]
			} else if (x$order ==  "DCob23-7823") {
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
	

#		j <- which(inf$type[i] == "reference") + 3
#		ibs_ref <- ibs[, c(1:3, j)]
#		ibs_fld <- ibs[, -j]
#		utils::write.csv(ibs_ref, file.path(outpath, paste0(oname, "_IBS-ref.csv")), na="-", row.names=FALSE)
#		utils::write.csv(ibs_fld, file.path(outpath, paste0(oname, "_IBS-fld.csv")), na="-", row.names=FALSE)

		bname <- file.path(outpath, x$order)

		utils::write.csv(x$info, paste0(bname, "_variety-info.csv"), row.names=FALSE)
#		utils::write.csv(x$marker, paste0(bname, "_marker-info.csv"), row.names=FALSE)

		cropff <- list.files(croppath, pattern=paste0("^Report_", x$order, ".*.csv$"), ignore.case=TRUE, full.names=TRUE)
		outff <- gsub("Report_", "", basename(cropff))
		ok <- file.copy(cropff, file.path(outpath, outff))
		
		x$order <- NULL
		writexl::write_xlsx(x, paste0(bname, ".xlsx"), format_headers=FALSE)
	}
}


