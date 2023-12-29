
make_unique_ids <- function(x) {
	tab <- table(x)
	dups <- tab[tab > 1]
	if (length(dups) > 0) {
		message("adding _Dn to duplicates")
		for (i in 1:length(dups)) {
			n <- dups[i]
			nm <- names(n)
			x[x == nm] <- paste0(nm, "_D", 1:n)
		}
	}
	x
}


assign_info <- function(filename) {
	info <- read.csv(filename)
	out <- info["TargetID"]
	out$Genotype <- make_unique_ids(info$sample)
	out$RefType <- info$variety
	out$SampleType <- ifelse(info$reference, NA, "field")
	out
}	


assign_write_excel <- function(x, info, filename) {
	nms <- names(x$res_full)
	resfull = do.call(rbind, lapply(1:length(x$res_full), 
				\(i) data.frame(field_TID=nms[i], x$res_full[[i]])))
	colnames(resfull)[1:3] <- c("field_Tid", "ref_Tid", "ref_id")
	info <- info[,c("TargetID", "Genotype")]
	colnames(info)[2] <- "field_id"
	full <- merge(resfull, info, by=1, all.x=TRUE)
	full$rank <- with(full, ave(Probability, field_id, FUN=\(x) rank(110 - x)))
	x[[2]] <- full
	colnames(x$res_summary)[1:4] <- c("field_Tid", "field_id", "ref_Tid", "ref_id variety")
	writexl::write_xlsx(x[1:2], paste0(filename, ".xlsx"))
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
	fld <- as.data.frame(readxl::read_excel(ffld, 2))
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	if ("order.id" %in% names(fld)) {
		fld <- fld[, c("crop", "order.id", "dart.id", "planting.barcode")]
		colnames(fld)[2:4] <- c("order", "sample", "inventory")
	} else {
		fld <- fld[, c("crop", "dart.id", "planting.barcode")]
		colnames(fld)[2:3] <- c("sample", "inventory")	
	}
	fld$variety <- ""
	fld$reference <- FALSE
	
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	fref <- fref[substr(basename(fref), 1, 2) != "~$"]
	ref <- as.data.frame(readxl::read_excel(fref, 2))
	if ("order.id" %in% names(ref)) {
		ref <- ref[, c("crop", "order.id", "dart.id", "var.name.full")]
		colnames(ref)[2:4] <- c("order", "sample", "variety")
	} else {
		ref <- ref[, c("crop", "dart.id", "var.name.full")]
		colnames(ref)[2:3] <- c("sample", "variety")	
	}
	ref$inventory <- ""
	ref$reference <- TRUE
	
	m <- matrix(byrow=TRUE, ncol=2, c(
# TZA
		"BORA", "Bora",
		"PASI", "Pasi",
		"PESA", "Pesa",
		"ROJO", "Rojo",
		"DK 90  89", "DK 9089",
		"DK 80 -31", "DK 80-31",
		"KALALU", "Kalalu",
		"KIPAPI", "Kipapi",
		"KIzimbani", "Kizimbani",
		"Nerika 4", "NERICA 4",
		'SUPA', 'Supa',  
		'TAI', 'Tai',
		"Meru HB 623", "MERU HB 623",
		"NATA -K 60 2023", "NATA K 60 2023",
		"NATAH 104", "NATA H 104",
		"PAN 4m 23", "PAN 4M-23",
		"SC419", "SC 419",
		"SC513", "SC 513",
		"SC627", "SC 627",   
		"SC719", "SC 719",
		"Selian  97", "Selian 97",
		"SITUKA", "Situka",  # M1?
		"SITUKA M1", "Situka M1",
		"Sy 514", "SY 514",
		"T105", "T 105",
		"TMV1", "TMV 1",
		"TMV2", "TMV 2",
		# "TMV 2", "TMV 2(FUH 6105)", ??
		"UYOLE 84", "Uyole 84",  		
		"WE2109", "WE 2109",
		"WE2113", "WE 2113",
		"WE3102", "WE 3102",
		"WE5135", "WE 5135",
		"WE5141", "WE 5141",
		"WE7118", "WE 7118",
		"WE7133", "WE 7133",
		"UH 5350", "UHS 5350", 
		"UH 6303", "UHS 6303",
		"UH5350", "UHS 5350",
		"UH615",  "UHS 615",
		"UH6303", "UHS 6303",
		"UHS 5210", "UHS 5210 (F UH 6303)", 
		"UHS5350", "UHS 5350",
		"UL49(FUH5350)", "UL 49 (FUH5350)", 
		"UL49/UL5095(FUH5350)", "UL 49/UL 5095 (FUH5350)",
		"UL5095(F6303)", "UL 5095(F6303)",
		"UL5218", "UL 5218",
#NGA 
		"SC 645", "SC645"
	))

	h <- cbind(1:nrow(ref), match(ref$variety, m[,1])) 
	h <- h[!is.na(h[,2]), ]
	if (nrow(h) > 0) {
#TZA
		ref$variety[h[,1]] <- m[h[,2], 2]
		ref$variety <- gsub("Nerica", "NERICA", ref$variety)
		ref$variety <- gsub("TARI CASS", "TARI CASS ", ref$variety)
		ref$variety <- gsub("TARICASS", "TARI CASS ", ref$variety)
#NGA 
		ref$variety <- gsub("SAMPEA-", "SAMPEA ", ref$variety)
	}	
	#s = sort(tolower(unique(ref$variety)))
	#b = table(gsub("-", "", gsub(" ", "", s)))
	#b[b>1]

	fldref <- rbind(fld, ref)
	fldref$crop <- tolower(fldref$crop)
	unique(fldref)	
}


copy_dart_files <- function(path, outpath, ordernr) {
	ff <- list.files(path, pattern=paste0("^Report_", ordernr, ".*.csv$"), ignore.case=TRUE, full.names=TRUE)
	outff <- gsub("Report_", "", basename(ff))
	dir.create(outpath, FALSE, TRUE)
	file.copy(ff, file.path(outpath, outff))
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
		matchpoint:::copy_dart_files(croppath, outpath, x$order)

		cnts0 <- matchpoint::read_dart(countf) 
		cnts <- matchpoint:::dart_make_unique(cnts0)
		if (!identical(cnts, cnts0)) {
			fout <- gsub("Report_", "", basename(countf))
			matchpoint:::write_dart(cnts, file.path(outpath, fout))
			x <- matchpoint:::dart_make_unique(x)
			fout <- gsub("Report_", "", basename(cropf))
			matchpoint:::write_dart(x, file.path(outpath, fout))
		}

		if (is.null(cnts$geno$TargetID)) {  # ETH teff
			x$geno$TargetID <- NA
		} else {
			i <- match(cnts$geno$genotype, x$geno$genotype)
			x$geno$TargetID <- cnts$geno$TargetID[i]
		}
		colnames(x$marker)[colnames(x$marker) == "AlleleID"] <- "MarkerName"
		colnames(x$snp)[colnames(x$snp) == "AlleleID"] <- "MarkerName"
		#colnames(x$snp) <- matchpoint:::make_unique_ids(colnames(x$snp))
		
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
		
		x$geno$geno_match <- gsub("_D.$", "", x$geno$genotype)
		
		if ("order" %in% names(inf)) {
			x$info <- merge(inf, x$geno, by.x=c("order", "sample"), by.y=c("order", "geno_match"))
			if (nrow(unique(x$info[, c("order", "sample")])) != nrow(x$info)) {
				message("info-file order/sample numbers were not unique")
			}
		} else {
			x$info <- merge(inf, x$geno, by.x="sample", by.y="geno_match")
			if (length(unique(x$info$sample)) != nrow(x$info)) {
				message("info-file sample numbers were not unique")
			}		
		}
		x$info$sample <- x$info$genotype
		x$info$genotype <- NULL
		x$geno$geno_match <- NULL
		
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
		
		x$order <- NULL
		writexl::write_xlsx(x, paste0(bname, ".xlsx"), format_headers=FALSE)
	}
}


