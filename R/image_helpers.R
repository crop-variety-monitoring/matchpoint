
collapse_names <- function(x, sep="_#_") {
	
	lcPrefix <- function (x) {
		x <- as.character(x)
		nc <- nchar(x, type = "char")
		for (i in 1:min(nc)) {
			ss <- substr(x, 1, i)
			if (any(ss != ss[1])) {
				return(substr(x[1], 1, i - 1))
			}
		}
		substr(x[1], 1, i)
	}
	
	nms <- strsplit(x, sep)
	sapply(nms, \(n) { 
		if (length(n) > 1) {
			u <- strsplit(n, sub("\\D*(\\d+).*", "\\1", n))
			u <- sapply(u, \(i) i[1]) |> unique()
			u <- u[u!=""]
			if (length(u) > 0) {
				out <- ""[0]
				for (i in 1:length(u)) {
					j <- grep(u[i], n, fixed=TRUE)
					if (length(j) == 1) {
						out <- c(out, n[j])
					} else {
						out <- c(out, paste(c(n[j][1], gsub(u[i], "", n[j][-1])), collapse="*"))
					}
				}
				paste(out, collapse=sep)
			} else {
				paste(n, collapse=sep)
			}
		} else {
			paste(n, collapse=sep)
		}
	})		
}
	


order_names <- function() {	

	d <- data.frame(
		name = c('DCob22-7521_moreOrders', 'DEra22-7523_1_moreOrders', 'DMz23-8837_DMz23-8838', 'DCas23-7954', 'DCpea23-7956_moreOrders', 'DMz23-7957_7228', 'DRi23-7955_moreOrders', 'DCas23-7816', 'DCob23-7823', 'DMz23-7824', 'DRi23-7825'),
		crop = c("bean", "teff", "maize", "cassava", "cowpea", "maize", "rice", "cassava", "bean", "maize", "rice"),
		iso3 = c("ETH", "ETH", "ETH","NGA", "NGA", "NGA", "NGA", "TZA", "TZA", "TZA", "TZA"),
		country = c("Ethiopia", "Ethiopia", "Ethiopia", "Nigeria", "Nigeria", "Nigeria", "Nigeria", "Tanzania", "Tanzania", "Tanzania", "Tanzania")
	)
	d$cc <- paste(d$country, d$crop, sep=" / ")
	d
}
	


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


DAP_info <- function(filename) {
	info <- utils::read.csv(filename)
	out <- list()
	out$TargetID <- info["targetID"]
	out$Genotype <- make_unique_ids(info$sample)
	out$RefType <- gsub("\\*$| \\*$", "", info$variety)
	out$SampleType <- ifelse(info$reference, NA, "field")
	out
}	


DAP_write_excel <- function(x, info, filename) {

	nms <- names(x$res_full)
	resfull = do.call(rbind, lapply(1:length(x$res_full), 
				\(i) data.frame(field_TID=nms[i], x$res_full[[i]])))
	colnames(resfull)[1:4] <- c("field_Tid", "ref_Tid", "ref_id", "variety")
	info <- info[,c("targetID", "variety")]
	colnames(info)[1:2] <- c("targetID", "field_id")
	full <- merge(resfull, info, by=1, all.x=TRUE)
	full$var_rank <- with(full, stats::ave(Probability, field_id, FUN=\(x) rank(1000 - x, ties.method="min")))
	full$Probability <- full$Probability/100
	x[[2]] <- full
	colnames(x$res_summary)[1:7] <- c("field_Tid", "field_id", "ref_Tid", "ref_id", "variety", "NA.perc", "Probability")
	x$res_summary$Probability <- x$res_summary$Probability/100
	x[[3]] <- as.data.frame(x$ref_distance, check.names=FALSE)
	writexl::write_xlsx(x[1:3], paste0(filename, ".xlsx"))
}


bindr <- function( ...) {
	x <- list(...)
	nms <- unique(unlist(lapply(x, names)))
	x <- lapply(x, function(x) data.frame(c(x, sapply(setdiff(nms, names(x)), function(y) NA))))
	x$make.row.names <- FALSE
	do.call(rbind, x)
}


fix_varnames <- function(vars) {

	m <- matrix(byrow=TRUE, ncol=2, c(
# TZA
		"Kablent Progeny", "Kablanket Progeny",
		"CML 442(M WE 4106)", "CML442(MWE4106)",
		"CML442(MWE410)", "CML442(MWE4106)",
		"BORA", "Bora",
		"PASI", "Pasi",
		"PESA", "Pesa",
		"ROJO", "Rojo",
		"DK 90  89", "DK 9089",
		"DK 80 -31", "DK 80-31",
		"KALALU", "Kalalu",
		"KIPAPI", "Kipapi",
		"KIzimbani", "Kizimbani",
		"Kibanda Memo", "Kibanda meno",
		"Kibanda memo", "Kibanda meno",
		"Kibanda Meno", "Kibanda meno",
		"Nerika 4", "NERICA 4",
		"(OBA-98)", "OBA-98", 
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
		"SC645", "SC 645",
		"SC419", "SC 419",
		"SAMMAZ-11", "SAMMAZ 11",
		"WAC 42PVEE", "WAC42PVEE",
		"OBASUPER 1", "OBA SUPER 1",
		"OBASUPER 2", "OBA SUPER 2",
		"ART- 98-SW-1", "ART-98-SW-1",
		"ARICA 14", "ARICA-14",
		"ARICA 4", "ARICA-4",
		"IITA-TMS-IBA00070", "IITA TMS IBA00070",
		"NR-8083", "NR 8083"
	))

	h <- cbind(1:length(vars), match(vars, m[,1])) 
	h <- h[!is.na(h[,2]), ]
	if (nrow(h) > 0) {
#TZA
		vars[h[,1]] <- m[h[,2], 2]
		vars <- gsub("Nerica", "NERICA", vars)
		vars <- gsub("NERICA-", "NERICA ", vars)
		vars <- gsub("TARI CASS", "TARI CASS ", vars)
		vars <- gsub("TARICASS", "TARI CASS ", vars)
#NGA 
		vars <- gsub("SAMPEA-", "SAMPEA ", vars)
		"SAMMAZ-11"
		vars <- gsub("ARICA-", "SAMPEA ", vars)
		vars <- gsub(" \\(DROUGHTTEGO®WE52..\\)", "", vars)
	}	
	
	vars
}	



get_varinfo <- function(path, fix_vars=TRUE) {

	ffld <- list.files(path, pattern="identifier.xlsx$", full.names=TRUE)
	ffld <- ffld[substr(basename(ffld), 1, 2) != "~$"]
	if (length(ffld) == 0) {
		stop("no identifier file found. Check path?")
	}
	fld <- as.data.frame(readxl::read_excel(ffld, 2))
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	wantnms <- c("crop", "order.id", "dart.id", "tosci.id", "planting.barcode", "var.name.full", "reference", "plate.id", "well.row", "well.col", "plate_row", "plate_col", "well.id", "target.id")
	newnms <-  c("crop", "order", "sample", "inventory", "inventory", "variety", "reference", "plate.id", "well.row", "well.col", "well.row", "well.col", "well.id", "target.id")
	i <- which(wantnms %in% names(fld))
	fld <- fld[, wantnms[i]]
	names(fld) <- newnms[i]
	fld$variety <- ""
	fld$reference <- FALSE
	if (!is.null(fld$well.id)) {
		fld$well.row <- substr(fld$well.id, 1, 1)
		fld$well.col <-	substr(fld$well.id, 2, 3)
		fld$well.id <- NULL
	}
	
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	fref <- fref[substr(basename(fref), 1, 2) != "~$"]
	ref <- as.data.frame(readxl::read_excel(fref, 2))
	j <- which(wantnms %in% names(ref))
	ref <- ref[, wantnms[j]]
	names(ref) <- newnms[j]
	ref$reference <- TRUE
	if (fix_vars) {
		ref$variety <- fix_varnames(ref$variety)
		nuv <- length(unique(ref$variety))
		uv <- gsub(" ", "", (tolower(unique(ref$variety))))
		uv <- gsub("-", "", uv)
		stopifnot(length(unique(uv)) == nuv)
	}
	#b = table(gsub("-", "", gsub(" ", "", s)))
	#b[b>1]
	fldref <- bindr(ref, fld)
	fldref$crop <- tolower(fldref$crop)
	unique(fldref)	
}


copy_dart_files <- function(path, outpath, ordernr, overwrite=TRUE, counts_only=FALSE) {
	ff <- list.files(path, pattern=paste0("^Report_", ordernr, ".*.csv$"), ignore.case=TRUE, full.names=TRUE)
	if (counts_only) ff <- grep("_Counts.csv$", ff, value=TRUE)
	outff <- gsub("Report_", "", basename(ff))
	dir.create(outpath, FALSE, TRUE)
	file.copy(ff, file.path(outpath, outff), overwrite=overwrite)
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
		message(crop); utils::flush.console()
		croppath <- file.path(path, crop)

		inf <- varinfo[varinfo$crop == crop, ]

		cropf1 <- list.files(croppath, pattern="SNP.csv$", full.names=TRUE, include.dirs=TRUE, ignore.case=TRUE) |> normalizePath()
		cropf2 <- gsub("SNP.csv$", "SNP_2row.csv", cropf1)
		countf <- gsub("SNP.csv$", "Counts.csv", cropf1)

		if (length(cropf1) > 1) {  
			# first combine everything into one file
			message(paste0("   combining dart files"))
			xin <- lapply(cropf1, matchpoint::read_dart) 
			x <- matchpoint:::dart_combine(xin)
			cropf1 <- fout1 <- file.path(outpath, paste0(x$order, "_SNP.csv"))
			matchpoint:::write_dart(x, fout1)
			
			y <- matchpoint:::make_dart_2row(x)
			cropf2 <- fout2 <- file.path(outpath, paste0(x$order, "_SNP_2row.csv"))
			matchpoint:::write_dart(y, fout2)

			z <- lapply(countf, matchpoint::read_dart, useTID=FALSE) 
			cnts <- matchpoint:::dart_combine(z)
			stopifnot(x$order == cnts$order)
			foutc <- countf <- file.path(outpath, paste0(cnts$order, "_Counts.csv"))
			matchpoint:::write_dart(cnts, foutc)

		} else {
		
			foutc <- file.path(outpath, gsub("Report_", "", basename(countf))) |> normalizePath(mustWork=FALSE) 
			fout1 <- file.path(outpath, gsub("Report_", "", basename(cropf1))) |> normalizePath(mustWork=FALSE)
			fout2 <- gsub("SNP.csv$", "SNP_2row.csv", fout1) |> normalizePath(mustWork=FALSE)
			if (any(cropf1 == fout1)) {
				stop("same input as output")
			}
		}

		
		row1 <- matchpoint::read_dart(cropf1) 
		if (length(unique(row1$geno$ID)) < nrow(row1$geno)) {
#		if (length(unique(colnames(row1$snp))) < ncol(row1$snp)) {
			cnts <- matchpoint::read_dart(countf) 
			copycounts <- TRUE
			if (is.null(cnts$geno$targetID)) {
				cnts$geno$ID <- cnts$geno$targetID <- 1:nrow(cnts$geno)
				message("   IDs not unique, but Counts file has no targetID; assigning 1,2,3,...")
				copycounts <- FALSE
				matchpoint:::write_dart(cnts, foutc)
			}
			ids <- cnts$geno[, c("genotypeID", "targetID")]
			
			message("   transferring targetID")

			row2 <- matchpoint::read_dart(cropf2) 
			stopifnot(all(colnames(row2 == colnames(row1))))
			cns <- colnames(row1$snp)
			#mtch <- lapply(1:length(cns), function(i) cbind(i, which(ids$genotypeID %in% cns[i])) )
			#mtch <- lapply(mtch, function(m) cbind(m, ids[m[,2], ]))
			#mtch <- do.call(rbind, mtch)

			u <- unique(cns)
			out <- rep(NA, length(cns))
			gout <- rep(as.integer(NA), length(cns))
			for (i in u) {
				k <- which(cns==i)
				# 1:length(k) because teff has more genotypes in cnts and inf than in 1_row
				out[k] <- rep_len(cnts$geno$targetID[which(cnts$geno$genotypeID == i)], length(k))
				gout[k] <- rep_len(which(inf$sample==i), length(k))
			}
			
			if (any(is.na(out))) {
				stop("1_row file does not have unique IDs and they do not match targetID of Counts file")
			}
			colnames(row1$snp) <- row1$geno$targetID <- row1$geno$ID <- out
			colnames(row2$snp) <- row2$geno$targetID <- row2$geno$ID <- out
			inf <- inf[gout, ]
			inf$targetID <- row1$geno$targetID
			inf <- inf[!is.na(inf$targetID), ]
			
			# new refs with duplicate sample numbers may have been added
			inf$reference[inf$sample %in% inf$sample[inf$reference]] <- TRUE
						
			matchpoint:::write_dart(row1, fout1)
			matchpoint:::write_dart(row2, fout2)

			if (copycounts) matchpoint:::copy_dart_files(croppath, outpath, row1$order, counts_only=TRUE)
		} else {
			matchpoint:::copy_dart_files(croppath, outpath, row1$order)
			if (!file.exists(cropf2)) {
				message("   creating 2 row file")
				y <- make_dart_2row(row1)
				matchpoint:::write_dart(y, fout2)
			}
		}
		utils::flush.console()

		pos <- matchpoint::marker_positions(crop)
		i <- pos$MarkerName %in% row1$markers$marker
		pos <- pos[i, ]
		if (nrow(pos) != nrow(row1$markers)) {
			message(paste("   SNP data has", nrow(row1$markers), "markers. Panel has", nrow(pos), "matches")) 
		}
		m <- match(row1$markers$MarkerName, pos$MarkerName)	
		if (any(is.na(m))) message(paste0(sum(is.na(m)), " unknown markers found"))
		row1$markers$order <- 1:nrow(row1$markers)
		m <- merge(pos, row1$markers, by="MarkerName", all.y=TRUE)
		m <- m[m$order, ]
		m$order <- NULL
		row1$markers <- m

		if (!is.null(loc)) {
			n <- nrow(inf)
			inf <- merge(inf, loc, by="inventory", all.x=TRUE)
			stopifnot(nrow(inf) == n)
			inf <- inf[, c(names(varinfo), c("longitude", "latitude"))]
		}
		
		# cowpea mixtures
		inf <- inf[inf$variety != "MIXTURES", ]
		
		if ("order" %in% names(inf)) {
			test <- merge(inf, row1$geno, by.x=c("order", "sample"), by.y=c("order", "ID"))
			if (nrow(unique(test[, c("order", "sample")])) != nrow(test)) {
				message("   info-file order/sample numbers are not unique")
			}
		} else {
			test <- merge(inf, row1$geno, by.x="sample", by.y="ID")
			if (length(unique(test$sample)) != nrow(test)) {
				message("   info-file sample numbers are not unique")
			}
		}
		#dartnms <- row1$geno$genotype
		#refnms <- as.character(inf[inf$crop == crop, "sample"])
		#i <- match(dartnms, refnms)

		bname <- file.path(outpath, row1$order)


		utils::write.csv(inf, paste0(bname, "_variety-info.csv"), row.names=FALSE)
#		utils::write.csv(row1$marker, paste0(bname, "_marker-info.csv"), row.names=FALSE)
		
		row1$order <- row1$type <- NULL
		row1$snp <- data.frame(SNP=rownames(row1$snp), row1$snp)
		writexl::write_xlsx(row1, paste0(bname, ".xlsx"), format_headers=FALSE)
	}
}


