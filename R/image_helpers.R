
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
	

split_lump <- function(mdst, lump=0.05, split=0.2) {

	oldnms <- colnames(mdst)

	mdst <- matchpoint::split_groups(mdst, split, TRUE)
	mdst <- matchpoint::lump_similar(mdst, lump) 
	newnms <- colnames(mdst)

	i <- !grepl("_#_", newnms)
	j <- (newnms[i] != oldnms[i]) & (!grepl("_.$", newnms[i]))
	newnms[i][j] <- paste0(newnms[i][j], "*")

	i <- grep("_#_", newnms)
	newnms[i] <- gsub("_.$", "", newnms[i])
	newnms[i] <- gsub("_._#_", "_#_", newnms[i])

	simplify <- function(names, pattern) {
		nms <- names[!grepl("_#_", names)] |> unique()
		p <- unique(grep(pattern, nms, value=TRUE))
		if (length(p) > 0) {
			for (i in 1:length(p)) {
				n <- grep(gsub(pattern, "", p[i]), nms, value=TRUE) |> unique() |> length() 
				if (n == 1) {
					names[names==p[i]] <- gsub(pattern, "", p[i])
				}
			}
		}
		names
	}
	newnms <- simplify(newnms, "_a$")
	newnms <- simplify(newnms, "_b$")
	
	newnms <- collapse_names(newnms)
	newnms <- gsub("_#_", " # ", newnms)

	list(old=oldnms, new=newnms, split=split, lump=lump)
}

split_lump_old <- function(mdst, maxlump=0.05, minsplit=0.2) {

	oldnms <- colnames(mdst)
	v <- matchpoint::self_dist(mdst)
	qth <- stats::quantile(v$value)
	#outh <- max(minsplit, qth[4] + 1.5 * (qth[4] - qth[2]))
	outh <- minsplit
	splitln <- Inf
	if (outh < qth[5]) {
		splitln <- outh
		mdst <- matchpoint::split_groups(mdst, outh, TRUE)
		v <- matchpoint::self_dist(mdst)
		qth <- stats::quantile(v$value)
#		outh <- max(minsplit, qth[4] + 1.5 * (qth[4] - qth[2]))
#		if ((outh < qth[5]) && (outh < splitln) && (outh > minsplit)) {
#			splitln <- outh
#			mdst <- matchpoint::split_groups(mdst, outh, TRUE)
#			v <- matchpoint::self_dist(mdst)
#			qth <- quantile(v$value)
#		}
	}

	mdst <- matchpoint::lump_similar(mdst, max(.01, qth[2])) 

	v <- matchpoint::self_dist(mdst)
	cut <- min(maxlump, max(0.01, stats::quantile(v$value)[3]))
	mdst <- matchpoint::lump_similar(mdst, cut) 

	v <- matchpoint::self_dist(mdst)
	cut <- min(maxlump, max(0.01, stats::quantile(v$value)[4]))
	mdst <- matchpoint::lump_similar(mdst, cut) 
	
	newnms <- colnames(mdst)

	i <- !grepl("_#_", newnms)
	j <- (newnms[i] != oldnms[i]) & (!grepl("_.$", newnms[i]))
	newnms[i][j] <- paste0(newnms[i][j], "*")

	i <- grep("_#_", newnms)
	newnms[i] <- gsub("_.$", "", newnms[i])
	newnms[i] <- gsub("_._#_", "_#_", newnms[i])

	simplify <- function(names, pattern) {
		nms <- names[!grepl("_#_", names)] |> unique()
		p <- unique(grep(pattern, nms, value=TRUE))
		if (length(p) > 0) {
			for (i in 1:length(p)) {
				n <- grep(gsub(pattern, "", p[i]), nms, value=TRUE) |> unique() |> length() 
				if (n == 1) {
					names[names==p[i]] <- gsub(pattern, "", p[i])
				}
			}
		}
		names
	}
	newnms <- simplify(newnms, "_a$")
	newnms <- simplify(newnms, "_b$")
	
	newnms <- collapse_names(newnms)
	newnms <- gsub("_#_", " # ", newnms)

	list(old=oldnms, new=newnms, split=splitln, lump=cut)
}



order_names <- function() {	

	d <- data.frame(
		name = c('DCob22-7521_moreOrders', 'DEra22-7523_1_moreOrders', 'DCas23-7954', 'DCpea23-7956_moreOrders', 'DMz23-7957_7228', 'DRi23-7955_moreOrders', 'DCas23-7816', 'DCob23-7823', 'DMz23-7824', 'DRi23-7825'),
		crop = c("bean", "teff", "cassava", "cowpea", "maize", "rice", "cassava", "bean", "maize", "rice"),
		iso3 = c("ETH", "ETH", "NGA", "NGA", "NGA", "NGA", "TZA", "TZA", "TZA", "TZA"),
		country = c("Ethiopia", "Ethiopia", "Nigeria", "Nigeria", "Nigeria", "Nigeria", "Tanzania", "Tanzania", "Tanzania", "Tanzania")
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
	out <- info["TargetID"]
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
	info <- info[,c("TargetID", "variety")]
	colnames(info)[2] <- "field_id"
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

	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	ffld <- ffld[substr(basename(ffld), 1, 2) != "~$"]
	if (length(ffld) == 0) {
		stop("no identifier file found. Check path?")
	}
	fld <- as.data.frame(readxl::read_excel(ffld, 2))
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	wantnms <- c("crop", "order.id", "dart.id", "tosci.id", "planting.barcode", "var.name.full", "reference", "plate.id", "well.row", "well.col", "plate_row", "plate_col", "well.id")
	newnms <-  c("crop", "order", "sample", "inventory", "inventory", "variety", "reference", "plate.id", "well.row", "well.col", "well.row", "well.col", "well.id")
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
	fldref <- bindr(fld, ref)
	fldref$crop <- tolower(fldref$crop)
	unique(fldref)	
}


copy_dart_files <- function(path, outpath, ordernr, overwrite=FALSE) {
	ff <- list.files(path, pattern=paste0("^Report_", ordernr, ".*.csv$"), ignore.case=TRUE, full.names=TRUE)
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
		print(crop); utils::flush.console()
		croppath <- file.path(path, crop)

		inf <- varinfo[varinfo$crop == crop, ]

		cropf <- list.files(croppath, pattern="SNP.csv$", full.names=TRUE, ignore.case=TRUE)
		cropf2 <- gsub("SNP.csv$", "SNP_2row.csv", cropf)
		stopifnot(length(cropf) == 1)
		countf <- gsub("SNP.csv$", "Counts.csv", cropf)

		x <- matchpoint::read_dart(cropf) 
		matchpoint:::copy_dart_files(croppath, outpath, x$order)

		foutc <- file.path(outpath, gsub("Report_", "", basename(countf)))
		fout1 <- file.path(outpath, gsub("Report_", "", basename(cropf)))
		fout2 <- gsub("SNP.csv$", "SNP_2row.csv", fout1)
		if (!file.exists(cropf2)) {
			message("creating 2 row file")
			y <- make_dart_2row(x)
			matchpoint:::write_dart(y, fout2)
		}

		cnts0 <- matchpoint::read_dart(countf) 
		cnts <- matchpoint:::dart_make_unique(cnts0)		

		if (!identical(cnts, cnts0)) {  # Ethiopia
			message("not identical cnts cnts0")
			i <- duplicated(cnts0$geno)
			if (any(i)) {
				message("removing duplicate count data")
				cnts0$geno <- cnts0$geno[!i, ]
				cnts0$snp  <- cnts0$snp[, c(1, 2, which(!i)+2)]
			}
			cnts <- matchpoint:::dart_make_unique(cnts0)		
			
			x2 <- matchpoint::read_dart(cropf2) 
			
			patrn <- " \\[Original species: Eragrostis tef\\]"
			x$geno$plate_id <- gsub(patrn, "", x$geno$plate_id)
			x2$geno$plate_id <- gsub(patrn, "", x2$geno$plate_id)
			cnts$geno$plate_id <- gsub(patrn, "", cnts$geno$plate_id)

			if (is.null(cnts$geno$TargetID)) {
				message("no TargetID")
				## no TargetID
				cnts$geno$TargetID <- cnts$geno$genotype
				colnames(cnts$snp)[-(1:2)] <- cnts$geno$genotype
			}

			mflds <- c("plate_barcode", "plate_row", "plate_col", "genotype")
			tom <- cnts$geno[, mflds]
			tom$new_genotype <- tom$genotype
			tom$genotype <- gsub("_D.$", "", tom$genotype)

			matchm <- function(d, y=1) {
				d$geno$id <- 1:nrow(d$geno)
				m <- merge(d$geno, tom, by=mflds, all=TRUE)
				m <- m[order(m$id), ]				
				stopifnot(all(m$genotype == gsub("_D.$", "", m$new_genotype)))
				m$genotype <- NULL
				d$geno$id <- NULL
				colnames(m) <- gsub("new_genotype", "genotype", colnames(m))
				d$geno <- m[, colnames(d$geno)]
				colnames(d$snp)[-y] <- d$geno$genotype
				d
			}
			
			x <- matchm(x)
			x2 <- matchm(x2, 1:2)
			
			inf <- inf[, c("crop", "sample", "inventory", "variety", "reference")]

			m <- merge(inf, tom, by.x="sample", by.y="genotype")
			m$sample <- m$new_genotype
			m$new_genotype <- NULL
			cns <- gsub("plate_barcode", "plate.id", colnames(m))
			cns <- gsub("plate_row", "well.row", cns)
			colnames(m) <- gsub("plate_col", "well.col", cns)
			inf <- m[, colnames(varinfo)]

			matchpoint:::write_dart(cnts, foutc)
			matchpoint:::write_dart(x, fout1)
			matchpoint:::write_dart(x2, fout2)
			
#			head(cnts0$geno)		
#			cnts$snp[1:4, 1:6]
#			head(cnts$geno)
#			x$snp[1:4, 1:6]
#			head(x$geno)			
		} 
		
		stopifnot(length(unique(colnames(x$snp))) == ncol(x$snp))
		stopifnot(length(unique(colnames(cnts$snp))) == ncol(cnts$snp))

		i <- match(x$geno$genotype, cnts$geno$genotype)
		x$geno$TargetID <- cnts$geno$TargetID[i]

		colnames(x$marker)[colnames(x$marker) == "AlleleID"] <- "MarkerName"
		colnames(x$snp)[colnames(x$snp) == "AlleleID"] <- "MarkerName"
		#colnames(x$snp) <- make_unique_ids(colnames(x$snp))
		
		x$marker$MarkerName <- toupper(x$marker$MarkerName)
		x$snp$MarkerName <- toupper(x$snp$MarkerName)

		#x <- matchpoint::make_dart_1row(x)
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

	
		if (!is.null(loc)) {
			n <- nrow(inf)
			inf <- merge(inf, loc, by="inventory", all.x=TRUE)
			stopifnot(nrow(inf) == n)
			inf <- inf[, c(names(varinfo), c("longitude", "latitude"))]
		}
		
	#	x$geno$geno_match <- gsub("_D.$", "", x$geno$genotype)
		
		if ("order" %in% names(inf)) {
			x$info <- merge(inf, x$geno, by.x=c("order", "sample"), by.y=c("order", "genotype"))
			if (nrow(unique(x$info[, c("order", "sample")])) != nrow(x$info)) {
				message("info-file order/sample numbers were not unique")
			}
		} else {
			x$info <- merge(inf, x$geno, by.x="sample", by.y="genotype")
			if (length(unique(x$info$sample)) != nrow(x$info)) {
				message("info-file sample numbers were not unique")
			}		
		}
	#	stopifnot(nrow(x$geno) == nrow(x$info))
	#	x$info$sample <- x$info$genotype
	#	x$geno$geno_match <- NULL
	#	x$info$genotype <- NULL
		
#		ibs <- merge(x$marker[, 1:3], x$snp, all=TRUE)
#		dartnms <- gsub("_D.$", "", colnames(x$snp)[-1])

		dartnms <- colnames(x$snp)[-1]
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

		# cowpea mixtures
		x$info <- x$info[x$info$variety != "MIXTURES", ]

		utils::write.csv(x$info, paste0(bname, "_variety-info.csv"), row.names=FALSE)
#		utils::write.csv(x$marker, paste0(bname, "_marker-info.csv"), row.names=FALSE)
		
		x$order <- NULL
		writexl::write_xlsx(x, paste0(bname, ".xlsx"), format_headers=FALSE)
	}
}


