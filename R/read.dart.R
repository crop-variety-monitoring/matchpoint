
read_dart <- function(f) {
	r <- read.csv(f, header=FALSE)
	x <- r[-1, ]
	x[x != "*"] <- NA
	srow <- which(apply(x, 1, \(i) all(is.na(i))))[1] + 1
	scol <- which(apply(x, 2, \(i) all(is.na(i))))[1]

	marker <- r[(srow+1):nrow(r), 1:(scol-1), ]
	colnames(marker) <- r[srow, 1:(scol-1)]

	hdr <- data.frame(t(r[1:srow, scol:ncol(r)]))

	d <- r[(srow+1):nrow(r), scol:ncol(r), ]
	d[d == "-"] <- NA
	
	d <- data.frame(r[(srow+1):nrow(r), 1:2], lapply(d, as.integer))
	if (isTRUE(any(d[,-c(1:2)] > 2))) { 
		ids <- paste(trimws(hdr[,1]), trimws(hdr[, ncol(hdr)-1]), sep="_")
		colnames(d) <- c(colnames(marker)[1:2], ids)
		return(list(snp=d, marker=marker, geno=hdr, type="counts"))
	} 
	i <- seq(1, nrow(d), 2)
	if (all(d[i,1] == d[i+1,1])) {
		ids <- paste(trimws(hdr[,1]), trimws(hdr[,ncol(hdr)]), sep="_")
		colnames(d) <- c(colnames(marker)[1:2], ids)
		list(snp=d, marker=marker, geno=hdr, type="2_row")
	} else if (nrow(d) == length(unique(d[, 1]))) {
		d <- d[,-2]
		ids <- paste(trimws(hdr[,1]), trimws(hdr[,ncol(hdr)]), sep="_")
		colnames(d) <- c(colnames(marker)[1], ids)
		list(snp=d, marker=marker, geno=hdr, type="1_row")	
	} else {
		warning("don't know if this is a Count, 1_row or a 2_row file")
		list(snp=d, marker=marker, geno=hdr, type="???")		
	}
}


make_dart_1row <- function(x) {
	if (x$type=="1_row") return(x)
	stopifnot(x$type=="2_row")
#	||(x$type=="counts"))
	d <- x$snp
	i <- seq(1, nrow(d), 2)
	dd <- as.matrix(d[,-c(1:2)])

	if (x$type == "counts") {
		# not OK as we also need marker specific thresholds for pres/abs
		threshold = 0.2
		f1 <- dd[i,] / dd[i+1,]
		f1[] <- pmin(1, f1)
		f1[dd[i,] == 0] <- 0
		f1[dd[i+1,] == 0 & dd[i,] != 0] <- 1
		f1 <- ifelse(f1 < threshold, 0, 1)

		f2 <- dd[i+1,] / dd[i,]
		f2[] <- pmin(1, f2)
		f2[dd[i+1,] == 0] <- 0
		f2[dd[i,] == 0 & dd[i+1,] != 0] <- 1
		f2 <- ifelse(f2 < threshold, 0, 1)
		snp = 10 * f1 + f2
	} else {
		snp = 10 * dd[i,] + dd[i+1,]
	}
	
	snp[snp== 0] <- NA
	snp[snp==10] <- 0
	snp[snp==11] <- 2
	
	x$snp <- cbind(d[i, 1, drop=FALSE], snp)

	m <- x$marker[i+1,]
	m <- cbind (m[1], x$marker[i,2], m[-1])
	colnames(m)[2:3] <- paste0(colnames(x$marker)[2], c("Ref", "Alt"))
	x$marker <- m
	x$type="1_row"
	x
}


make_dart_2row <- function(x) {
	if (x$type=="2_row") return(x)
	stopifnot(x$type=="1_row")
	d <- x$snp
	e <- x$marker 
	first  <- c(0,1,0,1)
	second <- c(0,0,1,1)
	dd <- as.matrix(d[,-1])
	a = lapply(1:nrow(d), \(i) {
			v <- dd[i, ] + 2
			v[is.na(v)] <- 1
			rbind(first[v], second[v])
		})
	a <- do.call(rbind, a)
	a <- data.frame(rep(d[,1], each=2), as.vector(t(e[,2:3])), a)
	names(a) <- c("MarkerName", "AlleleSequence", names(d)[-1])

}



read_dart_folder <- function(path) {
#	ff <- list.files(path, pattern="^metadata_.*xlsx$", full=TRUE)
#	orders <- gsub("^metadata_", "", basename(ff))
#	orders <- gsub("\\.xlsx$", "", orders)
#	marker <- lapply(ff, \(f) suppressMessages(as.data.frame(readxl::read_excel(f))))

#	ff <- list.files(path, pattern="Counts.csv$", full.names=TRUE, ignore.case=TRUE)
#	orders <- gsub("^Report_", "", basename(ff))
#	orders <- gsub("\\_Counts.csv$", "", orders)

	ff <- list.files(path, pattern="SNP.csv$", full.names=TRUE, ignore.case=TRUE)
	ff2 <- list.files(path, pattern="2row.csv$", full.names=TRUE, ignore.case=TRUE)
	if (length(ff2) > length(ff)) {
		ff <- ff2
	}	
	if (length(ff) == 0) {
		stop("no files")
	}
	
	if (length(ff) == 1) {
		out <- read_dart(ff) 
	} else {
		x <- lapply(ff, read_dart)
		n <- sapply(x, \(i) nrow(i$snp))
		if (length(unique(n)) != 1) {
			warning(paste("number of markers does not match between files:", unique(n)))
		}
		snp <- x[[1]]$snp
		snp$sortID <- 1:nrow(snp)
		hdr <- x[[1]]$geno
		for (i in 2:length(x)) {
			snp <- merge(snp, x[[i]]$snp, by=1:2, all=TRUE)
			hdr <- rbind(hdr, x[[i]]$geno)
		}
		out <- x[[1]]
		snp <- snp[order(snp$sortID), ]
		snp$sortID <- NULL
		out$snp <- snp
		out$geno <- hdr
	}
	out
}


prepare_dart <- function(path) {
	outpath <- gsub("raw/dart", "clean", path)
	dir.create(outpath, FALSE, FALSE)
	crops <- list.dirs(path, full.names=FALSE)
	crops <- crops[crops != ""]
	crops <- crops[!grepl("/", crops)]
	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	fld <- as.data.frame(readxl::read_excel(ffld, 2))[, c("crop", "dart.id", "planting.barcode")]
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	colnames(fld)[3] <- "field.id"
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	ref <- as.data.frame(readxl::read_excel(fref, 2))[, c("crop", "dart.id", "var.name.full")]
	ref$crop <- tolower(ref$crop)
	intpath <- gsub("raw/dart", "intermediate", path)
	floc <- list.files(intpath, pattern="location", full.names=TRUE)
	loc <- NULL
	if (length(floc) == 1) {
		loc <- readxl::read_excel(floc)
	}
	fpan <- file.path(path, "DarTAG Mid density panel info.xlsx")
	
	for (crop in crops) {
		print(crop); flush.console()
		x <- try(matchpoint::read_dart_folder(file.path(path, crop)))
		if (inherits(x, "try-error")) {
			next
		}
		x <- matchpoint::make_dart_1row(x)
		x$ref_id <- ref[ref$crop == crop, ]
		x$field_id <- fld[fld$crop == crop, ]
		x$type <- NULL
		if (crop == "cassava") {
			x$position <- data.frame(do.call(rbind, strsplit(x$snp[,1], "_")))
			colnames(x$position) <- c("dart.id", "chr", "pos")
		} else {
			x$position <- as.data.frame(readxl::read_excel(fpan, crop))[, c("DART.ID", "Chr", "Pos")]
			colnames(x$position) <- tolower(colnames(x$position))
		}
		x$geo <- loc
		writexl::write_xlsx(x, file.path(outpath, paste0(crop, ".xlsx")), format_headers=FALSE)
	}
}




read_dart_order <- function(path, ordnr) {
	f <- list.files(path, pattern=paste0("^Report_", ordnr), full=TRUE)
	d <- lapply(f, read_dart)
	nms <- sapply(d, \(i) i$type)
	j1 <- grepl("1_row", nms)
	j2 <- grepl("2_row", nms)
	j3 <- grepl("counts", nms)
	if (any(j1)) {
		x <- d[[which(j1)]]
		names(x)[1] <- "row1"
		if (any(j2)) {
			x$row2 <- d[[which(j2)]]$snp
		} else {
			# make 2row
		}
	} else if (any(j2)) {
		row2 <- d[[which(j2)]]
		x <- make_dart_1row(row2)
		names(x)[1] <- "row1"
		x$row2 <- d[[which(j2)]]$snp
	} else {
		stop("need 1 or 2row data")
	}
	if (any(j3)) {
		x$counts <- d[[which(j3)]]$snp
	}
	x
}



prepare_dart2 <- function(path) {
	outpath <- gsub("raw/dart", "clean", path)
	dir.create(outpath, FALSE, FALSE)
	crops <- list.dirs(path, full.names=FALSE)
	crops <- crops[crops != ""]
	crops <- crops[!grepl("/", crops)]
	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	fld <- as.data.frame(readxl::read_excel(ffld, 2))[, c("crop", "dart.id", "planting.barcode")]
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	colnames(fld)[3] <- "field.id"
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	ref <- as.data.frame(readxl::read_excel(fref, 2))[, c("crop", "dart.id", "var.name.full")]
	ref$crop <- tolower(ref$crop)
	intpath <- gsub("raw/dart", "intermediate", path)
	floc <- list.files(intpath, pattern="location", full.names=TRUE)
	loc <- NULL
	if (length(floc) == 1) {
		loc <- readxl::read_excel(floc)
	}
	fpan <- file.path(path, "DarTAG Mid density panel info.xlsx")
	
	for (crop in crops) {
		print(crop); flush.console()
		croppath <- file.path(path, crop)
		ff <- list.files(croppath, pattern="^Report", full=TRUE)
		orders <- do.call(rbind, lapply(strsplit(ff, "_"), \(i) i[1:2]))[,2] |> unique()
		x <- lapply(orders, \(i) read_dart_order(croppath, i))
		names(x) <- orders
		
		for (i in 1:length(x)) {	
			x[[i]]$ref_id <- ref[ref$crop == crop, ]
			x[[i]]$field_id <- fld[fld$crop == crop, ]
			x[[i]]$type <- NULL
			if (crop == "cassava") {
				x[[i]]$position <- data.frame(do.call(rbind, strsplit(x[[i]]$row1[,1], "_")))
				colnames(x[[i]]$position) <- c("dart.id", "chr", "pos")
			} else {
				x[[i]]$position <- as.data.frame(readxl::read_excel(fpan, crop))[, c("DART.ID", "Chr", "Pos")]
				colnames(x[[i]]$position) <- tolower(colnames(x[[i]]$position))
			}
			x[[i]]$geo <- loc
			x[[i]]$type <- NULL
			nms <- c("counts", "row2", "row1", "geno", "marker", "position", "geo")
			j <- match(nms, names(x[[i]]))
			x[[i]] <- x[[i]][na.omit(j)]
			writexl::write_xlsx(x[[i]], file.path(outpath, paste0(crop, "_", orders[i], ".xlsx")), format_headers=FALSE)
		}
	}
}


make_infofiles <- function(path) {
	outpath <- gsub("raw/dart", "clean", path)
	dir.create(outpath, FALSE, FALSE)
	ffld <- list.files(path, pattern="identifier", full.names=TRUE)
	fld <- as.data.frame(readxl::read_excel(ffld, 2))[, c("crop", "dart.id", "planting.barcode")]
	fld$planting.barcode <- gsub(".jpg", "", fld$planting.barcode)
	colnames(fld)[3] <- "field.id"
	fld$variety <- ""
	fld$type <- "survey"
	fref <- list.files(path, pattern="inventory", full.names=TRUE)
	ref <- as.data.frame(readxl::read_excel(fref, 2))[, c("crop", "dart.id", "var.name.full")]
	colnames(ref)[3] <- "variety"
	ref$field.id <- NA
	ref$type <- "reference"
	z <- rbind(fld, ref)
	
	z$crop <- tolower(z$crop)
	
	for (crop in unique(z$crop)) {
		zc <- z[z$crop == crop, ]
		write.csv(zc, file.path(outpath, paste0(crop, "_info.csv")))
	}
}
