
read_dart <- function(filename, useTID=TRUE) {

	r <- data.frame(data.table::fread(filename, header=FALSE), check.names=FALSE)
	ordname <- gsub("^Report_", "", basename(filename))
	ordname <- gsub("\\_SNP_2row.csv$|\\_SNP.csv$|\\_Counts.csv$", "", ordname)

	srow <- sum(r[1:30, 1] == "*") + 1
	scol <- sum(r[1, 1:50] == "*") + 1

	markers <- r[(srow+1):nrow(r), 1:(scol-1), ]
	colnames(markers) <- r[srow, 1:(scol-1)]

	hdr <- data.frame(t(r[1:srow, scol:ncol(r)]))
	cns <- c("order", "plate_barcode", "sample_reproducibility", "plate_row", "plate_col", "plate_id", "genotypeID")
	if (ncol(hdr) == 7) {
		colnames(hdr) <- cns
		hdr$ID <- hdr$genotypeID
	} else 	if (ncol(hdr) == 8) {
		colnames(hdr) <- c(cns, "targetID")
		if (useTID & (length(unique(hdr$genotypeID)) < nrow(hdr))) {
			hdr$ID <- hdr$targetID
		} else {
			hdr$ID <- hdr$genotypeID	
		}
	} else {
		stop("I do not understand the header of this file")
	}
	d <- as.matrix(r[(srow+1):nrow(r), scol:ncol(r)])
	v <- as.vector(d)
	v[v == "-"] <- NA
	d <- matrix(as.integer(v), nrow=nrow(d), ncol=ncol(d))	
	colnames(d) <- hdr$ID

#	d <- data.frame(r[(srow+1):nrow(r), 1:2], d, check.names=FALSE)
	mname <- function(m) {
		names(m) <- tolower(names(m))
		if (!is.null(m$variant)) {
			i <- m$variant != ""
			m[i,1] <- paste0(m[i,1], "|", m$variant[i])
			m[,1]
		} else {
			paste0(m[,1], "_", 1:2)
		}
	}

	colnames(markers)[1] <- "MarkerName"
	markers$MarkerName <- toupper(markers$MarkerName)

	if (grepl("_Counts", filename) || isTRUE(any(d > 2))) {
		marker <- mname(markers)
		ftype <- "counts"
	} else {
		i <- seq(2, nrow(d), 2)
		if (grepl("_2row", filename) || (all(markers[i,1] == markers[i-1,1]))) {
			marker <- mname(markers)
			ftype <- "2_row"
		} else if (nrow(d) == length(unique(markers[, 1]))) {
			marker <- markers[,1]
			ftype <- "1_row"
		} else {
			warning("don't know if this is a 'counts', '1_row' or a '2_row' file")
		}
	}
	markers <- data.frame(marker=marker, markers)
	rownames(d) <- marker
	out <- list(snp=d, markers=markers, geno=hdr, type=ftype, order=ordname)	
	class(out) <- "darter"
	out
}

show.darter <-  \(object) print(object)

print.darter <-  \(x, ...) {
	cat("class       :" , class(x), "\n")
	cat("order       :" , x$order, "\n")
	cat("type        :" , x$type, "\n")
	if (x$type == "1_row") {
		cat("markers     :" , nrow(x$snp), "\n")
	} else {
		cat("markers     :" , nrow(x$snp) / 2, "\n")	
	}
	cat("genotypes   :" , nrow(x$geno), "\n")
}

head.darter <-  \(x, n=6, ...) {
	cat("snp\n")
	print(x$snp[1:n, 1:12])
}




make_dart_1row <- function(x) {
	if (x$type=="1_row") return(x)
	stopifnot(x$type=="2_row")
#	||(x$type=="counts"))
	d <- x$snp
	i <- seq(1, nrow(d), 2)
	#dd <- d

	if (x$type == "counts") {
		# not OK as we also need marker specific thresholds for pres/abs
		threshold = 0.2
		f1 <- d[i,] / d[i+1,]
		f1[] <- pmin(1, f1)
		f1[d[i,] == 0] <- 0
		f1[d[i+1,] == 0 & d[i,] != 0] <- 1
		f1 <- ifelse(f1 < threshold, 0, 1)

		f2 <- d[i+1,] / d[i,	]
		f2[] <- pmin(1, f2)
		f2[d[i+1,] == 0] <- 0
		f2[d[i,] == 0 & d[i+1,] != 0] <- 1
		f2 <- ifelse(f2 < threshold, 0, 1)
		snp = 10 * f1 + f2
	} else {
		snp = 10 * d[i,] + d[i+1,]
	}
	
	snp[snp== 0] <- NA
	snp[snp==10] <- 0
	snp[snp==11] <- 2
	
	x$snp <- snp

	m <- x$markers[i+1,]
	m <- cbind (m[1], x$markers[i,2], m[-1])
	colnames(m)[2:3] <- paste0(colnames(x$markers)[2], c("Ref", "Alt"))
	x$markers <- m
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
	a <- lapply(1:nrow(d), \(i) {
			v <- d[i, ] + 2
			v[is.na(v)] <- 1
			rbind(first[v], second[v])
		})
	a <- do.call(rbind, a)
	colnames(a) <- colnames(d)
	rownames(a) <- paste0(rep(rownames(d), each=2), "_", 1:2)

	e$id <- 1:nrow(e)
	e <- rbind(e, e)
	e <- e[order(e$id), ]
	i <- seq(2, nrow(e), 2)
	e[i,2] <- e[i,3]
	e <- e[, -3]
	colnames(e)[2] <- "AlleleSequence"
	rownames(e) <- NULL

	x$snp <- a
	x$marker <- e
	x$type <- "2_row"
	x
}




old_dart_make_unique <- function(x) {
	n <- length(x$geno$genotype)
	if (n != length(unique(x$geno$genotype))) {
		x$geno$genotype <- make_unique_ids(x$geno$genotype)
		colnames(x$snp) <- make_unique_ids(colnames(x$snp))
	}
	if (x$type == "counts") {
		if (!is.null(x$geno$TargetID)) {
			x$geno$TargetID <- make_unique_ids(x$geno$TargetID)
		}	
	}
	x
}


write_dart <- function(x, filename) {
	s <- x$snp	
	g <- x$geno

	if (!is.null(g$targetID)) {
		colnames(s) <- x$geno$targetID
	} else {
		colnames(s) <- x$geno$ID
	}
	g$ID <- g$targetID <- NULL

	s <- cbind(x$marker, s)	
	s <- rbind(colnames(s), s)

	g <- t(g)
	d <- abs(dim(g) - c(0, ncol(s)))
#	if (x$type == "counts") {
#		#s[1, (d[2]+1):ncol(s)] <- g[nrow(g), ]
#		g <- cbind(matrix("*", d[1]-1, d[2]), g[-nrow(g), ])
#	} else {
		g <- cbind(matrix("*", d[1], d[2]), g)
#	}
	colnames(s) <- paste0("X", 1:ncol(s))
	colnames(g) <- paste0("X", 1:ncol(g))
	gs <- rbind(g, s)
	utils::write.table(gs, filename, na="-", col.names=FALSE, row.names=FALSE, sep=",")
}



dart_combine <- function(x, make_unique=FALSE) {

# needs to check if TID was used in read.dart
	stopifnot(is.list(x))
	out <- x[[1]]
	if (length(x) == 1) return(out)
	cns <- colnames(out$snp)
	colnames(out$snp) <- paste0("A", 1:nrow(out$geno))
	for (i in 2:length(x)) {
		y <- x[[i]]
		s <- nrow(out$geno)
		cns <- c(cns, colnames(y$snp))
		colnames(y$snp) <- paste0(LETTERS[i], (s+1):(s + nrow(y$geno))) 
		out$snp <- merge(out$snp, y$snp, by="row.names", all=TRUE)
		rownames(out$snp) <- out$snp[, 1]
		out$snp <- as.matrix(out$snp[, -1])
		out$geno <- rbind(out$geno, y$geno)
		out$order <- paste0(out$order, "_", y$order)
	}
	colnames(out$snp) <- cns
	out
}

