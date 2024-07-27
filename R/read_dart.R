
read_dart <- function(filename) {

	r <- data.frame(data.table::fread(filename, header=FALSE), check.names=FALSE)
	ordname <- gsub("^Report_", "", basename(filename))
	ordname <- gsub("\\_SNP_2row.csv$|\\_SNP.csv$|\\_Counts.csv$", "", ordname)

	srow <- sum(r[1:30, 1] == "*") + 1
	scol <- sum(r[1, 1:50] == "*") + 1

	markers <- r[(srow+1):nrow(r), 1:(scol-1), ]
	colnames(markers) <- r[srow, 1:(scol-1)]

	hdr <- data.frame(t(r[1:srow, scol:ncol(r)]))
	cns <- c("order", "plate_barcode", "sample_reproducibility", "plate_row", "plate_col", "plate_id", "genotype")
	if (ncol(hdr) == 7) {
		colnames(hdr) <- cns
	} else 	if (ncol(hdr) == 8) {
		colnames(hdr) <- c(cns, "TargetID")
	}
	ids <- paste0("G", 1:nrow(hdr))
	hdr <- data.frame(ID=ids, hdr)
	d <- as.matrix(r[(srow+1):nrow(r), scol:ncol(r)])
	v <- as.vector(d)
	v[v == "-"] <- NA
	d <- matrix(as.integer(v), nrow=nrow(d), ncol=ncol(d))	
	colnames(d) <- ids

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
	d <- data.frame(marker=marker, d)
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
	a <- lapply(1:nrow(d), \(i) {
			v <- dd[i, ] + 2
			v[is.na(v)] <- 1
			rbind(first[v], second[v])
		})
	a <- do.call(rbind, a)
	a <- data.frame(rep(d[,1], each=2), as.vector(t(e[,2:3])), a)
	names(a) <- c("MarkerName", "AlleleSequence", names(d)[-1])

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
	x$type=="2_row"
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
	s <- x$snp[, -1]
	colnames(s) <- x$geno$genotype
	s <- cbind(x$marker, s)	
	s <- rbind(colnames(s), s)
	
	g <- x$geno
	
	if (is.null(g$TargetID)) g$genotype <- NULL
	g$ID <- g$TargetID <- NULL
	g <- t(g)

	d <- abs(dim(g) - c(0, ncol(s)))
	if (x$type == "counts") {
		s[1, (d[2]+1):ncol(s)] <- g[nrow(g), ]
		g <- cbind(matrix("*", d[1]-1, d[2]), g[-nrow(g), ])
	} else {
		g <- cbind(matrix("*", d[1], d[2]), g)
	}
	colnames(s) <- paste0("X", 1:ncol(s))
	colnames(g) <- paste0("X", 1:ncol(g))
	gs <- rbind(g, s)
	utils::write.table(gs, filename, na="-", col.names=FALSE, row.names=FALSE, sep=",")
}



dart_combine <- function(x) {
	stopifnot(is.list(x))
	out <- x[[1]]
	if (length(x) == 1) return(out)
	
	out$geno$ID <- paste0("G", 1:nrow(out$geno))
	names(out$snp)[-1] <- out$geno$ID 
	
	for (i in 2:length(x)) {
		y <- x[[i]]
		s <- nrow(out$geno)
		y$geno$ID <- paste0("G", (s+1):(s + nrow(y$geno)))
		names(y$snp)[-1] <- y$geno$ID 
		
		out$snp <- merge(out$snp, y$snp, by="marker", all=TRUE)
		out$geno <- rbind(out$geno, y$geno)
		out$order <- paste0(out$order, "_", y$order)
	}
	out
}

