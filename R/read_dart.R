
read_dart <- function(filename) {

	r <- data.frame(data.table::fread(filename, header=FALSE), check.names=FALSE)
	ordname <- gsub("^Report_", "", basename(filename))
	ordname <- gsub("\\_SNP.csv$", "", ordname)

	srow <- sum(r[1:30, 1] == "*") + 1
	scol <- sum(r[1, 1:50] == "*") + 1

	marker <- r[(srow+1):nrow(r), 1:(scol-1), ]
	colnames(marker) <- r[srow, 1:(scol-1)]

	hdr <- data.frame(t(r[1:srow, scol:ncol(r)]))
	cns <- c("order", "plate_barcode", "sample_reproducibility", "plate_row", "plate_col", "plate_id", "genotype")
	if (ncol(hdr) == 7) {
		colnames(hdr) <- cns
	} else 	if (ncol(hdr) == 8) {
		colnames(hdr) <- c(cns, "TargetID")
	}
	
	d <- as.matrix(r[(srow+1):nrow(r), scol:ncol(r)])
	v <- as.vector(d)
	v[v == "-"] <- NA
	d <- matrix(as.integer(v), nrow=nrow(d), ncol=ncol(d))	
	d <- data.frame(r[(srow+1):nrow(r), 1:2], d, check.names=FALSE)

	if (grepl("_Counts", filename) || isTRUE(any(d[,-c(1:2)] > 2))) { 
		ids <- trimws(hdr[, ncol(hdr)-1])
		colnames(d) <- c(colnames(marker)[1:2], ids)
		out <- list(snp=d, marker=marker, geno=hdr, type="counts")
	} else {
		i <- seq(1, nrow(d), 2)
		if (grepl("_2row", filename) || (all(d[i,1] == d[i+1,1]))) {
			if (ncol(hdr) == 8) {
				ids <- trimws(hdr[,ncol(hdr)-1])		
			} else {
				ids <- trimws(hdr[,ncol(hdr)])
			}
			colnames(d) <- c(colnames(marker)[1:2], ids)
			out <- list(snp=d, marker=marker, geno=hdr, type="2_row", order=ordname)
		} else if (nrow(d) == length(unique(d[, 1]))) {
			d <- d[,-2]
			ids <- trimws(hdr[,ncol(hdr)])
			colnames(d) <- c(colnames(marker)[1], ids)
			out <- list(snp=d, marker=marker, geno=hdr, type="1_row", order=ordname)	
		} else {
			warning("don't know if this is a Count, 1_row or a 2_row file")
			out <- list(snp=d, marker=marker, geno=hdr, type="???", order=ordname)		
		}
	}
	class(out) <- "darter"
	out
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
	x$snp <- a
	x$type=="2_row"
	x
}




dart_make_unique <- function(x) {
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
	if (x$type == "1_row") {
		idcols <- 1
	} else {
		idcols <- 1:2	
	}
	s <- cbind(x$marker, x$snp[, -idcols])
	s <- rbind(colnames(s), s)
	g <- t(x$geno)
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
