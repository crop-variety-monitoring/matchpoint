
read_dart <- function(filename) {

	r <- data.frame(data.table::fread(filename, header=FALSE), check.names=FALSE)
	ordname <- gsub("^Report_", "", basename(filename))
	ordname <- gsub("\\_SNP.csv$", "", ordname)

	srow <- sum(r[1:30, 1] == "*") + 1
	scol <- sum(r[1, 1:50] == "*") + 1

	marker <- r[(srow+1):nrow(r), 1:(scol-1), ]
	colnames(marker) <- r[srow, 1:(scol-1)]

	hdr <- data.frame(t(r[1:srow, scol:ncol(r)]))
	if (ncol(hdr) == 7) {
		colnames(hdr) <- c("order", "plate_barcode", "sample_reproducibility", "extract_row", "extract_col", "extract_plate", "TargetID")
	} else 	if (ncol(hdr) == 8) {
		colnames(hdr) <- c("order", "plate_barcode", "sample_reproducibility", "extract_row", "extract_col", "extract_plate", "TargetID", "genotype")
	}
	
	d <- as.matrix(r[(srow+1):nrow(r), scol:ncol(r)])
	v <- as.vector(d)
	v[v == "-"] <- NA
	d <- matrix(as.integer(v), nrow=nrow(d), ncol=ncol(d))	
	d <- data.frame(r[(srow+1):nrow(r), 1:2], d, check.names=FALSE)
	if (isTRUE(any(d[,-c(1:2)] > 2))) { 
		ids <- trimws(hdr[, ncol(hdr)-1])
		colnames(d) <- c(colnames(marker)[1:2], ids)
		return(list(snp=d, marker=marker, geno=hdr, type="counts"))
	} 
	i <- seq(1, nrow(d), 2)
	if (all(d[i,1] == d[i+1,1])) {
		ids <- trimws(hdr[,ncol(hdr)])
		colnames(d) <- c(colnames(marker)[1:2], ids)
		list(snp=d, marker=marker, geno=hdr, type="2_row", order=ordname)
	} else if (nrow(d) == length(unique(d[, 1]))) {
		d <- d[,-2]
		ids <- trimws(hdr[,ncol(hdr)])
		colnames(d) <- c(colnames(marker)[1], ids)
		list(snp=d, marker=marker, geno=hdr, type="1_row", order=ordname)	
	} else {
		warning("don't know if this is a Count, 1_row or a 2_row file")
		list(snp=d, marker=marker, geno=hdr, type="???", order=ordname)		
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


