

remove_sparse_records <- function(x, sample_mr=0.2, snp_mr=0.2, rows, verbose=TRUE) {
	# keep these 
	i <- (rowSums(is.na(x[,-1])) / ncol(x)) <= snp_mr
	if (any(!i)) {
		if (rows != 1) {
			j <- which(!i)
			k <- rep_len(1:2, length.out=length(j)) 
			if (!all(j[k==1] == (j[k==2]-1))) {
				stop("check this")
			}
		}
		x <- x[i, ]
		if (verbose) {
			message(paste("   removed", sum(!i), "marker(s) with too many missing values"))
		}
	}

	i <- apply(x[,-1], 1, \(i) length(na.omit(unique(i))))
	i <- which(i < 2) # could be 0 or 1
	if (length(i) > 0) {
		if (rows != 1) {
			j <- which(!i)
			k <- rep_len(1:2, length.out=length(j)) 
			if (!all(j[k==1] == (j[k==2]-1))) {
				stop("check this")
			}
		}
		x <- x[-i, ]
		if (verbose) {
			nr <- ifelse(x$type == "1_row", length(i), length(i)/2) 
			message(paste("   removed", nr, "markers with no variation"))
		}
	}

	# keep these 
	i <- (colSums(is.na(x[,-1])) / nrow(x)) <= sample_mr
	if (any(!i)) {
		x <- x[, c(TRUE, i)]
		if (verbose && any(!i)) {
			message(paste("   removed", sum(!i), "sample(s) with too many missing values"))
		}
	}
	x
}


fix_filename <- function(filename) {
	filename <- trimws(filename)
	if (filename != "") {
		filename <- tools::file_path_sans_ext(filename)
		outdir <- dirname(filename)
		dir.create(outdir, FALSE, TRUE)
	}
	filename
}


fix_duplicate_names <- function(x, suf="_D", verbose=FALSE) {
	cn <- colnames(x)
	tab <- table(cn)
	dups <- tab[tab > 1]
	if (length(dups) > 0) {
		if (verbose) message(paste0("   adding ", suf, "* to duplicates"))
		for (i in 1:length(dups)) {
			n <- dups[i]
			sfx <- paste0(suf, 1:n)
			nm <- names(n)
			cn[cn == nm] <- paste0(nm, sfx)
		}
		colnames(x) <- cn
	}
	x
}

remove_unknown_samples <- function(x, sample.id, verbose=FALSE) {
	cns <- colnames(x)
	i <- match(cns[-1], sample.id)
	j <- is.na(i)
	if (any(j)) {
		if (verbose) message(paste("  ", sum(j), "unmatched genotype(s) removed"))
		k <- which(j) + 1
		x <- x[, -k]
		colnames(x) <- cns[-k]
	}
	x
}


prepare_data <- function(x, genotypes, markers=NULL, filename, missing_rate=NULL, dupsuf="_D", verbose=FALSE) {

	filename <- matchpoint:::fix_filename(filename)
	if (x$type != "1_row") {
		nr <- 1
		x$snp <- x$snp[,-2]
	} else {
		nr <- 2
	}
	
	snp <- matchpoint:::remove_unknown_samples(x$snp, genotypes$sample, verbose=verbose)
	snp <- matchpoint:::fix_duplicate_names(snp, suf=dupsuf, verbose=verbose)

	if (!(is.null(missing_rate) || is.na(missing_rate))) {
		snp <- remove_sparse_records(snp, missing_rate, missing_rate, rows=nr, verbose)
	}
	
	cns <- colnames(snp)[-1]
	i <- match(cns, genotypes$sample)
	ref.id <- cns[genotypes$reference[i]]
	field.id <- cns[!genotypes$reference[i]]
	gds <- NULL

	if (!is.null(markers)) {
		markers <- markers[, c("MarkerName", "Chr", "Pos")]
		imark <- match(toupper(snp[,1]), toupper(markers$MarkerName))
		if (any(is.na(imark))) {
			unk <- snp[is.na(imark), 1]
			message("   unknown markers in snp file:\n   ", paste(unk, collapse=", "))
		}
		markers <- markers[imark, ]
		markers$Chr[is.na(markers$Chr)] <- ""
		markers$Pos[is.na(markers$Pos)] <- ""
		if (filename == "") {
			gds <- paste0(tempfile(), "_geno.gds")	
		} else {
			gds <- paste0(filename, "_geno.gds")
		}
	}
	list(snp=snp, markers=markers, ref.id=ref.id, field.id=field.id, filename=filename, gds=gds)
}

