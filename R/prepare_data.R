

fix_filename <- function(filename) {
	filename <- trimws(filename)
	if (filename != "") {
		filename <- tools::file_path_sans_ext(filename)
		outdir <- dirname(filename)
		dir.create(outdir, FALSE, TRUE)
	}
	filename
}


fix_duplicate_names <- function(x, verbose=FALSE) {
	cn <- colnames(x)
	tab <- table(cn)
	dups <- tab[tab > 1]
	if (length(dups) > 0) {
		if (verbose) message("   adding _D1 _D2 to duplicates")
		for (i in 1:length(dups)) {
			n <- dups[i]
			sfx <- paste0("_D", 1:n)
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


prepare_data <- function(x, genotypes, markers=NULL, filename, missing_rate=NULL, verbose=FALSE) {
	filename <- fix_filename(filename)
	snp <- remove_unknown_samples(x$snp, genotypes$sample, verbose=verbose)
	snp <- fix_duplicate_names(snp, verbose=verbose)

	if (!(is.null(missing_rate) || is.na(missing_rate))) {
		# keep these 
		i <- round(colSums(is.na(x$snp[,-1])) / nrow(x$snp), 2) <= missing_rate
		x$snp <- x$snp[, c(TRUE, i)]
		if (verbose && any(!i)) {
			message(paste("   removed", sum(!i), "sample(s) that had too many missing values"))
		}
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

