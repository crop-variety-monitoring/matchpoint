

remove_sparse_records <- function(x, sample_mr=0.2, snp_mr=0.2, rows, verbose=TRUE) {

	if (!is.null(snp_mr)) {
		# keep these 
		i <- (rowSums(is.na(x)) / ncol(x)) <= snp_mr
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
		i <- apply(x, 1, \(i) length(stats::na.omit(unique(i))))
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
				nr <- if(x$type == "1_row") length(i) else length(i)/2
				message(paste("   removed", nr, "markers with no variation"))
			}
		}
	}
	
	if (!is.null(sample_mr)) {
		# keep these 
		i <- (colSums(is.na(x)) / nrow(x)) <= sample_mr
		if (any(!i)) {
			x <- x[, c(TRUE, i)]
			if (verbose && any(!i)) {
				message(paste("   removed", sum(!i), "sample(s) with too many missing values"))
			}
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


..fix_duplicate_names <- function(x, suf="_D", verbose=FALSE) {
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
	i <- match(colnames(x$snp), sample.id)
	j <- is.na(i)
	if (any(j)) {
		if (verbose) message(paste("  ", sum(j), "unmatched genotype(s) removed"))
		x$snp <- x$snp[, !j]
		x$geno <- x$geno[!j, ]
	}
	x
}


prepare_data <- function(x, genotypes, match_field = "", markers=NULL, filename, sample_mr=NULL, snp_mr=NULL, verbose=FALSE) {

	if (match_field=="" || is.na(match_field)) {
		stop("supply a value for 'match_field'")
	}

	filename <- matchpoint:::fix_filename(filename)
	if (x$type != "1_row") {
		nr <- 1
	} else {
		nr <- 2
	}

	if ("order" %in% names(genotypes)) {
		genotypes <- merge(x$geno, genotypes, by.x=c("order", "ID"), by.y=c("order", match_field))
		if (nrow(unique(genotypes[, c("order", "ID")])) != nrow(genotypes)) {
			message("genotypes file order/sample numbers are not unique")
		}
	} else {
		genotypes <- merge(x$geno, genotypes, by.x="ID", by.y=match_field)
		if (length(unique(genotypes$ID)) != nrow(genotypes)) {
			message("genotypes file sample numbers are not unique")
		}		
	}

	x <- matchpoint:::remove_unknown_samples(x, genotypes$ID, verbose=verbose)
	#snp <- matchpoint:::fix_duplicate_names(snp, suf=dupsuf, verbose=verbose)


	if (!(is.null(sample_mr) && is.null(snp_mr))) {
		x$snp <- remove_sparse_records(x$snp, sample_mr, snp_mr, rows=nr, verbose)
	}
	
	i <- match(x$geno$ID, genotypes$ID)
	nomatch <- which(is.na(i))
	if (length(nomatch) > 0) {
		nms <- x$geno$ID[nomatch]
		match(nms, names(x$snp))
		x$geno$genotype <- x$geno$genotype[-nomatch, ]
	}

	cns <- colnames(x$snp)
	ref.id <- cns[genotypes$reference[i]]
	field.id <- cns[!genotypes$reference[i]]
	gds <- NULL

	if (!is.null(markers)) {
		markers <- markers[, c("MarkerName", "Chr", "Pos")]
		imark <- match(toupper(rownames(x$snp)), toupper(markers$MarkerName))
		if (any(is.na(imark))) {
			unk <- x$snp[is.na(imark), 1]
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
	list(snp=x$snp, markers=markers, ref.id=ref.id, field.id=field.id, filename=filename, gds=gds)
}

