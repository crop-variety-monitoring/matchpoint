..clean_data_dist <- function(snp, genotypes, filename, verbose=TRUE) {
#	filename <- fix_filename(filename)
	smps <- gsub("_D.$", "", genotypes$sample)
	snp <- remove_unknown_samples(snp, smps, verbose=verbose)
#?	snp <- fix_duplicate_names(snp, verbose=verbose)
	cns <- colnames(snp)[-1]
	snp <- as.matrix(snp[, -1])
	colnames(snp) <- cns
	i <- match(cns, gsub("_D.$", "", genotypes$sample))
	ref.id   <- cns[genotypes$reference[i]]
	field.id <- cns[!genotypes$reference[i]]
	i <- cns %in% ref.id
	ref <- snp[,i]
	field <- snp[,!i]
	list(ref=ref, field=field, ref.id=ref.id, field.id=field.id, filename=filename)
}


hamming_distance <- function(ref, fld=NULL) {
	if (is.null(fld)) {
		out <- sapply(1:ncol(ref), \(i) colSums((ref - ref[,i]) != 0, na.rm=TRUE)) / nrow(ref)
		colnames(out) <- rownames(out)
	} else {
		out <- t(sapply(1:ncol(fld), \(i) colSums((ref - fld[,i]) != 0, na.rm=TRUE)) / nrow(ref))
		rownames(out) <- colnames(fld)
	}
	out
}	


counts_distance <- function(ref, fld=NULL) {
	if (ref) {
		out <- sapply(1:ncol(ref), \(i) colSums(abs(ref - ref[,i]), na.rm=TRUE))
		colnames(out) <- rownames(out)
	} else {
		out <- t(sapply(1:ncol(fld), \(i) colSums(abs(ref - fld[,i]), na.rm=TRUE)) / nrow(ref))
		rownames(out) <- colnames(fld)
	}
	out
}



match_distance <- function(x, genotypes, compare="ref2fld", missing_rate=0.25, filename="", verbose=TRUE) {

	compare <- match.arg(tolower(compare), c("ref2fld", "ref2ref", "fld2fld", "all"))
	
	input <- prepare_data(x, genotypes, missing_rate=missing_rate, filename=filename, verbose=verbose)

	if (compare == "fld2fld") {
		ref <- input$snp[, colnames(input$snp) %in% input$field.id]	
	} else if (compare == "all") {
		ref <- input$snp
	} else {
		ref <- input$snp[, colnames(input$snp) %in% input$ref.id]
	}
	
	if (compare == "ref2fld") {
		fld <- input$snp[, colnames(input$snp) %in% input$field.id]
	} else {
		fld <- NULL
	}
	if (x$type == "counts") {
		d <- counts_distance(ref, fld)
	} else if (x$type == "2_row") {
		d <- hamming_distance(ref, fld)
	} else { # make 2_row?
		stop("need a counts or 2_row object")
	}
	out = list(dist=d)
	if (is.null(fld)) {
		diag(d) <- NA
	}
	i <- apply(d, 1, which.min)

	if (substr(compare, 1, 3) == "ref") {
		v <- genotypes$variety[match(colnames(d), genotypes$sample)]
		out$best_match <- data.frame(sample=rownames(d), ref=colnames(d)[i], variety=v[i])
	} else {
		out$best_match <- data.frame(sample=rownames(d), ref=colnames(d)[i])	
	}
	out$meta <- data.frame(metric="distance", value=ifelse(x$type == "counts", "euclidean", "hamming"))

	if (filename != "") {
		writexl::write_xlsx(out, paste0(filename, "xlsx"), format_headers=FALSE)
		invisible(out)
	} else {
		out
	}
}


