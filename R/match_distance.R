

clean_data_dist <- function(snp, genotypes, filename) {
#	filename <- fix_filename(filename)
	smps <- gsub("_D.$", "", genotypes$sample)
	snp <- remove_unknown_samples(snp, smps, verbose=FALSE)
#	snp <- fix_duplicate_names(snp, verbose=verbose)
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


hamming_distance <- function(x, genotypes, ref=TRUE) {
	stopifnot(x$type == "2_row")
	input <- clean_data_dist(x$snp, genotypes, "")
	d <- input$ref
	if (ref) {
		out <- sapply(1:ncol(d), \(i) colSums((d - d[,i]) != 0, na.rm=TRUE)) / nrow(d)
		colnames(out) <- rownames(out)
	} else {
		e <- input$field
		out <- sapply(1:ncol(e), \(i) colSums((d - e[,i]) != 0, na.rm=TRUE)) / nrow(d)
		rownames(out) <- input$ref.id
		colnames(out) <- input$field.id
	}
	out
}	


count_distance <- function(x, genotypes) {
	stopifnot(x$type == "counts")
	input <- clean_data_dist(x$snp, genotypes, "")
	d <- input$ref
	if (ref) {
		out <- sapply(1:ncol(d), \(i) colSums(abs(d - d[,i]), na.rm=TRUE))
		colnames(out) <- rownames(out)
	} else {
		e <- input$field
		out <- sapply(1:ncol(e), \(i) colSums(abs(d - e[,i]), na.rm=TRUE)) / nrow(d)
		rownames(out) <- input$ref.id
		colnames(out) <- input$field.id
	}
	out
}


