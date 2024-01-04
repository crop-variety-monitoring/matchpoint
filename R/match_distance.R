
hamming_distance <- function(ref, fld=NULL) {
	if (is.null(fld)) {
		out <- matrix(NA, ncol=ncol(ref), nrow=ncol(ref))
		colnames(out) <- rownames(out) <- colnames(ref)
		nc <- ncol(ref)
		for (i in 1:nc) {
			j <- i:nc
			out[i, j] <- out[j, i] <- colSums(ref[,j,drop=FALSE] != ref[,i], na.rm=TRUE)
		}
	} else {
		out <- t(sapply(1:ncol(fld), \(i) colSums(ref != fld[,i], na.rm=TRUE)) / nrow(ref))
		rownames(out) <- colnames(fld)
	}
	out
}	


counts_distance <- function(ref, fld=NULL) {
	if (is.null(fld)) {
		out <- sapply(1:ncol(ref), \(i) colSums(abs(ref - ref[,i]), na.rm=TRUE))
		colnames(out) <- rownames(out)
	} else {
		out <- t(sapply(1:ncol(fld), \(i) colSums(abs(ref - fld[,i]), na.rm=TRUE)) / nrow(ref))
		rownames(out) <- colnames(fld)
	}
	out
}



get_output <- function(x, genotypes, input, name, meta=NULL, comp_all=TRUE) {

	if (comp_all) {
		m <- x[rownames(x) %in% input$field.id, colnames(x) %in% input$ref.id]
	} else {
		m <- x
	}
	
	d <- data.frame(field_id=rownames(m), 
			ref_id=rep(colnames(m), each=nrow(m)),
			variety="", dist=as.vector(m))
	i <- match(d$ref_id, genotypes$sample)
	d$variety <- genotypes$variety[i]
	d$ref_Tid <- genotypes$TargetID[i]
	i <- match(d$field_id, genotypes$sample)
	d$field_Tid <- genotypes$TargetID[i]
	d[with(d, order(field_id, -dist)),]
	best <- d[!duplicated(d$field_id),]
	best$dist <- round(best$dist, 6)


	d$id_rank <- ave(d[[name]], d$field_id, FUN=\(x) rank(1-x, ties.method="min"))
	a <- stats::aggregate(d[, "dist", drop=FALSE], d[, c("field_id", "variety")], max, na.rm=TRUE)
	a$var_rank <- ave(a[[name]], a$field_id, FUN=\(x) rank(1-x, ties.method="min"))
	a[[name]] <- NULL
	d <- merge(d, a, by=c("field_id",  "variety"))
	d <- d[order(d$field_id, d$id_rank), ]
	d[[name]] <- round(d[[name]], 6)
	dups <- duplicated(d[, c("field_id", "var_rank")])

	if (comp_all) {
		out <- list(meta=meta, best_match=best, d, d[!dups, ], distance=x, similarity=1-m)
	} else {
		out <- list(meta=meta, best_match=best, d, d[!dups, ], similarity=1-m)	
	}
	names(out)[3:4] <- paste0(name, c("_id", "_variety"))
	out
}




match_distance <- function(x, genotypes, missing_rate=0.25, comp_all=FALSE, filename="", verbose=TRUE) {

	input <- matchpoint:::prepare_data(x, genotypes, missing_rate=missing_rate, filename=filename, verbose=verbose)


	if (x$type == "counts") {
		if (comp_all) {
			dst <- counts_distance(input$snp[, -1])
		} else {
			fld <- input$snp[, colnames(input$snp) %in% input$field.id]
			ref <- input$snp[, colnames(input$snp) %in% input$ref.id]
			dst <- counts_distance(ref, fld)
		}
	} else if (x$type == "2_row") {
		if (comp_all) {
			dst <- matchpoint:::hamming_distance(input$snp[, -1])
		} else {
			fld <- input$snp[, colnames(input$snp) %in% input$field.id]
			ref <- input$snp[, colnames(input$snp) %in% input$ref.id]
			dst <- matchpoint:::hamming_distance(ref, fld)
		}
	} else { # make 2_row?
		stop("need a counts or 2_row object")
	}

	meta <- data.frame(metric="distance", 
				value=ifelse(x$type == "counts", "euclidean", "hamming"))
	output <- get_output(dst, genotypes, input, "dist", meta, comp_all)

	if (filename != "") {
		output$distance <- as.data.frame(output$distance, check.names=FALSE)
		output$similarity <- as.data.frame(output$similarity, check.names=FALSE)
		writexl::write_xlsx(output, paste0(filename, "_HAM.xlsx"), format_headers=FALSE)
		invisible(output)
	} else {
		output
	}
}


