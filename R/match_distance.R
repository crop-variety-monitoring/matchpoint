
ibs_distance <- function(ref, fld=NULL) {
	if (is.null(fld)) {
		ref <- as.matrix(ref)
		ref[] <- match(ref, c(0,2,1)) - 1

		out <- matrix(NA, ncol=ncol(ref), nrow=ncol(ref))
		colnames(out) <- rownames(out) <- colnames(ref)
		nc <- ncol(ref)
		for (i in 1:nc) {
			j <- i:nc
			out[i, j] <- out[j, i] <- colSums(2 - abs(ref[,j,drop=FALSE] - ref[,i]), na.rm=TRUE)
		}
	} else {
		ref <- as.matrix(ref)
		fld <- as.matrix(fld)
		ref[] <- match(ref, c(0,2,1)) - 1
		fld[] <- match(fld, c(0,2,1)) - 1
		out <- t(sapply(1:ncol(fld), \(i) colSums(2 - abs(ref[,j,drop=FALSE] - ref[,i]), na.rm=TRUE)))
		rownames(out) <- colnames(fld)
	}
	1 - (out / (2 * nrow(ref)))
}	


hamming_distance <- function(ref, fld=NULL) {
	if (is.null(fld)) {
		ref <- as.matrix(ref)
		out <- matrix(NA, ncol=ncol(ref), nrow=ncol(ref))
		colnames(out) <- rownames(out) <- colnames(ref)
		nc <- ncol(ref)
		for (i in 1:nc) {
			j <- i:nc
			out[i, j] <- out[j, i] <- colSums(ref[,j,drop=FALSE] != ref[,i], na.rm=TRUE)
		}
	} else {
		ref <- as.matrix(ref)
		fld <- as.matrix(fld)
		out <- t(sapply(1:ncol(fld), \(i) colSums(ref != fld[,i], na.rm=TRUE)))
		rownames(out) <- colnames(fld)
	}
	out 
	# / nrow(ref)
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


	d$id_rank <- stats::ave(d[[name]], d$field_id, FUN=\(x) rank(1-x, ties.method="min"))
	a <- stats::aggregate(d[, "dist", drop=FALSE], d[, c("field_id", "variety")], max, na.rm=TRUE)
	a$var_rank <- stats::ave(a[[name]], a$field_id, FUN=\(x) rank(1-x, ties.method="min"))
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




match_distance <- function(x, genotypes, missing_rate=0.2, comp_all=FALSE, filename="", verbose=TRUE) {

	input <- matchpoint:::prepare_data(x, genotypes, missing_rate=missing_rate, filename=filename, verbose=verbose)

	if (x$type == "counts") {
		dst <- poppr::nei.dist(as.matrix(input$snp[, -1]))
		if (!comp_all) {
			fld <- colnames(input$snp) %in% input$field.id
			ref <- colnames(input$snp) %in% input$ref.id
			dst <- as.matrix(dst)[ref, fld]
		}
		if (filename != "") filename <- paste0(filename, "_NEI.xlsx")
	} else if (x$type == "2_row") {
		if (comp_all) {
			dst <- matchpoint:::hamming_distance(input$snp[, -1])
		} else {
			fld <- input$snp[, colnames(input$snp) %in% input$field.id]
			ref <- input$snp[, colnames(input$snp) %in% input$ref.id]
			dst <- matchpoint:::hamming_distance(ref, fld)
		}
		if (filename != "") filename <- paste0(filename, "_HAM.xlsx")
	} else if (x$type == "1_row") {
		if (comp_all) {
			dst <- matchpoint:::ibs_distance(input$snp[, -1])
		} else {
			fld <- input$snp[, colnames(input$snp) %in% input$field.id]
			ref <- input$snp[, colnames(input$snp) %in% input$ref.id]
			dst <- matchpoint:::ibs_distance(ref, fld)
		}
		if (filename != "") filename <- paste0(filename, "_IBS.xlsx")
	}

	meta <- data.frame(metric="distance", 
				value=ifelse(x$type == "counts", "euclidean", "hamming"))
	output <- get_output(dst, genotypes, input, "dist", meta, comp_all)

	if (filename != "") {
		output$distance <- as.data.frame(output$distance, check.names=FALSE)
		output$similarity <- as.data.frame(output$similarity, check.names=FALSE)
		
		writexl::write_xlsx(output, filename, format_headers=FALSE)
		invisible(output)
	} else {
		output
	}
}


