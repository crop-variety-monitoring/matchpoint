

count_fractions <- function(counts, mincounts=NULL) {
	n <- nrow(counts)
	i <- seq(1, n, 2)
	a <- counts[i,-1]
	ab <- a + counts[i+1,-1]
	f <- a / ab
	if (!is.null(mincounts)) {
		f[ab < mincounts] <- NA
	}
	data.frame(counts[i, 1, drop=FALSE], f, check.names=FALSE)
}


match_CDS <- function(x, genotypes, markers, method = "cor", snp_missing_rate=0.2, CDS_cutoff=0.5,
		mincounts=NULL, assign_threshold=NULL, filename, verbose=FALSE) {

	input <- matchpoint:::prepare_data(x, genotypes, markers, filename=filename, verbose=verbose)

	fract <- matchpoint:::count_fractions(input$snp, mincounts=mincounts)
	mr <- colSums(is.na(fract[,-1])) / nrow(fract)
	fract <- fract[, c(TRUE, mr <= snp_missing_rate)]

	if (method == "cor") {
		dmat <- 1 - cor(fract[,-1], use="pairwise.complete.obs")
	} else {
		stop("unknown method")
	}

	out_all <- 1 - dmat

#placeholder
	output <- list(metadata=data.frame(metric=c("distance method", "snp_missing_rate"), value=c(method, snp_missing_rate)))

	i <- which(colnames(out_all) %in% input$ref.id)
	out_match <- out_all[-i, i]
	ref_match <- out_all[i, i]
	
	if (is.null(assign_threshold)) {
		refnms <- genotypes$variety[match(colnames(ref_match), genotypes$sample)]
		dimnames(ref_match) <- list(refnms, refnms)
		pun <- matchpoint:::punity(1-ref_match, seq(0, 0.5, .01))
		# we want the last which.max
		assign_threshold <- 1 - pun[nrow(pun) - which.max(rev(pun[,"mean"])) + 1, "threshold"] 
		output[["metadata"]] <- rbind(output[["metadata"]], data.frame(metric="assign_threshold (computed)", value=assign_threshold)) 
	} else {
		output[["metadata"]] <- rbind(output[["metadata"]], data.frame(metric="assign_threshold (user input)", value=assign_threshold)) 	
	}
	
	d <- data.frame(field_id=rownames(out_match), 
					ref_id=rep(colnames(out_match), each=nrow(out_match)),
					variety="", CDS=as.vector(out_match))

	out_match <- data.frame(FieldSample=colnames(out_all)[-i], out_match, check.names=FALSE)
	rownames(out_match) <- NULL

#	mref_id <- gsub("_D.$", "", d$ref_id)
	mref_id <- d$ref_id
	i <- match(mref_id, genotypes$sample)
	if (any(is.na(i))) stop("check ID comparison")
	d$variety <- genotypes$variety[i]
	d$ref_Tid <- genotypes$TargetID[i]
#	mfld_id <- gsub("_D.$", "", d$field_id)
	mfld_id <- d$field_id
	i <- match(mfld_id, genotypes$sample)
	if (any(is.na(i))) stop("check ID comparison")
	d$field_Tid <- genotypes$TargetID[i]
	
	d <- d[order(d$field_id, -d$CDS), ]
# matches
	best <- d[!duplicated(d$field_id),]
	best$CDS <- round(best$CDS, 6)
	output[["best_match"]] <- best
	
#	vg <- matchpoint:::var_groups(out_all, assign_threshold, input$ref.id, genotypes$variety[match(input$ref.id, genotypes$sample)])
#	vg$variety <- sapply(strsplit(vg$ref_sample, ";"), 
#				\(i) paste(unique(genotypes$variety[match(i, genotypes$sample)]), collapse="; ")
#			)
	vg <- matchpoint:::old_var_groups(out_all, assign_threshold, input$ref.id)

	output[["all_match"]] <- vg

	ib <- d[d$CDS > CDS_cutoff[1], ] 
	ib$id_rank <- with(ib, stats::ave(CDS, field_id, FUN=\(x) rank(1-x, ties.method="min")))
	
	a <- stats::aggregate(ib[, "CDS", drop=FALSE], ib[, c("field_id", "variety")], max, na.rm=TRUE)
	a$var_rank <- with(a, stats::ave(CDS, field_id, FUN=\(x) rank(1-x, ties.method="min")))
	a$CDS <- NULL
	ib <- merge(ib, a, by=c("field_id",  "variety"))
	ib <- ib[order(ib$field_id, ib$id_rank), ]
	ib$CDS <- round(ib$CDS, 6)
	
	output[[paste0("CDS_id")]] <- ib
	
	dups <- duplicated(ib[, c("field_id", "var_rank")])
	output[[paste0("CDS_variety")]] <- ib[!dups, ]

	nr <- as.data.frame(table(field_id=d$field_id))
	rept <- stats::aggregate(d[, "CDS", drop=FALSE], d[, "field_id", drop=FALSE], 
			\(i) c(mean(i, na.rm=TRUE), stats::sd(i, na.rm=TRUE), length(stats::na.omit(i))))
	rept <- data.frame(rept[1], rept[[2]])
	colnames(rept)[-1] <- c("avg", "sd", "nr")

	output[["similarity"]] <- out_match
	out_all <- data.frame(ID=colnames(out_all), 1-out_all, check.names=FALSE)
	rownames(out_all) <- NULL
	output[["distance"]] <- out_all

	if (input$filename != "") {
		xlsx <- paste0(input$filename, "_CDS.xlsx")
		writexl::write_xlsx(output, xlsx, format_headers=FALSE)
		invisible(output)
	} else {
		output
	}

}


