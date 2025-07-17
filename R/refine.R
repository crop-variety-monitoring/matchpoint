

refine_reference <- function(dst, lump=.01, split=.05) {

	dcn <- colnames(dst)
	nothr1 <- matchpoint:::nngb_dist(dst)
	nsame1 <- matchpoint:::min_self_dist(dst)
	splmp <- matchpoint:::split_lump(dst, lump, split) 
	out <- list(varieties=data.frame(old_name=splmp$old, new_name=splmp$new), parameters=data.frame(splmp[c("split", "lump")]))
	colnames(dst) <- rownames(dst) <- splmp$new
	
	out$near_other_old <- nothr1
	out$near_other_new <- matchpoint:::nngb_dist(dst)

	nsame2 <- matchpoint:::min_self_dist(dst)
	near_same <- data.frame(nsame1, nsame2)
	colnames(near_same) <- c("name_old", "value_old", "name_new", "value_new")
	out$near_same <- near_same

	out
}




finish_refine <- function(x, dstm, genotypes, match_field, ref_lump, ref_split, filename) {

#	gtypes <- genotypes[genotypes$reference, ]

	ids <- colnames(dstm)
	nms <- genotypes$variety[match(colnames(dstm), genotypes[,match_field])]
	d2 <- dstm
	colnames(d2) <- rownames(d2) <- nms

#	gtype <- x$geno$gtype[match(colnames(dstm), x$geno$ID)]
#	refnms <- gtypes$variety[match(gtype, gtypes$sample)]
#	dimnames(dstm) <- list(refnms, refnms)
	pun <- matchpoint:::punity(d2, seq(0, 0.5, .01))
		# we want the last which.max
	if (is.null(ref_split) || (is.na(ref_split))) {
		ref_split <- pun[nrow(pun) - which.max(rev(pun[,"mean"])) + 1, "threshold"] 
	}
	if (is.null(ref_lump) || is.na(ref_lump)) ref_lump <- ref_split / 2
	
	output <- matchpoint:::refine_reference(d2, lump=ref_lump, split=ref_split)

	output$distance <- data.frame(newvar=output$varieties$new, dstm, check.names=FALSE)

	colnames(d2) <- rownames(d2) <- output$varieties$new_name
	pun2 <- matchpoint:::punity(d2, seq(0, 0.5, .01))
	output$punity_original <- data.frame(pun)
	output$punity_refined <- data.frame(pun2)
	ref_split2 <- pun2[nrow(pun2) - which.max(rev(pun2[,"mean"])) + 1, "threshold"] 
	ref_lump2 <- ref_split2 / 2
	output$parameters <- rbind(output$parameters, data.frame(split=ref_split2, lump=ref_lump2))
	
	newnms <- data.frame(ID=ids, variety=output$varieties$new_name)
	colnames(genotypes)[colnames(genotypes) == "variety"] <- "old_variety"
	output$genotypes <- merge(genotypes, newnms, by.x=match_field, by.y="ID", all.x=TRUE)
	#output$genotypes$variety[is.na(output$variety)] <- ""
	# to avoid "Coercing column plate_barcode from int64 to double" warning
	if (!is.null(output$genotypes$plate_barcode)) output$genotypes$plate_barcode <- as.character(output$genotypes$plate_barcode)

	filename <- trimws(filename)
	if (filename == "") {
		return(output)
	} else {
		xlsx <- paste0(filename, ".xlsx")
		dir.create(dirname(xlsx), FALSE, TRUE)
		writexl::write_xlsx(output, xlsx, format_headers=FALSE)
		utils::write.csv(output$genotypes, paste0(filename, "_variety-info.csv"), row.names=FALSE)
		saveRDS(output, paste0(filename, ".rds"))
		invisible(output)
	}
}


old_finish_refine <- function(x, dstm, genotypes, match_field, ref_lump, ref_split, filename) {

#	gtypes <- genotypes[genotypes$reference, ]

	ids <- colnames(dstm)
	nms <- genotypes$variety[match(colnames(dstm), genotypes[,match_field])]
	d2 <- dstm
	colnames(d2) <- rownames(d2) <- nms

#	gtype <- x$geno$gtype[match(colnames(dstm), x$geno$ID)]
#	refnms <- gtypes$variety[match(gtype, gtypes$sample)]
#	dimnames(dstm) <- list(refnms, refnms)
	pun <- matchpoint:::punity(d2, seq(0, 0.5, .01))
		# we want the last which.max
	if (is.null(ref_split) || (is.na(ref_split))) {
		ref_split <- pun[nrow(pun) - which.max(rev(pun[,"mean"])) + 1, "threshold"] 
	}
	if (is.null(ref_lump) || is.na(ref_lump)) ref_lump <- ref_split / 2
	
	output <- matchpoint:::refine_reference(d2, lump=ref_lump, split=ref_split)

	output$distance <- data.frame(newvar=output$varieties$new, dstm, check.names=FALSE)

	colnames(d2) <- rownames(d2) <- output$varieties$new_name
	pun2 <- matchpoint:::punity(d2, seq(0, 0.5, .01))
	output$punity_original <- data.frame(pun)
	output$punity_refined <- data.frame(pun2)
	ref_split2 <- pun2[nrow(pun2) - which.max(rev(pun2[,"mean"])) + 1, "threshold"] 
	ref_lump2 <- ref_split2 / 2
	output$parameters <- rbind(output$parameters, data.frame(split=ref_split2, lump=ref_lump2))
	
	newnms <- data.frame(ID=ids, variety=output$varieties$new_name)
	colnames(genotypes)[colnames(genotypes) == "variety"] <- "old_variety"
	output$genotypes <- merge(genotypes, newnms, by.x=match_field, by.y="ID", all.x=TRUE)
	#output$genotypes$variety[is.na(output$variety)] <- ""
	# to avoid "Coercing column plate_barcode from int64 to double" warning
	if (!is.null(output$genotypes$plate_barcode)) output$genotypes$plate_barcode <- as.character(output$genotypes$plate_barcode)

	filename <- trimws(filename)
	if (filename == "") {
		return(output)
	} else {
		xlsx <- paste0(filename, ".xlsx")
		dir.create(dirname(xlsx), FALSE, TRUE)
		writexl::write_xlsx(output, xlsx, format_headers=FALSE)
		utils::write.csv(output$genotypes, paste0(filename, "_variety-info.csv"), row.names=FALSE)
		saveRDS(output, paste0(filename, ".rds"))
		invisible(output)
	}
}

