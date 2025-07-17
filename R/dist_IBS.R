
run_IBS <- function(input, MAF_cutoff=0.05, SNP_Missing_Rate=0.2, Sample_Missing_Rate=0.2, Sample_Heterozygosity_Rate=1,
				Inb_method="mom.visscher", threads=4, verbose=FALSE) {

	
	meta1 <- data.frame(metric=c("n references", "n samples", "n markers"), 
						value=c(length(input$ref.id), length(input$field.id), nrow(input$snp)))
	
	### recode SNP to be the number of A alleles
	# There are four possible values stored in the variable genotype: 0, 1, 2 and 3.
	# For bi-allelic SNP sites, 
	# 0 = two B alleles
	# 1 = one A and one B allele
	# 2 = two A alleles
	# 3 = missing 
	# For multi-allelic sites, it is a count of the reference allele (3 meaning no call). 
	# “Bit2” indicates that each byte encodes up to four SNP genotypes 
	# since one byte consists of eight bits.

	dmat <- input$snp
	dmat[] <- match(dmat, c(0,2,1,NA)) - 1
	
	## create a gds file
	SNPRelate::snpgdsCreateGeno(
					paste0(tempfile(), ".gds"), 
					genmat=dmat,
					sample.id = colnames(dmat), 
					snp.id = input$markers$MarkerName,
					snp.chromosome = input$markers$Chr,
					snp.position = input$markers$Pos,
					snp.allele = NULL, 
					snpfirstdim=TRUE)
	
		
	# open the gds file
	genofile <- SNPRelate::snpgdsOpen(input$gds)
	
	# calculate SNP maf and missing rate
	frq <- SNPRelate::snpgdsSNPRateFreq(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
	
	snpinfo <- data.frame(snpid=frq$snp.id, maf=frq$MinorFreq, mr=frq$MissingRate)
	d0 <- snpinfo[snpinfo$maf >= MAF_cutoff & snpinfo$mr <= SNP_Missing_Rate, ]
	
	# reformatting SNP info
	snpinfo$drop_maf <- (snpinfo$maf < MAF_cutoff) | is.na(snpinfo$maf)
	snpinfo$drop_mr <- snpinfo$mr > SNP_Missing_Rate
	snpinfo$drop <- snpinfo$drop_maf | snpinfo$drop_mr

	# calculate het rate
	h <- apply(dmat, 2, \(x) sum(x==1) / sum(x!=3))
	h <- data.frame(sid=colnames(dmat), h=h)


	clamp_merge_drop <- function(het, sid, inb, mr, mr_thold, het_thold, name) {
		x <- data.frame(sid=sid, inb=inb)
		mr <- data.frame(mr=mr)
		# normalize the range between 0 and 1
		x$inb[x$inb < 0] = 0
		x$inb[x$inb > 1] = 1
		x$het <- 1 - x$inb
		ic <- merge(x, het, by="sid")
		out <- merge(mr, ic, by.x="row.names", by.y="sid")
		out <- out[, c("Row.names", "mr", "h")]
		names(out)[c(1,3)] <- c(name, "het")
		out$drop_miss <- out$mr > mr_thold
		out$drop_het <- out$het > het_thold
		out$drop <- out$drop_miss | out$drop_het
		out
	}
 	
	# calculate sample missing rate
	if (length(input$field.id) > 0) {
		smp_mr <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=input$field.id, snp.id=d0$snpid, with.id=TRUE)
	# calculate inbreeding coefficient and then heterozygosity for field data
		inb_fld <- SNPRelate::snpgdsIndInb(genofile, sample.id=input$field.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN,missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)
		fld <- clamp_merge_drop(h, inb_fld$sample.id, inb_fld$inbreeding, smp_mr, Sample_Missing_Rate, Sample_Heterozygosity_Rate, "fld_id") 
	} else {
		smp_mr <- NA
		inb_fld <- NA
		fld <- data.frame(fld_id=0, het=0, drop_miss=TRUE)
	}
	ref_mr <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=input$ref.id, snp.id=d0$snpid, with.id=TRUE)
	# calculate inbreeding coefficient and then heterozygosity for ref data
	inb_ref <- SNPRelate::snpgdsIndInb(genofile, sample.id=input$ref.id, snp.id=d0$snp.id, autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)
	ref <- clamp_merge_drop(h, inb_ref$sample.id, inb_ref$inbreeding, ref_mr, Ref_Missing_Rate, Ref_Heterozygosity_Rate, "ref_id") 


	meta2 <- data.frame(
		metric=c("MAF_cutoff", "SNP_Missing_Raste Cutoff", "Marker Final", "Field Sample Total", "Field Sample Missing Cutoff", "Field Final", "Reference Sample Total", "Reference Sample Missing Cutoff", "Reference Final"), 
		value=c(MAF_cutoff, SNP_Missing_Rate, nrow(d0), length(input$field.id), Sample_Missing_Rate, sum(smp_mr <= Sample_Missing_Rate), length(input$ref.id), Ref_Missing_Rate, sum(ref_mr <= Sample_Missing_Rate))
	)

	meta3 <- data.frame(
		metric=c("Sample Het Avg", "Sample Het SD", "Sample Het Max", "Sample Het Min", "Ref Heterozygosity Cutoff", "Sample Heterozygosity Cutoff"), 
		value=c(mean(fld$het), 	stats::sd(fld$het), max(fld$het), min(fld$het), Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate)
	)
	meta <- rbind(meta1, meta2, meta3)
	
	# calculate pair-wise IBS

### RH instead of
#	ibs <- SNPRelate::snpgdsIBS(genofile, num.thread=threads, autosome.only=FALSE,
#		maf= MAF_cutoff, missing.rate=.5, verbose=verbose)

### RH omitting genotypes / markers with too much missing data 
	use_ids <- c(fld$fld_id[!fld$drop_miss], ref$ref_id[!ref$drop_miss])

	ibs <- SNPRelate::snpgdsIBS(genofile, sample.id=use_ids, num.thread=threads, autosome.only=FALSE, verbose=verbose,
		maf=MAF_cutoff, missing.rate=SNP_Missing_Rate)

	SNPRelate::snpgdsClose(genofile)

#	ref$variety <- genotypes$variety[match(gsub("_D.$", "", ref$ref_id), genotypes$sample)]
	ref$variety <- genotypes$variety[match(ref$ref_id, genotypes$ID)]

	list(metadata=meta, ref=ref, fld=fld, snpinfo = snpinfo, ibs=ibs, smp_mr=smp_mr, ref=ref)
}


match_IBS <- function(x, genotypes, markers, match_field, MAF_cutoff=0.05, SNP_Missing_Rate=0.2, 
				Ref_Missing_Rate=0.2, Sample_Missing_Rate=0.2, 
				Ref_Heterozygosity_Rate=1, Sample_Heterozygosity_Rate=1,
				IBS_cutoff=0.5, Inb_method="mom.visscher", assign_threshold=NULL,
				threads=4, verbose=FALSE, filename="") {


	input <- matchpoint:::prepare_data(x, genotypes, match_field, markers, filename=filename, verbose=verbose)
	
# x, genotypes, markers, 
	rIBS <- matchpoint:::run_IBS(input, MAF_cutoff=MAF_cutoff, SNP_Missing_Rate=SNP_Missing_Rate, 
				Ref_Missing_Rate=Ref_Missing_Rate, Sample_Missing_Rate=Sample_Missing_Rate, 
				Ref_Heterozygosity_Rate=Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate=Sample_Heterozygosity_Rate,
				Inb_method=Inb_method, threads=threads, verbose=verbose) 

	ibs <- rIBS$ibs
	smp_mr <- rIBS$smp_mr
	output <- rIBS[c("metadata", "ref", "fld", "snpinfo")]
	
	out_all <- ibs[[3]]
	colnames(out_all) <- rownames(out_all) <- ibs[[1]]
	i <- which(colnames(out_all) %in% input$ref.id)
	ref_match <- out_all[i, i]
	out_match <- out_all[-i, i]	

	if (is.null(assign_threshold)) {
		gtype <- x$geno$genotype[match(colnames(ref_match), x$geno$ID)]
		refnms <- genotypes$variety[match(gtype, genotypes$sample)]
		dimnames(ref_match) <- list(refnms, refnms)
		pun <- matchpoint:::punity(1-ref_match, seq(0, 0.5, .01))
		# we want the last which.max
		assign_threshold <- 1-pun[nrow(pun) - which.max(rev(pun[,"mean"])) + 1, "threshold"] 
		output$metadata <- rbind(output$metadata, data.frame(metric="assign_threshold (computed)", value=assign_threshold)) 
		output$punity <- data.frame(pun)
	} else {
		output$metadata <- rbind(output$metadata, data.frame(metric="assign_threshold (user input)", value=assign_threshold)) 	
	}

	
	d <- data.frame(GID=rownames(out_match), field_id=x$geno$genotype[match(rownames(out_match), x$geno$ID)], 
					ref_id=rep(x$geno$genotype[match(colnames(out_match), x$geno$ID)], each=nrow(out_match)),
					variety="", IBS=as.vector(out_match))

	out_match <- data.frame(FieldSample=colnames(out_all)[-i], out_match, check.names=FALSE)
	rownames(out_match) <- NULL

	i <- match(d$ref_id, genotypes$sample)
	if (any(is.na(i))) stop("check ID comparison")
	d$variety <- genotypes$variety[i]
	d$ref_Tid <- genotypes$TargetID[i]
#	mfld_id <- gsub("_D.$", "", d$field_id)
	mfld_id <- d$field_id
	i <- match(mfld_id, genotypes$sample)
	if (any(is.na(i))) stop("check ID comparison")
	d$field_Tid <- genotypes$TargetID[i]
	
	d <- d[order(d$field_id, -d$IBS), ]
# matches
	best <- d[!duplicated(d$field_id),]
	best <- merge(best, smp_mr, by.x="GID", by.y="row.names", all.y=TRUE)
	names(best)[ncol(best)] <- "sample_SNP_missing_rate"
	best$IBS <- round(best$IBS, 6)
	output[["best_match"]] <- best

	vg <- matchpoint:::old_var_groups(out_all, assign_threshold, input$ref.id)

	vg$variety <- sapply(strsplit(vg$ref_sample, ";"), 
				\(i) paste(unique(genotypes$variety[match(i, genotypes$sample)]), collapse="; ")
			)
	output[["all_match"]] <- vg

	ib <- d[d$IBS > IBS_cutoff[1], ] 
	ib$id_rank <- with(ib, stats::ave(IBS, field_id, FUN=\(x) rank(1-x, ties.method="min")))
	
	a <- stats::aggregate(ib[, "IBS", drop=FALSE], ib[, c("field_id", "variety")], max, na.rm=TRUE)
	a$var_rank <- with(a, stats::ave(IBS, field_id, FUN=\(x) rank(1-x, ties.method="min")))
	a$IBS <- NULL
	ib <- merge(ib, a, by=c("field_id",  "variety"))
	ib <- ib[order(ib$field_id, ib$id_rank), ]
	ib$IBS <- round(ib$IBS, 6)
	
	output[[paste0("IBS_id")]] <- ib
	
	dups <- duplicated(ib[, c("field_id", "var_rank")])
	output[[paste0("IBS_variety")]] <- ib[!dups, ]

	nr <- as.data.frame(table(field_id=d$field_id))
	rept <- stats::aggregate(d[, "IBS", drop=FALSE], d[, "field_id", drop=FALSE], 
			\(i) c(mean(i, na.rm=TRUE), stats::sd(i, na.rm=TRUE), length(stats::na.omit(i))))
	rept <- data.frame(rept[1], rept[[2]])
	colnames(rept)[-1] <- c("avg", "sd", "nr")

	meta4 <- data.frame(
		metric = paste0(c("Samples with match using IBS=", "Avg number matches per sample using IBS=", "Avg of IBS value with IBS=", "SD of IBS value with IBS="), IBS_cutoff), 
		value = c(nrow(rept), mean(rept$nr), mean(rept$avg, na.rm=TRUE), mean(rept$sd, na.rm=TRUE) ))

	output$metadata <- rbind(output$metadata, meta4)
	output$metadata$value <- round(output$metadata$value, 5)
	
	output[["similarity"]] <- out_match
	out_all <- data.frame(ID=colnames(out_all), 1-out_all, check.names=FALSE)
	rownames(out_all) <- NULL
	output[["distance"]] <- out_all


	if (input$filename != "") {
		xlsx <- paste0(input$filename, "_IBS.xlsx")
		dir.create(dirname(xlsx), FALSE, TRUE)
		writexl::write_xlsx(output, xlsx, format_headers=FALSE)
		saveRDS(ibs, paste0(input$filename, "_IBS.rds"))
		invisible(output)
	} else {
		output$gds <- input$gds
		output
	}
}


refine_IBS <- function(x, genotypes, match_field, markers, ref_split=0.1, ref_lump=0.05, 
				MAF_cutoff=0.05, SNP_Missing_Rate=0.2, 
				Ref_Missing_Rate=0.2, Sample_Missing_Rate=0.2, 
				Ref_Heterozygosity_Rate=1, Sample_Heterozygosity_Rate=1,
				Inb_method="mom.visscher",	threads=4, verbose=FALSE, filename="") {


	gtypes <- genotypes[genotypes$reference, ]
	input <- matchpoint:::prepare_data(x, gtypes, match_field, markers, filename=filename, verbose=verbose)
	
	rIBS <- matchpoint:::run_IBS(input, MAF_cutoff=MAF_cutoff, SNP_Missing_Rate=SNP_Missing_Rate, 
				Ref_Missing_Rate=Ref_Missing_Rate, Sample_Missing_Rate=Sample_Missing_Rate, 
				Ref_Heterozygosity_Rate=Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate=Sample_Heterozygosity_Rate,
				Inb_method=Inb_method, threads=threads, verbose=verbose) 

	ibs <- rIBS$ibs

	dstm <- 1-ibs[[3]]
	colnames(dstm) <- rownames(dstm) <- ibs[[1]]

	finish_refine(x, dstm, genotypes, match_field, ref_lump, ref_split, filename)

}

