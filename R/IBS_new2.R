
dart_IBS <- function(dart_file, info_file, MAF_cutoff, SNP_Missing_Rate, 
				Ref_Missing_Rate, Sample_Missing_Rate, 
				Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate,
				IBS_cutoff, outdir, out_prefix, Inb_method = "mom.visscher", cpus=1,
				verbose=FALSE){


	Inb_method <- match.arg(Inb_method, c("mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"))
 
	dart <- data.frame(data.table::fread(dart_file), check.names=FALSE)
	info <- data.table::fread(info_file)

	stopifnot(length(unique(colnames(dart))) == ncol(dart))
	cns <- colnames(dart)[-c(1:3)]
	mcns <- gsub("_D1$|_D2$", "", cns)
	i <- match(mcns, info$dart.id)
	if (any(is.na(i))) {
		stop("cannot match all dart.id")
	}
	ir <- info$type[i] == "reference"
	ref.id <- cns[which(ir)]
	field.id <- cns[which(!ir)]
#	hdr.id <- colnames(dart)[1:3]
#	ref <- dart[ , c(hdr.id, ref.id)]
#	field <- dart[ , c(hdr.id, fld.id)]
	
	meta1 <- data.frame(Metric=c("Marker Total Ref", "Marker Total Sample",	"Marker Common"), 
						Value=c(length(ref.id), length(field.id), nrow(dart)))
	
	### recode SNP to be the number of A alleles
	# There are four possible values stored in the variable genotype: 0, 1, 2 and 3. 
	# For bi-allelic SNP sites, “0” indicates two B alleles, “1” indicates one A allele and one B allele, “2” indicates 
	# two A alleles, and “3” is a missing genotype. 
	# For multi-allelic sites, it is a count of the reference allele (3 meaning no call). 
	# “Bit2” indicates that each byte encodes up to four SNP genotypes since one byte consists of eight bits.
	tem <- as.matrix(dart[, -c(1:3)])
	tem[] <- match(tem, c(0,2,1,NA)) - 1

	df <- apply(tem, 2, as.numeric)
	df0 <- cbind(dart[, 1:3], tem)
	
	# calculate het rate
	h <- apply(df, 2, matchpoint:::het_rate)
	
	dir.create(outdir, FALSE, TRUE)
	
	obj_gds <- file.path(outdir, paste0(out_prefix, "_geno.gds"))
	## create a gds file
	SNPRelate::snpgdsCreateGeno(obj_gds, genmat = df,
					sample.id = names(df0)[-c(1:3)], snp.id = df0$MarkerName,
					snp.chromosome = df0$Chr,
					snp.position = df0$Pos,
					snp.allele = NULL, snpfirstdim=TRUE)
	
	
	# open the gbs obj
	genofile <- SNPRelate::snpgdsOpen(obj_gds)
	
	# calculate SNP maf and missing rate
	frq <- SNPRelate::snpgdsSNPRateFreq(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
	d <- data.frame(snpid=frq$snp.id, maf=frq$MinorFreq, mr=frq$MissingRate)
	d0 <- subset(d, maf >= MAF_cutoff & mr <= SNP_Missing_Rate)
	
	# reformatting SNP info
	snpinfo <- d
	snpinfo$drp_maf <- snpinfo$drp_mr <- snpinfo$drp <- 0
	snpinfo$drp_maf[(snpinfo$maf < MAF_cutoff) | is.na(snpinfo$maf)] <- 1
	snpinfo$drp_mr[snpinfo$mr > SNP_Missing_Rate] <- 1
	snpinfo$drp[snpinfo$drp_maf | snpinfo$drp_mr] <- 1
 	
	# calculate sample missing rate
	sm <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=field.id, snp.id=d0$snpid, with.id=TRUE)
	fld <- as.data.frame(sm)
	#s0 <- subset(fld, sm <= Sample_Missing_Rate)
	
	rm <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=ref.id, snp.id=d0$snpid, with.id=TRUE)
	rf <- as.data.frame(rm)
	#r0 <- subset(rf, rm <= Ref_Missing_Rate)
	meta2 <- data.frame(
		Metric=c("Marker MAF Cutoff", "Marker Missing Cutoff", "Marker Final", "Marker Final Coverage", "Field Sample Total", "Field Sample Missing Cutoff", "Field Final", "Reference Sample Total", "Reference Sample Missing Cutoff", "Reference Final"), 
		Value=c(MAF_cutoff, SNP_Missing_Rate, nrow(d0), ".", length(field.id), Sample_Missing_Rate, nrow(subset(fld, sm <= Sample_Missing_Rate)), length(ref.id), Ref_Missing_Rate, nrow(subset(fld, sm <= Sample_Missing_Rate)))
	)
	
	# calculate inbreeding coefficient and then heterozygosity for field data
	inb <- SNPRelate::snpgdsIndInb(genofile, sample.id=field.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN,missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)
	
	# normalize the range between 0 and 1
	ic <- data.frame(sid=inb$sample.id, inb=inb$inbreeding)
	ic$inb[ic$inb < 0] = 0
	ic$inb[ic$inb >1] = 1
	ic$het <- 1 - ic$inb

	#range(ic$inb)
	h1 <- as.data.frame(h)
	h1$sid <- row.names(h1)
	ic <- merge(ic, h1, by="sid")
	fld <- merge(fld, ic, by.x="row.names", by.y="sid")
	fld <- fld[, c("Row.names", "sm", "h")]
	names(fld) <- c("fld_smp_id", "mr", "het")
	fld$drp <- fld$drp_miss <- fld$drp_het <- 0
	fld$drp_miss[fld$mr > Sample_Missing_Rate] <- 1
	## should this be "fld$het < Sample_Heterozygosity_Rate" ?
	## either way, it seems this is not used
	fld$drp_miss[fld$het > Sample_Heterozygosity_Rate] <- 1
	fld$drp[fld$drp_miss | fld$drp_het] <- 1
	
	# calculate inbreeding coefficient and then heterozygosity for ref data
	inb_ref <- SNPRelate::snpgdsIndInb(genofile, sample.id=ref.id, snp.id=d0$snp.id, autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE,verbose=verbose)
	
	# normalize the range between 0 and 1
	ic1 <- data.frame(sid=inb_ref$sample.id, inb=inb_ref$inbreeding)
	ic1$inb[ic1$inb < 0] <- 0
	ic1$inb[ic1$inb > 1] <- 1
	ic1$het <- 1 - ic1$inb
	#range(ic$inb)
	h1 <- as.data.frame(h)
	h1$sid <- row.names(h1)
	ic1 <- merge(ic1, h1, by="sid")
	rf <- merge(rf, ic1, by.x="row.names", by.y="sid")
	rf <- rf[, c("Row.names", "rm", "h")]
	names(rf) <- c("ref_smp_id", "mr", "het")
	rf$drp_miss <- rf$drp_het <- rf$drp <- 0
	rf$drp_miss[rf$mr > Ref_Missing_Rate] <- 1
	## should this be "fld$het < Sample_Heterozygosity_Rate" ?
	## either way, it seems this is not used
	rf$drp_het[rf$het > Ref_Heterozygosity_Rate] <- 1
	rf$drp[rf$drp_miss | rf$drp_het] <- 1

	meta3 <- data.frame(
		Metric=c("Sample Het Avg", "Sample Het SD", "Sample Het Max", "Sample Het Min", "Ref Heterozygosity Cutoff", "Sample Heterozygosity Cutoff"), 
		Value=c(mean(ic1$h), sd(ic1$h), max(ic1$h), min(ic1$h),								Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate)
	)
	meta <- rbind(meta1, meta2, meta3)
	outlist <- list(metadata=meta)
	
	# calculate pair-wise IBS
	ibs <- SNPRelate::snpgdsIBS(genofile, num.thread=cpus, autosome.only=FALSE, maf= MAF_cutoff, missing.rate=SNP_Missing_Rate, verbose=verbose)
	SNPRelate::snpgdsClose(genofile)
	
	out <- round(ibs[[3]], 4)
	colnames(out) <- rownames(out) <- names(df0)[-c(1:3)]
	i <- colnames(out) %in% ref.id
	out2 <- out[-i, i]
	
	d <- data.frame(field_id=rownames(out2), 
					ref_id=rep(colnames(out2), each=nrow(out2)),
					variety="", IBS=as.vector(out2))

	out2 <- data.frame(FieldSample=colnames(out)[-i], out2, check.names=FALSE)
	rownames(out2) <- NULL

	mref_id <- d$ref_id
	mref_id <- gsub("_D1$|_D2$", "", mref_id)

	d$variety <- info$variety[match(mref_id, info$dart.id)]
	d <- d[with(d, order(field_id, -IBS)),]

# matches
	best <- d[!duplicated(d$field_id),]
	best <- merge(best, sm, by.x="field_id", by.y="row.names", all.y=TRUE)
	names(best)[5] <- "Sample_SNP_Missing_Rate"
	outlist[["best_match"]] <- best

	IBS_cutoff <- IBS_cutoff[1]
	d <- d[d$IBS > IBS_cutoff, ]
	outlist[[paste0("IBS_", IBS_cutoff)]] <- d
		
	report <- plyr::ddply(d, plyr::`.`(field_id), plyr::summarise,
							avg = mean(IBS, na.rm=TRUE),
							sd = sd(IBS, na.rm=TRUE))
	nr <- plyr::ddply(d, plyr::`.`(field_id), nrow)


	meta4 <- data.frame(
		Metric=c(paste0("Samples with match using IBS=",IBS_cutoff), paste0("Avg number matches per sample using IBS=",IBS_cutoff), paste0("Avg of IBS value with IBS=", IBS_cutoff), paste0("SD of IBS value with IBS=", IBS_cutoff)), 
		Value=c(nrow(report), mean(nr$V1), mean(report$avg, na.rm = TRUE), mean(report$sd, na.rm = TRUE) ))
	meta <- rbind(meta, meta4)

	outlist[["similarity"]] <- out2
	out <- data.frame(ID=colnames(out), 1-out, check.names=FALSE)
	rownames(out) <- NULL
	outlist[["distance"]] <- out
	outlist[["snpinfo"]] <- snpinfo
	colnames(rf)[1] <- "ref_id"
	outlist[["ref"]] <- rf
	colnames(fld)[1] <- "field_id"
	outlist[["field"]] <- fld

	#file.remove(obj_gds)
	xlsx <- file.path(outdir, paste0(out_prefix, ".xlsx"))
	writexl::write_xlsx(outlist, xlsx)
	invisible(outlist)
}


