
ref_field_IBS2 <- function(ref_file, field_file, MAF_cutoff, SNP_Missing_Rate, 
				Ref_Missing_Rate, Sample_Missing_Rate, 
				Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate,
				IBS_cutoff, outdir, out_prefix, Inb_method = "mom.visscher", cpus=1,
				verbose=FALSE){
	
	Inb_method <- match.arg(Inb_method, c("mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"))
 
	ref <- data.table::fread(ref_file)
	field <- data.table::fread(field_file)
	stopifnot(nrow(ref) == nrow(field))
	
	# assign UID to ref samples
# should be fixed in input
	stopifnot(length(unique(colnames(ref))) == ncol(ref))
	stopifnot(length(unique(colnames(field))) == ncol(field))
	
	ref.id = paste0("ref",1:(ncol(ref)-3), "_", names(ref)[-c(1:3)])
	names(ref) <- c("SNPID", "Chr", "Pos", ref.id)
	
	# assign UID to field samples
	field.id = paste0("field",1:(ncol(field)-3), "_", names(field)[-c(1:3)])
	names(field) <- c("SNPID", "Chr", "Pos", field.id)
	
	df9 <- merge(ref, field[, -c(2:3)], by="SNPID")
	
	meta1 <- data.frame(Metric=c("Marker Total Ref", "Marker Total Sample",	"Marker Common"), 
							Value=c(nrow(ref), nrow(field), nrow(df9)))
	
	### recode SNP to be the number of A alleles
	# There are four possible values stored in the variable genotype: 0, 1, 2 and 3. 
	# For bi-allelic SNP sites, “0” indicates two B alleles, “1” indicates one A allele and one B allele, “2” indicates 
	# two A alleles, and “3” is a missing genotype. 
	# For multi-allelic sites, it is a count of the reference allele (3 meaning no call). 
	# “Bit2” indicates that each byte encodes up to four SNP genotypes since one byte consists of eight bits.
	tem <- as.matrix(df9[, -c(1:3)])
	tem[] <- match(tem, c(0,2,1,NA)) - 1

	df <- apply(tem, 2, as.numeric)
	df0 <- cbind(df9[, 1:3], tem)
	
	# calculate het rate
	h <- apply(df, 2, matchpoint:::het_rate)
	
	dir.create(outdir, FALSE, TRUE)
	
	obj_gds <- file.path(outdir, paste0(out_prefix, "_geno.gds"))
	## create a gds file
	snpgdsCreateGeno(obj_gds, genmat = df,
					sample.id = names(df0)[-c(1:3)], snp.id = df0$SNPID,
					snp.chromosome = df0$Chr,
					snp.position = df0$Pos,
					snp.allele = NULL, snpfirstdim=TRUE)
	
	
	# open the gbs obj
	genofile <- snpgdsOpen(obj_gds)
	
	# calculate SNP maf and missing rate
	frq <- snpgdsSNPRateFreq(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
	d <- data.frame(snpid=frq$snp.id, maf=frq$MinorFreq, mr=frq$MissingRate)
	d0 <- subset(d, maf >= MAF_cutoff & mr <= SNP_Missing_Rate)
	
	# reformatting SNP info
	snpinfo <- d
	snpinfo$drp_maf <- 0
	snpinfo$drp_mr <- 0
	snpinfo$drp <- 0
	snpinfo$drp_maf[(snpinfo$maf < MAF_cutoff) | is.na(snpinfo$maf)] <- 1
	snpinfo$drp_mr[snpinfo$mr > SNP_Missing_Rate] <- 1
	snpinfo$drp[snpinfo$drp_maf | snpinfo$drp_mr] <- 1
 	
	# calculate sample missing rate
	sm <- snpgdsSampMissRate(genofile, sample.id=field.id, snp.id=d0$snpid, with.id=TRUE)
	fld <- as.data.frame(sm)
	#s0 <- subset(fld, sm <= Sample_Missing_Rate)
	
	rm <- snpgdsSampMissRate(genofile, sample.id=ref.id, snp.id=d0$snpid, with.id=TRUE)
	rf <- as.data.frame(rm)
	#r0 <- subset(rf, rm <= Ref_Missing_Rate)
	meta2 <- data.frame(
		Metric=c("Marker MAF Cutoff", "Marker Missing Cutoff", "Marker Final", "Marker Final Coverage", "Field Sample Total", "Field Sample Missing Cutoff", "Field Final", "Reference Sample Total", "Reference Sample Missing Cutoff", "Reference Final"), 
		Value=c(MAF_cutoff, SNP_Missing_Rate, nrow(d0), ".", ncol(field)-3, Sample_Missing_Rate, nrow(subset(fld, sm <= Sample_Missing_Rate)), ncol(ref)-3, Ref_Missing_Rate, nrow(subset(fld, sm <= Sample_Missing_Rate)))
	)
	
	# calculate inbreeding coefficient and then heterozygosity for field data
	inb <- snpgdsIndInb(genofile, sample.id=field.id, snp.id=d0$snp.id, autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN,missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE, verbose=verbose)
	
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
	fld$drp_miss <- 0
	fld$drp_het <- 0
	fld$drp_miss[fld$mr > Sample_Missing_Rate] <- 1
	## ? should be "fld$het < Sample_Heterozygosity_Rate"
	## either way, it seems this is not used
	fld$drp_miss[fld$het > Sample_Heterozygosity_Rate] <- 1
	fld$drp[fld$drp_miss | fld$drp_het] <- 1
	
	# calculate inbreeding coefficient and then heterozygosity for ref data
	inb_ref <- snpgdsIndInb(genofile, sample.id=ref.id, snp.id=d0$snp.id, autosome.only=FALSE, remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, method=Inb_method, allele.freq=NULL, out.num.iter=TRUE,verbose=verbose)
	
	# normalize the range between 0 and 1
	ic1 <- data.frame(sid=inb_ref$sample.id, inb=inb_ref$inbreeding)
	ic1$inb[ic1$inb < 0] = 0
	ic1$inb[ic1$inb >1] = 1
	ic1$het <- 1 - ic1$inb
	#range(ic$inb)
	h1 <- as.data.frame(h)
	h1$sid <- row.names(h1)
	ic1 <- merge(ic1, h1, by="sid")
	rf <- merge(rf, ic1, by.x="row.names", by.y="sid")
	rf <- rf[, c("Row.names", "rm", "h")]
	names(rf) <- c("ref_smp_id", "mr", "het")
	rf$drp_miss <- 0
	rf$drp_het <-	0
	rf$drp <-	0
	rf$drp_miss[rf$mr > Ref_Missing_Rate] <- 1
	## ? should be "fld$het < Sample_Heterozygosity_Rate"
	## either way, it seems this is not used
	rf$drp_het[rf$het > Ref_Heterozygosity_Rate] <- 1
	rf$drp[rf$drp_miss | rf$drp_het] <- 1

	meta3 <- data.frame(
		Metric=c("Sample Het Avg", "Sample Het SD", "Sample Het Max", "Sample Het Min", "Ref Heterozygosity Cutoff", "Sample Heterozygosity Cutoff"), 
		Value=c(mean(ic1$h), sd(ic1$h), max(ic1$h), min(ic1$h),								Ref_Heterozygosity_Rate, Sample_Heterozygosity_Rate)
	)
	
	# calculate pair-wise IBS
	ibs <- snpgdsIBS(genofile, num.thread=cpus, autosome.only=FALSE, maf= MAF_cutoff, missing.rate=SNP_Missing_Rate, verbose=verbose)

	snpgdsClose(genofile)

	
	out <- ibs[[3]]
	dimnames(out) <- dimnames(out) <- list(names(df0)[-c(1:3)], names(df0)[-c(1:3)])
	xy <- t(combn(colnames(out), 2))
	
	# convert to data.frame
	d <- data.frame(xy, dist=out[xy])
	names(d) <- c("ref", "field", "IBS")
	
	# subset ref vs. field
	d$type1 <- "field"
	d$type2 <- "ref"
	d[d$ref %in% names(ref)[-c(1:3)], ]$type1 <- "ref"
	d[d$field %in% names(field)[-c(1:3)], ]$type2 <- "field"
	
	d <- subset(d, type1 == "ref" & type2 =="field")
	
	names(d)[1:3] <- c("FID1", "FID2", "IBS")
	
	meta <- rbind(meta1, meta2, meta3)
	outlist <- list()

# matches
	dmatch <- d[, c("FID2", "FID1", "IBS")]
	names(dmatch) <- c("field_id", "ref_id", "IBS")
	dmatch <- dmatch[with(dmatch, order(field_id, -IBS)),]
## restoring names
	dmatch$ref_id <- matchpoint:::remove_before(dmatch$ref_id)
	dmatch$field_id <- matchpoint:::remove_before(dmatch$field_id)

	best <- dmatch[!duplicated(dmatch$field_id),]
	best <- merge(best, sm, by.x="field_id", by.y="row.names", all.y=T)
	names(best)[4] <- "Sample_SNP_Missing_Rate"
	outlist[["best_match"]] <- best

	IBS_cutoff <- IBS_cutoff[1]
	dmatch <- dmatch[dmatch$IBS > IBS_cutoff, ]
	outlist[[paste0("IBS_", IBS_cutoff)]] <- dmatch
		
	report <- plyr::ddply(dmatch, plyr::`.`(field_id), plyr::summarise,
							avg = mean(IBS, na.rm=TRUE),
							sd = sd(IBS, na.rm=TRUE))
	nr <- plyr::ddply(dmatch, plyr::`.`(field_id), nrow)


	meta4 <- data.frame(Metric=c(paste0("Samples with match using IBS=",IBS_cutoff), paste0("Avg number matches per sample using IBS=",IBS_cutoff), paste0("Avg of IBS value with IBS=", IBS_cutoff), paste0("SD of IBS value with IBS=", IBS_cutoff)), 
	Value=c(nrow(report), mean(nr$V1), mean(report$avg, na.rm = TRUE), mean(report$sd, na.rm = TRUE) ))
	meta <- rbind(meta, meta4)

	out <- round(out, 4)
	i <- grep("^ref", colnames(out))
	rownames(out) <- colnames(out) <- matchpoint:::remove_before(colnames(out))

	out2 <- out[-i, i]
	out2 <- data.frame(FieldSample=colnames(out)[-i], out2, check.names=FALSE)
	rownames(out2) <- NULL
	outlist[["similarity"]] <- out2

	out <- data.frame(1-out, check.names=FALSE)
	rownames(out) <- NULL
	outlist[["distance"]] <- out

	outlist[["metadata"]] <- meta
	outlist[["snpinfo"]] <- snpinfo

	
	colnames(rf)[1] <- "ref_id"
	rf$ref_id <- matchpoint:::remove_before(rf$ref_id)
	outlist[["ref"]] <- rf

	colnames(fld)[1] <- "field_id"
	fld$field_id <- matchpoint:::remove_before(fld$field_id)
	outlist[["field"]] <- fld


	#file.remove(obj_gds)
	xlsx <- file.path(outdir, paste0(out_prefix, "_matches.xlsx"))
	writexl::write_xlsx(outlist, xlsx)
	invisible(outlist)
}


