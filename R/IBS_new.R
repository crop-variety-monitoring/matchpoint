
get_distances <- function(d, MAF_cutoff, SNP_mr, sample_mr, 
			IBS_cutoff, inb_method = "mom.visscher", cpus=1){
  
	inb_method <- match.arg(inb_method, c("mom.weir", "mom.visscher", "mle", "gcta1", "gcta2", "gcta3"))

##  create a gds file
	obj_gds <- paste0(tempfile(), "_geno.gds")
	SNPRelate::snpgdsCreateGeno(obj_gds, genmat = df,
                   sample.id = names(d)[-c(1:3)], snp.id = d$SNPID,
                   snp.chromosome = d$Chr,
                   snp.position = d$Pos,
                   snp.allele = NULL, snpfirstdim=TRUE)
 
  # open the gbs obj
	genofile <- SNPRelate::snpgdsOpen(obj_gds)
	on.exit(SNPRelate::snpgdsClose(genofile))

  # calculate SNP maf and missing rate
	frq <- SNPRelate::snpgdsSNPRateFreq(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)
	snpmr <- data.frame(snpid=frq$snp.id, maf=frq$MinorFreq, mr=frq$MissingRate)
	d0 <- snpmr[snpmr$maf >= MAF_cutoff & snpmr$mr <= SNP_mr, ]
	sm <- SNPRelate::snpgdsSampMissRate(genofile, sample.id=field.id, snp.id=d0$snpid, with.id=TRUE)

	meta123 <- describe_rf(d, d0, genofile, sm, MAF_cutoff, SNP_mr, sample_mr, inb_method)
    
  # calculate pair-wise IBS
	ibs <- SNPRelate::snpgdsIBS(genofile, num.thread=cpus, autosome.only=FALSE, maf= MAF_cutoff, missing.rate=SNP_mr)
	out <- ibs[[3]]
	rownames(out) <- colnames(out) <- ibs$sample.id

	xy <- t(utils::combn(colnames(out), 2))
  
  # convert to data.frame
	D <- data.frame(xy, dist=out[xy])
	names(D) <- c("ref", "field", "IBS")
  
  # subset ref vs. field
	D$type1 <- "field"
	D$type2 <- "ref"
	D[D$ref %in% ref.id, ]$type1 <- "ref"
	D[D$field %in% field.id, ]$type2 <- "field"
  
	D <- D[D$type1 == "ref" & D$type2 =="field", ]
  
	names(D)[1:3] <- c("FID1", "FID2", "IBS")
  
	outlist <- list()
	outlist[["metadata"]] <- meta123

	for(i in IBS_cutoff){
		lout <- get_match(D, cutoff=i)

		### best match
		lout2 <- merge(lout[[2]], sm, by.x="field_id", by.y="row.names", all.y=TRUE)
		names(lout2)[4] <- "Sample_SNP_mr"
		
		lout2$field_id <- gsub("^FLD_", "", lout2$field_id)
		lout2$ref_id <- gsub("^REF_", "", lout2$ref_id)

		name <- paste0("IBS_cutoff_", i, "_best_match")
		outlist[[name]] <- lout2
		
		name2 <- paste0("IBS_cutoff_", i, "_all_match")
		lout[[1]]$field_id <- gsub("^FLD_", "", lout[[1]]$field_id)
		lout[[1]]$ref_id <- gsub("^REF_", "", lout[[1]]$ref_id)
		outlist[[name2]] <- lout[[1]]

		
		rp <- stats::aggregate(lout[[1]]$IBS, lout[[1]][, "field_id", drop=FALSE], 
			function(i) c(mean(i, na.rm=TRUE), stats::sd(i, na.rm=TRUE)))
		rp <- apply(rp[,2], 2, mean, na.rm=TRUE)
								
		nr <- as.data.frame(table(lout[[1]]$field_id))
		
		meta4 <- data.frame(Metric=c(paste0("Samples with match using IBS=",i), paste0("Avg number matches per sample using IBS=",i), paste0("Avg of IBS value with IBS=", i), paste0("SD of IBS value with IBS=", i)), Value=c(nrow(nr), mean(nr$Freq), rp[1], rp[2]) )
		
		meta <- rbind(meta123, meta4)
	}
	outlist[["metadata"]] <- meta
	outlist[["snpinfo"]] <- snpmr
	outlist[["distance"]] <- out
	outlist
}



recode_dart <- function(d, biallelic=TRUE, missflags="-") {

	### recode SNP to be the number of A alleles
	# There are four possible values stored in the variable genotype: 0, 1, 2 and 3. 
	# For bi-allelic SNP sites, 
		#0 = two B alleles, 
		#1 = one A and one B allele, 
		#2 = two A alleles
		#3 = missing
		
	# For multi-allelic sites, it is a count of the reference allele (3 meaning no call).
	# “Bit2” indicates that each byte encodes up to four SNP genotypes 
	# since one byte consists of eight bits.

	d <- as.matrix(d)
	u <- unique(as.vector(tmp))
	unexpected <- u[!(u %in% c("0", "1", "2", missflags, NA))]
	if (length(unexpected) > 0) {
		stop(paste("Unexpected allelic value(s) found:", unexpected))
	}
	
	# missing
	tmp[is.na(tmp)] <- 3
	tmp[tmp %in% missflags] <- 3
	
	tmp <- apply(tmp, 2, as.numeric)
	if (biallelic) {
		tmp[tmp == 2] <- 9 # 2 het to 9 tmp
		tmp[tmp == 1] <- 2 # 1 homo alt to 2
		tmp[tmp == 9] <- 1 # 9 het to 1
	}
	
	tmp
}


new_workflow_rf  <- function(dart_file, ref_file, MAF_cutoff=0.05, SNP_mr=0.2, sample_mr=0.2, biallelic=TRUE, missflags="-", inb_method = "mom.visscher", cpus=1) {
  
	dartlst <- read_dart(dart_file)
	dart[, -1] <- recode(dartlst$snp[,-1], biallelic=biallelic, missflags=missflags)
	ref  <- data.table::fread(ref_file) 

	get_distances(dart, MAF_cutoff=MAF_cutoff, SNP_mr=SNP_mr, sample_mr=sample_mr, IBS_cutoff=IBS_cutoff, inb_method =inb_method, cpus=cpus)
	
   fxl <- file.path(outdir, paste0(c(country, crop, year, "matches.xlsx"), collapse="_"))
   writexl::write_xlsx(out, fxl)
   fxl
}
