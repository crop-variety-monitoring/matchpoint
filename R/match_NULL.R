
get_clean_data_NULL <- function(snp, genotypes, filename) {
	snp <- matchpoint:::remove_unknown_samples(snp, genotypes$sample, verbose=FALSE)
	i <- match(colnames(snp)[-1], genotypes$sample)
#	snp <- fix_duplicate_names(snp, verbose=verbose)
	cns <- colnames(snp)[-1]
	ref.id <- cns[genotypes$reference[i]]
	field.id <- cns[!genotypes$reference[i]]

	snp2 <- t(snp[,-1])
	colnames(snp2) <- paste0("m", rep(1:(ncol(snp2)/2), each=2), rep(c("a", "b"), ncol(snp2)/2))
	snp2 <- data.frame(ID=row.names(snp2), snp2) 
	rownames(snp2) <- NULL
	i <- snp2$ID %in% ref.id
	train <- snp2[i, ]
	train$ID <- as.factor(train$ID)
	pred  <- snp2[!i, ]

	filename <- matchpoint:::fix_filename(filename)
	list(train=train, pred=pred, ref.id=ref.id, field.id=field.id, filename=filename)
}



match_NULL <- function(SNPs, genotypes, filename="") {

	input <- get_clean_data_RF(SNPs, genotypes, filename=filename)

	train = t(input$train[, -1]) |> as.matrix()
	i <- seq(1, nrow(train), 2)
	trr <- train[i, ] / (train[i,] + train[i+1, ])
	i <- colSums(is.na(trr)) > 0
	pred = input$pred[, -1] |> as.matrix()
	train = t(input$train[, -1]) |> as.matrix()
	i <- apply(pred, 1, \(x) which.min(colMeans(abs(train - x))))
	r <- colnames(train)[i]

	if (input$filename != "") {
		xlsx <- paste0(input$filename, "_RF.xlsx")
		writexl::write_xlsx(output, xlsx, format_headers=FALSE)
		saveRDS(ibs, paste0(input$filename, "_RF.rds"))
		invisible(output)
	} else {
		output
	}
}


