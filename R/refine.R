

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