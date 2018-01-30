plotCoords <- function(dat, mst.tree) {

	dev.off()
	maps::map('world', col = "grey95", fill = TRUE, boundary = FALSE)
	points(dat$lngdec, dat$latdec, pch = 19, cex = 0.5, col = "#99000020")

	       pairwise.coords <- mst.tree[c("from","to")]

	       for (i in 1:nrow(mst.tree)) {

	       	x0 <- mean(subset(dat, collection_no == pairwise.coords[i,1])$lngdec)
	       	y0 <- mean(subset(dat, collection_no == pairwise.coords[i,1])$latdec)

	       	x1 <- mean(subset(dat, collection_no == pairwise.coords[i,2])$lngdec)
	       	y1 <- mean(subset(dat, collection_no == pairwise.coords[i,2])$latdec)

	       	segments(x0, y0, x1, y1, col = "blue")

	       }

}
