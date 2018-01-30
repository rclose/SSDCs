outersect <- function(x, y) {
	sort(c(x[!x%in%y],
	       y[!y%in%x]))
}
