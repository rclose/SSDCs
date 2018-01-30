#function to return "mean" value for range of data classes

my_mean <- function(x) {
	if (all(is.na(x))) {return(NA)}
	else if (class(x) == "numeric") {mean(x, na.rm = TRUE)}
	else if (class(x) == "integer") {round(mean(x, na.rm = TRUE))}
	else if (class(x) == "factor" | class(x) == "character") {names(sort(table(x),decreasing = TRUE)[1])}
}
