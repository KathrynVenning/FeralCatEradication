#' @export
max.lambda <- function(x) {
  Re((eigen(x)$values)[1])
}

return_one <- function() {
	return(1)
}
