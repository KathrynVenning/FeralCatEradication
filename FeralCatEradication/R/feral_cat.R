#' @export
max.lambda <- function(x) {
  Re((eigen(x)$values)[1])
}

max.r <- function(x) {
  log(max.lambda(x))
}

return_one <- function() {
	return(1)
}
