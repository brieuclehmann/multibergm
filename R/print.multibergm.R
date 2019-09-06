#' S3 print method for a multibergm object
#'
#' Automatically called when an object of class \code{\link{multibergm}} is printed.
#' Currently, summarizes the number of networks, the number of groups, the model
#' formula and the model constraints
#'
#' @param x A \code{\link{multibergm}} object
#' @param ... Further arguments passed to other methods
#'
#' @export
print.multibergm <- function(x, ...){
  nGroups <- length(unique(x$control$groupLabel))
  nNets   <- length(x$networks)

  cat("A multibergm fit of",nNets,"networks in",nGroups,"group(s).\n")
  cat("\nModel formula:\n")
  print(x$formula, ...)
  cat("\nConstraints:\n")
  print(x$constraints, ...)
}
