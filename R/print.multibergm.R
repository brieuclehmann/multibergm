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

  s_net <- ifelse(nNets > 1, "s", "")
  s_group <- ifelse(nGroups > 1, "s", "")
  cat("A multibergm fit of ", nNets, " network", s_net, " in ", 
      nGroups, " group", s_group, ".\n", sep = "")
  cat("\nModel formula:\n")
  print(x$formula, ...)
  cat("\nConstraints:\n")
  print(x$constraints, showEnv = FALSE, ...)
  
  invisible(x)
}
