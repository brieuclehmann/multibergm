#' S3 print method for a multibergm object
#'
#' Automatically called when an object of class \code{\link{multibergm}} is
#' printed. Currently, summarizes the number of networks, the number of groups,
#' the model formula and the model constraints.
#'
#' @param x A \code{\link{multibergm}} object
#' @param ... Further arguments passed to other methods
#'
#' @export
print.multibergm <- function(x, ...) {
  networks <- statnet.common::eval_lhs.formula(x$formula)

  n_groups <- length(unique(x$groups))
  n_nets   <- length(networks)

  s_net <- ifelse(n_nets > 1, "s", "")
  s_group <- ifelse(n_groups > 1, "s", "")
  cat("A multibergm fit of ", n_nets, " network", s_net, " in ",
      n_groups, " group", s_group, ".\n", sep = "")
  cat("\nModel formula:\n")
  print(x$formula, showEnv = FALSE, ...)
  cat("\nConstraints:\n")
  print(x$constraints, showEnv = FALSE, ...)

  invisible(x)
}
