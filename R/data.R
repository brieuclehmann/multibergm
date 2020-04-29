#' Ten synthetic networks
#' 
#' A dataset of ten networks simulated from the same exponential random graph
#' model with different model parameters. 
#' 
#' @format A list of networks with 40 nodes.
#' 
"nets"

set.seed(1)
nNets <- 10
nNodes <- 20

ergmFormula <- ~ edges + gwesp(0.7, fixed = TRUE)
meanCoef <- c(-3, 1)
covCoef <- matrix(c(0.3, -0.1, -0.1, 0.1), 2, 2)

indivCoef <- mvtnorm::rmvnorm(10, meanCoef, covCoef)

nets <- list()

for (n in 1:nNets){
  nets[[n]] <- simulate(ergmFormula, coef = indivCoef[n, ],
                        basis = network(nNodes, directed=FALSE))
}
