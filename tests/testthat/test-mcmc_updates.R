### mean_update tests

test_that("Zero prior covariance returns prior mean", {
  prior_cov <- diag(0, 2)
  prior_mean <- c(1, 1)
  obs_cov <- matrix(c(1, 0.1, 0.1, 2), nrow = 2)
  obs_data <- mvtnorm::rmvnorm(4, c(1,1), obs_cov)
  
  expect_equal(mean_update(obs_data, obs_cov, prior_mean, prior_cov), 
               prior_mean)
})


test_that("Identity covariance returns expected solution", {
  n <- 10
  
  prior_cov <- diag(1, 2)
  prior_mean <- c(1, 1)
  obs_cov <- diag(1, 2)
  obs_data <- mvtnorm::rmvnorm(n, prior_mean, obs_cov)
  
  post_mean <- (1/(n+1))*(prior_mean + colSums(obs_data))
  
  set.seed(1)
  out <- mvtnorm::rmvnorm(1, post_mean, diag(1/(n + 1), 2))[1, ]
  
  set.seed(1)
  expect_equal(mean_update(obs_data, obs_cov, prior_mean, prior_cov), 
               out)
})