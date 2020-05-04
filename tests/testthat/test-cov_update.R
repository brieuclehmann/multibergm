
test_that("Zero prior scale returns current covariance", {
  prior_scale <- diag(0, 2)
  prior_df <- 1
  obs_mean <- c(1, 2)
  obs_data <- mvtnorm::rmvnorm(4, obs_mean)
  curr_cov <- diag(2)
  
  expect_equal(cov_update(obs_data, 
                          prior_df,
                          prior_scale,
                          obs_mean,
                          curr_cov = curr_cov), 
               curr_cov)
})


test_that("Identical observations returns expected solution", {
  n <- 10
  
  prior_scale <- diag(2)
  prior_df <- 1
  obs_mean <- c(1, 2)
  obs_data <- matrix(obs_mean, n, 2, byrow = TRUE)
  
  set.seed(1)
  out <- MCMCpack::riwish(n + prior_df, prior_scale)
  
  set.seed(1)
  expect_equal(cov_update(obs_data, 
                          prior_df,
                          prior_scale,
                          obs_mean), 
               out)
})