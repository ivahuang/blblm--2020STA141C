#' @aliases NULL
#' @import purrr
#' @import stats
#' @import utils
#' @importFrom readr read_csv
#' @importFrom furrr future_map
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


#' @details
#' Regression models with Little Bag of Bootstraps
"_PACKAGE"
## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' global read_data function
#' @param folder directory
read_data = function(folder){
  file.path(folder, list.files(folder, pattern = "csv$")) %>%
    map(read_csv)
}


# function 1: linear regression
#' linear regression using LBB
#' @param formula regression formula
#' @param data data frame
#' @param m splitted data to m parts, default 10 splits
#' @param B numbers of bootstrap, default 5000
#' @export
blblm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  #regression
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  #residuals
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#function 2: linear regression wl parrael
#' linear regression using LLB, with parallelization
#' @param formula regression model
#' @param  data data frame
#' @param  m number of splits
#' @param  B number of bootstraps
#' @export

par_blblm <- function(formula, data, m = 10, B = 5000) {
  if(class(data) == "character"){
    data_list <- read_data(data)
  }
  else{
    data_list <- split_data(data, m)
  }
  
  #linear regression
  estimates <- future_map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(.), B = B))
  #store residuals
  res <- list(estimates = estimates, formula = formula)
  #assign class for further investigation
  class(res) <- "blblm"
  invisible(res)
}



# function 3: generalized linear regression
#' @param formula regression model
#' @param data data frame
#' @param m number of splits
#' @param B number of bootstrap
#' @param family some glm family
#' @export

blbglm <- function(formula, data, m = 10, B = 5000, family) {
  if(class(data) == "character"){
    data_list <- read_data(data)
  }
  else{
    data_list <- split_data(data, m)
  }
  # use glm here
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(.), B = B, family))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}

# function 4
#' generalized linear regression with parrel
#' @param formula regression model
#' @param data data frame
#' @param m number of splits
#' @param B number of bootstraps
#' @param family some glm family to use
#' @export

par_blbglm <- function(formula, data, m = 10, B = 5000, family) {
  if(class(data) == "character"){
    data_list <- read_data(data)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- future_map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(.), B = B, family))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#' @param data data frame
#' @param m number of splits
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

##########################the following three are about lm
#' lm for each subsample
#' @param formula regression model
#' @param data data frame
#' @param n how many vectors to use
#' @param B numbers of bootstrap
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute lm for each bootstrap
#' @param formula regression model
#' @param data data frame
#' @param n how many vectors to draw and use
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' lm for each bootstrap, specifying frequency
#' @param formula regression model
#' @param data data frame
#' @param freqs weights for each linear regressor
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

##########################the following three are about glm
#' lm for each subsample
#' @param formula regression model
#' @param data data frame
#' @param n how many vectors to use
#' @param B numbers of bootstrap
#' @param family glm family
glm_each_subsample <- function(formula, data, n, B, family) {
  replicate(B, glm_each_boot(formula, data, n,family), simplify = FALSE)
}


#' compute glm for each bootstrap
#' @param formula regression model
#' @param data data frame
#' @param n how many vectors to draw and use
#' @param family glm family
glm_each_boot <- function(formula, data, n, family) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, family)
}


#' glm for each bootstrap, specifying frequency
#' @param formula regression model
#' @param data data frame
#' @param freqs weights for each linear regressor
#' @param family glm family
glm1 <- function(formula, data, freqs, family) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, data,weights = freqs, family = family)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' compute the coefficients from fit
#' @param fit regression result fit
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit regression result fit
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
#' @param x regression result fit
#' @param ... other customized arguments
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' complute sigma for bootstrap regression fit
#' @export
#' @method sigma blblm
#' @param object LBB regression
#' @param confidence logical/boolean value
#' @param level overall intented confidence level for sigma's CI
#' @param ... other customized arguments
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}


#' coefficients for bootstrap lm
#' @export
#' @method coef blblm
#' @param object fit
#' @param ...  arguments
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' confidence interval for each terms
#' @export
#' @method confint blblm
#' @param object fit
#' @param parm boolean
#' @param level confidence level
#' @param ... arguments
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' @export
#' @method predict blblm
#' @param object fit
#' @param new_data data frame, list or environment
#' @param confidence boolean
#' @param level confidence level
#' @param ... customized arguments
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}

