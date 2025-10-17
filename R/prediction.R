#' Prediction
#'
#' @param fit A fitted `stan` model object returned from \code{dsgd()}.
#' @param data A dataframe of covariates and coordinates. The last two columns must be the spatial coordinates.
#' @param sweep The number of Gibbs sampling iterations.
#'
#' @return A vector of predicted binary outcomes.
#' @export
#'
#' @importFrom rstan summary
#' @examples \donttest{}


predict <- function(fit, data, sweep = 100) {
  
  n <- nrow(data)
  y_new <- sample(c(0,1), size = n, replace = TRUE)
  neighbor_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Extract coefficients
  coef_est <- rstan::summary(fit)$summary[, "mean"]
  P <- ncol(data) - 2
  beta <- coef_est[1 : P]
  eta <- coef_est[P + 1]
  
  
  for (s in 1 : sweep) {
    for (i in 1 : n) {
      for (j in 1 : n) {
        if (j != i && sqrt(((data[i,length(data)-1]-data[j,length(data)-1])^2)+((data[i,length(data)]-data[j,length(data)])^2))<=1) {
          neighbor_matrix[i, j] <- eta * y_new[j]
        } else {
          neighbor_matrix[i, j] <- 0
        }
      }
      logit_p <- sum(beta * data[i, 1:P]) + sum(neighbor_matrix[i, ])
      p <- plogis(logit_p)
      y_new[i] <- rbinom(1, size = 1, prob = p)
    }
  }
  
  return(y_new)
}


