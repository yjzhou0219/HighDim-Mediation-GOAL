## Bandwidth Diagnostics




### ESS Decomposition and CV estimation

#' @param w_final Numeric vector of the final generalized propensity score (GPS) weights corresponding to a specific potential outcome component estimation (e.g., E[Y(1,M(1))]).
#' @param w_kernel Numeric vector of the kernel weights derived from the kernel smoothing function (e.g., Epanechnikov kernel) for the target exposure evaluation point.
#' @param eps A small numeric scalar used as a tolerance threshold to handle floating-point inaccuracies and prevent division by zero. Default is \code{1e-12}.
#' @details This function performs a pre-estimation diagnostic evaluation of local empirical support and weight stability. It decomposes the effective sample size (ESS) into two components: the kernel ESS (reflecting local data density determined by the kernel window) and the final weighted ESS (further penalized by the GPS weights). It also computes the coefficient of variation (CV) to assess the dispersion and stability of the final weights, both globally and within the active kernel window.
#' @return A named numeric vector containing diagnostic metrics: \code{ESS_W} (final weighted ESS), \code{ESS_K} (kernel ESS), \code{ESS_ratio} (ratio of final ESS to kernel ESS, indicating additional information loss due to weighting), \code{CV_global} (global coefficient of variation for final weights), \code{CV_active} (coefficient of variation within the active kernel window), \code{n_active_final} (number of observations with non-zero final weights), and \code{n_active_kernel} (number of observations with non-zero kernel weights).


calc_ess_decomp <- function(w_final, w_kernel, eps = 1e-12) {
  w_final <- as.numeric(w_final)
  w_kernel <- as.numeric(w_kernel)
  
  w_final[!is.finite(w_final)] <- 0
  w_kernel[!is.finite(w_kernel)] <- 0
  
  if (sum(w_final) <= eps || sum(w_kernel) <= eps) {
    return(c(
      ESS_W = NA_real_,
      ESS_K = NA_real_,
      ESS_ratio = NA_real_,
      CV_w_global = NA_real_,
      CV_active = NA_real_,
      n_active_final = sum(w_final > 0),
      n_active_kernel = sum(w_kernel > 0)
    ))
  }
  
  w_final <- w_final / sum(w_final)
  w_kernel <- w_kernel / sum(w_kernel)
  
  ESS_W <- 1 / sum(w_final^2, na.rm = TRUE)
  ESS_K <- 1 / sum(w_kernel^2, na.rm = TRUE)
  ESS_ratio <- ESS_W / ESS_K
  
  CV_w_global <- sqrt(mean((w_final - mean(w_final))^2)) / mean(w_final)
  
  active_final <- w_final > 0
  
  if (sum(active_final) > 1 && mean(w_final[active_final]) > eps) {
    CV_active <- sd(w_final[active_final]) / mean(w_final[active_final])
  } else {
    CV_active <- NA_real_
  }
  
  c(
    ESS_W = ESS_W,
    ESS_K = ESS_K,
    ESS_ratio = ESS_ratio,
    CV_global = CV_w_global,
    CV_active = CV_active,
    n_active_final = as.integer(sum(active_final)),
    n_active_kernel = as.integer(sum(w_kernel > 0))
  )
}






### SMD metrics

#' @param x_trimmed A data frame or matrix of pre-treatment covariates that have been trimmed to remove observations with extreme generalized propensity scores. The covariate set should correspond to the specific potential outcome components being evaluated.
#' @param w Numeric vector of the relative weights corresponding to a specific potential outcome state (e.g., \code{relative_wgt_a1m1}). Length must match the number of rows in \code{x_trimmed}.
#' @param eps A small numeric scalar added to the denominator to prevent division by zero when unweighted standard deviations are extremely small. Default is \code{1e-10}.
#' @details This function assesses covariate balance across the exposure evaluation spectrum by calculating the Standardized Mean Difference (SMD). It compares the weighted mean of each covariate (representing the pseudo-population constructed by GPS weighting) to its unweighted mean (representing the original trimmed sample), standardized by the unweighted standard deviation.
#' @return A numeric vector of SMDs corresponding to each covariate column in \code{x_trimmed}.


calc_smd <- function(x_trimmed, w, eps = 1e-10) {
  x_trimmed <- as.data.frame(x_trimmed)
  w <- as.numeric(w)
  
  if (nrow(x_trimmed) != length(w)) {
    stop("Length of weights does not match number of rows in x_trimmed.")
  }
  
  if (abs(sum(w, na.rm = TRUE) - 1) > 1e-6) {
    warning("Weights do not sum to 1; normalizing weights.")
    w <- w / sum(w, na.rm = TRUE)
  }
  
  mean_unwt <- colMeans(x_trimmed, na.rm = TRUE)
  sd_unwt <- apply(x_trimmed, 2, sd, na.rm = TRUE)
  mean_wt <- colSums(sweep(as.matrix(x_trimmed), 1, w, `*`), na.rm = TRUE)
  
  smd <- abs(mean_wt - mean_unwt) / (sd_unwt + eps)
  smd[sd_unwt < eps] <- NA_real_
  return(smd)
}






