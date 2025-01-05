#' Checks
#' @noRd
check_inputs <- function(n, d, sd, threshold, fraction, power, BFtype) {
  if (is.null(n) & is.null(power) || !is.null(n) & !is.null(power)) {
    stop("Either the target power or the sample size must be specified, but not both.")
  } else {
    if (!is.null(n)) {
      if (n < 1) {
        stop("The sample size must be at least 1.")
      }
    } else {
      if (power <= 0 | power >= 1) {
        stop("The target power must be between 0 and 1.")
      }
    }
  }
  if (!is.numeric(d)) {
    stop("The effect size must be a numeric value.")
  }
  if (!is.numeric(sd)) {
    stop("The standard deviation must be a numeric value.")
  } else {
    if (sd <= 0) {
      stop("The standard deviation must be greater than 0.")
    }
  }
  if (!is.numeric(threshold)) {
    stop("The threshold must be a numeric value.")
  } else {
    if (threshold <= 0) {
      stop("The threshold must be greater than 0.")
    }
  }
  if(!is.numeric(fraction)) {
    stop("The fraction must be a numeric value.")
  } else {
    if (fraction < 1) {
      stop("The fraction of information in the data used to construct the prior distribution must be at least 1. The default value 1 denotes the minimal fraction, 2 denotes twice the minimal fraction, and so on.")
    }
  }
  if (!BFtype %in% c("AAFBF", "AFBF")) {
    stop("The BFtype must be either 'AAFBF' or 'AFBF'.")
  }
}

check_hypothesis <- function(hypmat) {

  # Check how many hypotheses should be evaluated against each other
  nhypos <- length(hypmat$hyp_mat)
  lapply(hypmat$hypmat, \(x) {
    if (nrow(x) > 2) {
      stop("The hypotheses of interest cannot have more than two range constraints.")
    }
  })
  if (sum(hypmat$n_ec) > 1) {
    stop("You can only evaluate one equality-constrained hypothesis at the time.")
  }

  if (all(hypmat$n_ec == 1)) {
    hyptype <- "eq"
  } else if (all(hypmat$n_ec == 0) & hypmat$n_constraints[2] == 1) {
    hyptype <- "one.sided"
  } else if (max(hypmat$n_ec) == 1) {
    hyptype <- "eq.one.sided"
  } else if (all(hypmat$n_ec == 0) & hypmat$n_constraints[2] == 2) {
    hyptype <- "range"
  } else {
    hyptype <- "mixed"
  }

  hypval <- sapply(hypmat$hyp_mat, \(x) {
    if (nrow(x) == 1) x[1,2]
    else if (nrow(x) == 2) mean(x[,1]*x[,2])
  })
  if (length(unique(hypval)) > 1) {
    stop("The hypotheses of interest must share the same prior center.")
  }
  hypmat$hyptype <- hyptype
  hypmat$hypval <- hypval

  hypmat
}
