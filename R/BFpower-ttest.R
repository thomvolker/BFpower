#' Calculate power for a one-sample Bayesian t-test with informative hypotheses
#'
#' @param n The sample size.
#' @param d The effect size.
#' @param sd The standard deviation.
#' @param threshold The threshold value for the Bayes factor (i.e., which Bayes
#' factor value is considered as convincing evidence for or against the hypothesis
#' of interest). By default, Bayes factors smaller than 1 are considered as evidence
#' against the hypothesis of interest.
#' @param fraction The fraction of information in the data used to construct the prior
#' distribution. The default value 1 denotes the minimal fraction, 2 denotes twice the
#' minimal fraction, and so on. See [bain::bain()] for more details.
#' @param power The target power of the test. If NULL, the function calculates the
#' power for the given sample size, effect size, and standard deviation.
#' @param hypothesis The hypothesis of interest. The default hypothesis denotes
#' that the effect size is equal to 0. See [bain::bain()] for more details.
#' @param BFtype The type of Bayes factor to calculate. Currently only the
#' approximate adjusted fractional Bayes factor is supported.
#'
#' @return The power of the test, or the sample size required to achieve said power.
#' @export
BFpower.one.sample.t <- function(n = NULL,
                                 d = 0.1,
                                 sd = 1,
                                 threshold = 3,
                                 fraction = 1,
                                 power = NULL,
                                 hypothesis = "d=0",
                                 BFtype = c("AAFBF", "AFBF")) {

  # Parse the hypothesis of interest using the bain functionality
  hypmat <- bain:::parse_hypothesis("d", hypothesis)
  # Add some additional checks and manipulations to the constraint-matrices
  hypmat <- check_hypothesis(hypmat)

  # TODO: allow for the adjusted fractional Bayes factor
  BFtype <- match.arg(BFtype, c("AAFBF", "AFBF"))
  # Check input arguments for anomalies
  check_inputs(n, d, sd, threshold, fraction, power, BFtype)
  # Calculate support against hypothesis of interest if the threshold value is smaller than 1.
  support <- ifelse(threshold > 1, "for", "against")

  # Initialize message for output
  message <- NULL

  if (!is.null(n)) {
    power <- one.sided.t.power(n, d, sd, threshold, fraction, BFtype, hypmat, support)
  } else {
    minpower <- one.sided.t.power(2, d, sd, threshold, fraction, BFtype, hypmat, support)[1]
    maxpower <- one.sided.t.power(1e6, d, sd, threshold, fraction, BFtype, hypmat, support)[1]
    if (isTRUE(all(c(minpower, maxpower) < 1e-16))) {
      message <- paste0("Power is and remains zero for all sample sizes. Perhaps you want support against the hypothesis of interest, or lower the threshold.")
    } else {
      n <- stats::uniroot(\(x) {
        one.sided.t.power(x, d, sd, threshold, fraction, BFtype, hypmat, support)[1] - power
      }, c(2, 1e6), extendInt = "upX")$root
    }
  }
  list(n = n, d = d, sd = sd, threshold = threshold, fraction = fraction, power = power[1], hypothesis = hypothesis)
}

#' Calculate the power of a one-sample t-test
#' @noRd
one.sided.t.power <- function(n, d, sd, threshold, fraction, BFtype, hypmat, support) {
  if (hypmat$hyptype == "eq") {
    std_mu <- (d-hypmat$hyp_mat[[1]][1,2])/(sd/sqrt(n))
    mu_critical2 <- 2 * log(sqrt(n/fraction)/threshold)
    mu_critical <- ifelse(mu_critical2 < 0, 0, sqrt(mu_critical2))
    if (support == "for") {
      p <- stats::pt(mu_critical, n-1, std_mu) - stats::pt(-mu_critical, n-1, std_mu)
    } else {
      p <- stats::pt(-mu_critical, n-1, std_mu) + stats::pt(mu_critical, n-1, std_mu, lower.tail = FALSE)
    }
    p_comp <- p
  } else if (hypmat$hyptype == "one.sided") {
    if (threshold > 2) {
      stop("For one-sided hypotheses, the Bayes factor cannot exceed 2, and thus the threshold should be smaller than 2.")
    }

    # bain only uses greater than inequalities, and thus transforms the effect size
    # x to -x if the hypothesis is "d<x". This corresponds to an equivalent
    # hypothesis stating -d > -x. Since the value in hypmat is transformed if
    # necessary, we also transform the effect size here if necessary.
    p <- stats::pt(stats::qnorm(threshold/2),
                   df = n-1,
                   ncp = (sign(hypmat$hyp_mat[[1]][1,1])*d - hypmat$hyp_mat[[1]][1,2])/(sd/sqrt(n)),
                   lower.tail = ifelse(support == "for", FALSE, TRUE))

    p_comp <- NULL
  }
  c(unconstrained = p, complement = p_comp)
}



