qq_plot <- function(p, title = "QQ plot", col_abline = "darkred", CB = TRUE, 
                    col_CB = "gray80", CB_level = 0.95, thinning = TRUE, ...) {
  
  if (any(is.na(p))) p <- p[!is.na(p)]
  if (any(p < 0 | p > 1)) stop("p-values should be in [0,1]")
  if (any(p == 0)) {
    warning("There are zero p-values that won't be displayed")
    p <- p[p > 0]
  }
  
  n <- length(p)
  observed <- sort(-log10(p))
  expected <- -log10((n:1) / (n + 1))
  
  # Confidence bands
  if (CB) {
    k <- n * 10^seq(-log10(n), 0, length = 200)
    hi <- -log10(qbeta(0.5 - CB_level / 2, k, n + 1 - k))
    lo <- -log10(qbeta(0.5 + CB_level / 2, k, n + 1 - k))
    conf_band <- data.frame(
      expected = c(-log10(k / (n + 1)), rev(-log10(k / (n + 1)))),
      observed = c(lo, rev(hi))
    )
  }
  
  # Prepare data for thinning
  if (thinning) {
    observed <- observed[seq(1, length(observed), length.out = min(10000, n))]
    expected <- expected[seq(1, length(expected), length.out = min(10000, n))]
  }
  
  qq_data <- data.frame(expected = expected, observed = observed)
  
  # Plot
  source("scripts/plotting_theme.R")
  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(size = 2) +
    geom_abline(slope = 2, intercept = 0, color = col_abline) +
    labs(
      title = title,
      x = expression(paste("Expected ", -log[10](p))),
      y = expression(paste("Observed ", -log[10](p)))
    ) 
  
  if (CB) {
    p <- p +
      geom_polygon(data = conf_band, aes(x = expected, y = observed), fill = col_CB, alpha = 0.5)
  }
  
  return(p)
}


