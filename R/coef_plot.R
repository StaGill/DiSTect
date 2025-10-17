#' Plot top beta coefficients
#'
#' @param fit   A stanfit object
#' @param n number of coefficients to be plotted
#' @param matrix_x a matrix of covariates, with the coordinates placed in the final two columns of the matrix
#' @return a ggplot object
#' @import ggplot2 dplyr
#' @export
#' @examples \donttest{}

plot_coef <- function(fit, matrix_x, n = 30) {
  
  fit_summary <- rstan::summary(fit)$summary
  beta_rows   <- grep("^beta\\[", rownames(fit_summary))
  beta_data   <- fit_summary[beta_rows, ]
  
  plot_data <- data.frame(
    gene = colnames(matrix_x[, 1:(ncol(matrix_x) - 2)]),
    coef = beta_data[, "mean"],
    sd   = beta_data[, "sd"]
  )
  
  # error bar = |mean / sd|
  plot_data$error     <- abs(plot_data$coef / plot_data$sd)
  plot_data$abs_effect <- abs(plot_data$coef)
  
  ## keep top n by absolute effect size
  plot_data <- plot_data |>
    arrange(desc(abs_effect)) |>
    slice_head(n = n)
  
  ## order factors for plotting
  plot_data$gene <- factor(plot_data$gene,
                           levels = plot_data$gene[order(-plot_data$abs_effect)])
  
  ggplot(plot_data, aes(x = gene, y = coef)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point() +
    geom_errorbar(aes(ymin = coef - error,
                      ymax = coef + error),
                  width = 0.2, size = 1.5) +
    theme_minimal() +
    theme_bw() +
    labs(x = "Gene Names", y = "Coefficient Values") +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_blank()
    )
}


