#' Network plot
#'
#' @param model A fitted `stan` model object returned from \code{dsgd()}.
#' @param interaction A vector specifying covariate names for interaction terms.
#' @param label_size Size of the node labels in the plot.
#' @param node_color Color of the nodes in the network plot.
#'
#' @return A ggplot object showing the interaction network.
#' @export
#'
#' @importFrom rstan summary
#' @importFrom GGally ggnet2
#' @importFrom network network
#' @examples \donttest{}

plot_network <- function(model, interaction, label_size = 13, node_color = '#C1D8E7') {
  
  coef_summary <- rstan::summary(model)$summary
  beta_rows <- grep("^beta\\[\\d+\\]$", rownames(coef_summary), value = TRUE)
  
  # Identify interaction term rows
  n_interactions <- choose(length(interaction), 2)
  interaction_rows <- tail(beta_rows, n_interactions)
  
  # Rename rows for interaction terms
  interaction_pairs <- combn(interaction, 2, simplify = FALSE)
  interaction_terms <- sapply(interaction_pairs, function(pair) paste(pair, collapse = "*"))
  rownames(coef_summary)[match(interaction_rows, rownames(coef_summary))] <- interaction_terms
  
  summary_interaction <- coef_summary[interaction_terms, ]
  
  # Create edges data frame
  edges <- data.frame(
    from = sapply(strsplit(rownames(summary_interaction), "\\*"), `[`, 1),
    to = sapply(strsplit(rownames(summary_interaction), "\\*"), `[`, 2),
    weight = abs(summary_interaction[, "mean"])
  )
  
  net <- network(edges, directed = FALSE)
  ggnet2(net, label = TRUE, edge.size = "weight", color = node_color, size = label_size)
}

