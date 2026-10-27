draw_network <- function(adjacency, title, subtitle = "", color = "#176B87") {
  p <- nrow(adjacency)
  angle <- pi / 2 - 2 * pi * (0:(p - 1L)) / p
  xy <- cbind(cos(angle), sin(angle))
  plot.new(); plot.window(c(-1.28, 1.28), c(-1.28, 1.28), asp = 1)
  edges <- which(adjacency & upper.tri(adjacency), arr.ind = TRUE)
  if (nrow(edges)) segments(xy[edges[, 1], 1], xy[edges[, 1], 2],
                           xy[edges[, 2], 1], xy[edges[, 2], 2],
                           col = adjustcolor(color, alpha.f = 0.7), lwd = 2)
  points(xy, pch = 21, bg = "white", col = color, cex = 3.6, lwd = 1.7)
  text(xy, labels = seq_len(p), cex = 0.85, col = "#172B3A")
  title(main = title, cex.main = 1.05, col.main = "#172B3A")
  mtext(subtitle, side = 1, line = 0.2, cex = 0.8, col = "#526472")
}

save_network_plot <- function(data, fit, filename, edge_fraction = 0.1) {
  rows <- 1L + ceiling(fit$K / 3L)
  png(filename, width = 1500, height = 450 * rows + 100, res = 150)
  on.exit(dev.off())
  par(mfrow = c(rows, 3L), mar = c(2, 0.5, 2.5, 0.5), oma = c(1, 0, 3, 0), bg = "white")
  for (k in 1:3) draw_network(data$graphs$adjacency[[k]], paste("True network", k),
                               paste(sum(data$graphs$adjacency[[k]]) / 2, "edges"))
  tab <- table(factor(fit$labels, levels = seq_len(fit$K)), factor(data$network, levels = 1:3))
  for (j in seq_len(fit$K)) {
    adjacency <- precision_edges(fit$precision[[j]], edge_fraction)
    draw_network(adjacency, paste("Estimated group", j),
      paste0("n = ", sum(tab[j, ]), " | true-group counts: ", paste(tab[j, ], collapse = "/")),
      color = "#987045")
  }
  unused <- 3L * ceiling(fit$K / 3L) - fit$K
  if (unused > 0) for (j in seq_len(unused)) plot.new()
  mtext(paste0(toupper(data$graphs$type), "  |  HMM-glasso  |  true K = 3; fitted K = ", fit$K,
               "; L = ", paste(unique(fit$L), collapse = ",")),
        side = 3, outer = TRUE, line = 1, font = 2, cex = 1.2)
  mtext("Estimated edges use |Omega[i,j]| > 0.1 x min(diag(Omega)); labels are arbitrary.",
        side = 1, outer = TRUE, line = -0.2, cex = 0.75)
}

save_diagnostics <- function(data, selected, filename) {
  png(filename, width = 1500, height = 600, res = 150)
  on.exit(dev.off())
  par(mfrow = c(1, 2), mar = c(4.5, 4, 3, 1), bg = "white")
  if (!is.null(selected$dendrogram)) {
    plot(selected$dendrogram, main = "Covariance dendrogram", xlab = "Initial MM-glasso components", sub = "", cex = 0.8)
    abline(h = selected$cut_height, col = "#A75131", lty = 2, lwd = 2)
  } else {
    plot.new(); text(0.5, 0.5, paste("Fixed K =", selected$fit$K))
  }
  tab <- table(factor(selected$fit$labels, levels = seq_len(selected$fit$K)), factor(data$network, levels = 1:3))
  barplot(t(tab), col = c("#176B87", "#72A98F", "#D4A15A"),
    xlab = "Estimated HMM group", ylab = "Observations", main = "True network composition",
    legend.text = paste("True", 1:3), args.legend = list(cex = 0.7, bty = "n"))
}
