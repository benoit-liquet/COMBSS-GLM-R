library(ggplot2)
# ================================================================
# Extract results from the COMBSS output object
# ================================================================
selected_features <- result$selected_models
test_acc          <- result$test_accuracy
Kmax              <- length(selected_features)

# ================================================================
# Order genes by first appearance in the inclusion path 
# "Gene 1955 is the first to enter the model (at k=1) and
# persists throughout, followed by Gene 509 at k=2 and
# Gene 187 at k=3" 
# ================================================================
all_genes <- c()
for (k in 1:Kmax) {
  feats <- selected_features[[k]]  # numeric index, not character
  new_genes <- setdiff(feats, all_genes)
  all_genes <- c(all_genes, sort(new_genes))
}

gene_order <- paste0("Gene ", all_genes)
n_genes <- length(all_genes)

# ================================================================
# Map each gene to a fixed x position
# ================================================================
gene_to_pos <- setNames(seq_along(all_genes), as.character(all_genes))

# ================================================================
# Build tile data
# ================================================================
df_tiles <- data.frame(k = integer(), x_pos = numeric(),
                       stringsAsFactors = FALSE)
for (k in 1:Kmax) {
  feats <- selected_features[[k]]
  for (g in feats) {
    pos <- gene_to_pos[as.character(g)]
    if (!is.na(pos)) {
      df_tiles <- rbind(df_tiles, data.frame(k = k, x_pos = pos))
    }
  }
}

# ================================================================
# Accuracy labels placed to the right of the rightmost tile
# Perfect accuracy (100%) highlighted in red 
# ================================================================
df_acc <- data.frame(k = integer(), x_pos = numeric(),
                     acc_label = character(), colour = character(),
                     stringsAsFactors = FALSE)

for (k in 1:Kmax) {
  feats <- selected_features[[k]]
  positions <- gene_to_pos[as.character(feats)]
  positions <- positions[!is.na(positions)]
  if (length(positions) == 0) next
  rightmost <- max(positions)
  acc_val <- test_acc[k]
  df_acc <- rbind(df_acc, data.frame(
    k = k,
    x_pos = rightmost + 0.7,
    acc_label = sprintf("%.0f%%", acc_val * 100),
    colour = ifelse(acc_val >= 1.0, "#d62728", "grey30")
  ))
}

# ================================================================
# Find first k with perfect accuracy for subtitle 
# ================================================================
perfect_k <- which(test_acc >= 1.0)
if (length(perfect_k) > 0) {
  subtitle_text <- sprintf(
    "COMBSS multinomial: perfect classification at k = %d with %d / 2308 genes",
    perfect_k[1], perfect_k[1])
} else {
  best_k <- which.max(test_acc)
  subtitle_text <- sprintf(
    "COMBSS multinomial: best accuracy %.0f%% at k = %d",
    test_acc[best_k] * 100, best_k)
}

# ================================================================
# Plot: Best-subset inclusion path 
# ================================================================
p <- ggplot() +
  geom_tile(data = df_tiles,
            aes(x = x_pos, y = factor(k)),
            fill = "#1f77b4", colour = "white", linewidth = 0.3,
            width = 0.9, height = 0.85) +
  geom_text(data = df_acc,
            aes(x = x_pos, y = factor(k), label = acc_label,
                colour = colour),
            size = 3, hjust = 0, fontface = "bold",
            show.legend = FALSE) +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(as.character(1:Kmax))) +
  scale_x_continuous(
    breaks = seq_len(n_genes),
    labels = gene_order,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Best-subset inclusion path (Khan SRBCT dataset)",
    subtitle = subtitle_text,
    x = "Gene",
    y = "Model size (k)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0, colour = "grey40"),
    axis.title = element_text(size = 11),
    plot.margin = margin(10, 30, 10, 10)
  )

print(p)

ggsave("khan_inclusion_path.pdf", p, width = 12, height = 7)
ggsave("khan_inclusion_path.png", p, width = 12, height = 7, dpi = 300)
