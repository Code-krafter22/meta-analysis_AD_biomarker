suppressPackageStartupMessages({
  library(here); library(igraph); library(ggraph); library(tidygraph); library(dplyr)
})
input_file <- here("data", "figure_inputs", "figure2", "string_interactions_short.tsv")
out_dir <- here("results", "figures", "intermediate", "figure2")
signature <- read.csv(file.path(out_dir, "consensus_37_direction.csv"), stringsAsFactors = FALSE)
edges <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
names(edges) <- sub("^#", "", names(edges))
stopifnot(all(c("node1", "node2", "combined_score") %in% names(edges)))
edge_df <- transmute(edges, from = node1, to = node2, combined_score = as.numeric(combined_score)) |>
  filter(from %in% signature$SYMBOL, to %in% signature$SYMBOL)
g <- graph_from_data_frame(edge_df, directed = FALSE,
                           vertices = data.frame(name = signature$SYMBOL, Direction = signature$Direction))
set.seed(42)
p <- ggraph(as_tbl_graph(g), layout = "fr") +
  geom_edge_link(aes(width = combined_score, alpha = combined_score), colour = "grey45", show.legend = FALSE) +
  geom_node_point(aes(fill = Direction, size = centrality_degree()), shape = 21, colour = "grey20") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "italic") +
  scale_fill_manual(values = c(Up = "#D84A4A", Down = "#3B78B5")) +
  scale_size_continuous(range = c(4, 10), name = "Node degree") +
  theme_void() + labs(title = "STRING network of 37 consensus genes", fill = "Direction") +
  theme(plot.title = element_text(face = "bold", hjust = .5))
ggsave(file.path(out_dir, "Figure2_STRING_PPI.png"), p, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "Figure2_STRING_PPI.pdf"), p, width = 10, height = 8)
write.csv(edge_df, file.path(out_dir, "Figure2_STRING_edges_used.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Figure2_PPI_sessionInfo.txt"))
message("Saved reproducible STRING PPI panel using ", nrow(edge_df), " edges.")
