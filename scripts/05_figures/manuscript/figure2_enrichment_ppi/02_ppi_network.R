# HISTORICAL PROVENANCE COPY — DO NOT RUN. Use
# 02_reproducible_ppi_network.R through run_figure2.R.
# ============================================================
# Publication-Quality PPI Network Plot (k=2 Clustering)
# Matches STRING DB layout with enhanced aesthetics
# ============================================================

# Install if needed (uncomment):
# install.packages(c("igraph", "ggraph", "tidygraph", "ggrepel", "scales"))

library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(scales)

# ---- 1) Read STRING short TSV ----
edges <- read.delim("C:\\Users\\kashvichirag\\Box\\Dr. ZHANG\\alzheimer_meta\\ad_meta\\FIGURES\\ppi\\string_interactions_short.tsv",
                    header = TRUE, stringsAsFactors = FALSE)

colnames(edges) <- gsub("^X\\.|^#", "", colnames(edges))

edge_df <- data.frame(
  from           = edges$node1,
  to             = edges$node2,
  combined_score = as.numeric(edges$combined_score),
  stringsAsFactors = FALSE
)

# ---- 2) All 37 consensus genes ----
all_consensus <- c(
  "CHML","KANK2","PRELP","HSPB1","PIK3R5","AEBP1","GFAP",
  "ITFG1","TCEA3","CCDC102A","ELOVL4","KLF15","NAP1L5",
  "NEUROD6","HPRT1","NRN1","MAS1","RPH3A","HSPB7","NUPR1",
  "ADAM33","COL25A1","RAB3C","SERPINI1","STAT4","SCG2",
  "TMPRSS5","RGS4","PRX","RAB3B","NCALD","OPN3","CLDN9",
  "TRIM36","GAD2","MRGPRF","GAD1"
)

# ---- 3) Correct cluster assignments from STRING k=2 ----
cluster1_genes <- c("GAD1", "GAD2", "GFAP", "NEUROD6", "HSPB1",
                    "HPRT1", "PRELP", "RGS4")

cluster2_genes <- c("CHML", "RAB3B", "RAB3C", "RPH3A", "OPN3")

disconnected_genes <- setdiff(all_consensus, c(cluster1_genes, cluster2_genes))

# ---- 4) Build graph ----
g <- graph_from_data_frame(edge_df, directed = FALSE)

connected_genes <- V(g)$name
to_add <- setdiff(all_consensus, connected_genes)
g <- add_vertices(g, length(to_add), name = to_add)

stopifnot("combined_score" %in% edge_attr_names(g))

# ---- 5) Annotate graph ----
tg <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(
    degree = centrality_degree(),
    cluster_name = case_when(
      name %in% cluster1_genes      ~ "Cluster 1: GABA Synthesis (n=8)",
      name %in% cluster2_genes      ~ "Cluster 2: Synaptic Vesicle Membrane (n=5)",
      name %in% disconnected_genes  ~ "No Interaction (n=24)",
      TRUE                          ~ "No Interaction (n=24)"
    ),
    connected = degree > 0
  )

# ---- 6) Custom layout ----
set.seed(42)
layout_fr <- create_layout(tg, layout = "fr")

disc_idx <- which(layout_fr$degree == 0)
conn_idx <- which(layout_fr$degree > 0)

if (length(disc_idx) > 0 && length(conn_idx) > 0) {
  cx <- mean(layout_fr$x[conn_idx])
  cy <- mean(layout_fr$y[conn_idx])
  spread <- max(
    max(layout_fr$x[conn_idx]) - min(layout_fr$x[conn_idx]),
    max(layout_fr$y[conn_idx]) - min(layout_fr$y[conn_idx])
  )
  r <- spread * 0.65
  
  n_disc <- length(disc_idx)
  angles <- seq(0, 2 * pi, length.out = n_disc + 1)[1:n_disc]
  angles <- angles + pi / 6
  
  layout_fr$x[disc_idx] <- cx + r * cos(angles)
  layout_fr$y[disc_idx] <- cy + r * sin(angles)
}

# ---- 7) Color palette ----
cluster_colors <- c(
  "Cluster 1: GABA Synthesis (n=8)"             = "#FC1D7B",
  "Cluster 2: Synaptic Vesicle Membrane (n=5)"  = "#69c2a1",
  "No Interaction (n=24)"                        = "#B3B8D8"
)

# ---- 8) Plot ----
ppi_plot <- ggraph(layout_fr) +
  
  # --- Edges (show legend for width) ---
  geom_edge_link(
    aes(width = combined_score, alpha = combined_score),
    colour      = "grey40",
    lineend     = "round"
  ) +
  scale_edge_width_continuous(
    name   = "Interaction Confidence\n(STRING Combined Score)",
    range  = c(0.5, 2.8),
    breaks = c(0.4, 0.6, 0.8, 1.0),
    labels = c("0.4 (Low)", "0.6 (Medium)", "0.8 (High)", "1.0 (Highest)"),
    limits = c(0.4, 1.0)
  ) +
  scale_edge_alpha_continuous(range = c(0.2, 0.75), guide = "none") +
  
  # --- Nodes: disconnected ---
  geom_node_point(
    data  = function(x) filter(x, !connected),
    aes(fill = cluster_name),
    shape  = 21,
    size   = 5.5,
    color  = "grey25",
    stroke = 0.9,
    alpha  = 0.6
  ) +
  
  # --- Nodes: connected (show legend for size) ---
  geom_node_point(
    data  = function(x) filter(x, connected),
    aes(fill = cluster_name, size = degree),
    shape  = 21,
    color  = "grey25",
    stroke = 0.7,
    alpha  = 0.8
  ) +
  scale_size_continuous(
    name   = "Node Degree\n(Number of Interactions)",
    range  = c(4,12),
    breaks = c(1,2,3,4,5,6)
  ) +
  
  # --- Labels: disconnected genes (lines pointing to circles) ---
  geom_node_text(
    data  = function(x) filter(x, !connected),
    aes(label = name),
    size          = 3,
    fontface      = "italic",
    color         = "grey1",
    repel         = TRUE,
    max.overlaps  = 100,
    point.padding = unit(0.8, "lines"),
    box.padding   = unit(0.5, "lines"),
    min.segment.length = 0,
    segment.color = "grey20",
    segment.size  = 0.4,
    force         = 3,
    force_pull    = 0.1,
    bg.color      = "white",
    bg.r          = 0.1
  ) +
  
  # --- Labels: connected genes (lines pointing to circles) ---
  geom_node_text(
    data  = function(x) filter(x, connected),
    aes(label = name),
    size          = 2.7,
    fontface      = "bold.italic",
    repel         = TRUE,
    max.overlaps  = 100,
    point.padding = unit(1.0, "lines"),
    box.padding   = unit(0.6, "lines"),
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size  = 0.5,
    force         = 2,
    force_pull    = 0.2,
    bg.color      = "white",
    bg.r          = 0.12
  ) +
  
  # --- Fill scale ---
  scale_fill_manual(
    values = cluster_colors,
    name   = "Functional Cluster\n(k=2)"
  ) +
  
  # --- Legend order and styling ---
  guides(
    fill = guide_legend(
      override.aes = list(size = 6, alpha = 0.4, stroke = 0.4),
      order = 1,
      direction = "vertical"
    ),
    size = guide_legend(
    size = guide_legend(
    override.aes = list(fill = "grey40", shape = 21, stroke = 0.5),
    order = 2
    ),
    order = 2,
    keyheight = unit(0.6, "cm"),
    keywidth  = unit(0.6, "cm")
    ),
    edge_width = guide_legend(order = 3)
  ) +
  
  # # --- Titles ---
  # labs(
  #   title    = "Protein\u2013Protein Interaction Network of Consensus Meta-DEGs",
  #   subtitle = "37 consensus genes | STRING v12.0 | k-means clustering (k = 2)",
  #   caption  = "Disconnected nodes (no STRING interactions at current confidence threshold) arranged in periphery"
  # ) +
  
  # --- Theme ---
  theme_void(base_size = 10) +
  theme(
    plot.title       = element_text(size = 20, face = "bold", hjust = 0.5,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 12, hjust = 0.5, color = "grey40",
                                    margin = margin(b = 14)),
    plot.caption     = element_text(size = 9, color = "grey50", hjust = 1,
                                    lineheight = 1.3, margin = margin(t = 14)),
    legend.position  = "right",
    legend.box       = "Vertical",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 10),
    # Add/modify these in your theme() block:
    legend.key.size  = unit(0.15, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.box.spacing = unit(0.3, "cm"),
    legend.spacing.x = unit(0.8, "cm"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(18, 20, 12, 20)
  )

# ---- 9) Save ----
ggsave("PPI_network_k2_publication.jpeg",  ppi_plot,
       width = 8.5, height = 6, dpi = 600, bg = "white")
ggsave("PPI_network_k2_publication.pdf",  ppi_plot,
       width = 16, height = 14, bg = "white")
ggsave("PPI_network_k2_publication.tiff", ppi_plot,
       width = 16, height = 14, dpi = 600, bg = "white",
       compression = "lzw")

print(ppi_plot)

message("\n=== Cluster summary ===")
message("Cluster 1 (GABA Synthesis): ", paste(cluster1_genes, collapse = ", "))
message("Cluster 2 (Synaptic Vesicle): ", paste(cluster2_genes, collapse = ", "))
message("No Interaction: ", paste(disconnected_genes, collapse = ", "))
message("\nFiles saved: PNG, PDF, TIFF (600 DPI)")
