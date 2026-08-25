suppressPackageStartupMessages({
  library(here); library(clusterProfiler); library(org.Hs.eg.db)
  library(dplyr); library(ggplot2); library(patchwork); library(stringr)
})
out_dir <- here("results", "figures", "intermediate", "figure2")
up <- readLines(file.path(out_dir, "consensus_up_16.txt"))
down <- readLines(file.path(out_dir, "consensus_down_21.txt"))

run_one <- function(symbols, direction) {
  mapping <- bitr(symbols, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db) |> distinct(SYMBOL, ENTREZID)
  if (nrow(mapping) != length(symbols)) warning(direction, ": not every symbol mapped uniquely.")
  ids <- unique(mapping$ENTREZID)
  frozen <- here("data", "figure_inputs", "figure2", "frozen_enrichment", direction, "KEGG.csv")
  use_live_kegg <- identical(tolower(Sys.getenv("USE_LIVE_KEGG", "false")), "true")
  kegg <- if (!use_live_kegg) {
    read.csv(frozen, stringsAsFactors = FALSE, check.names = FALSE)
  } else tryCatch(
    enrichKEGG(ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1),
    error = function(e) {
      if (!file.exists(frozen)) stop(e)
      message("KEGG live query unavailable; using frozen provenance table: ", frozen)
      read.csv(frozen, stringsAsFactors = FALSE, check.names = FALSE)
    }
  )
  ans <- list(
    GO_BP = enrichGO(ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
    GO_MF = enrichGO(ids, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
    GO_CC = enrichGO(ids, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
    KEGG = kegg
  )
  d <- file.path(out_dir, direction); dir.create(d, showWarnings = FALSE)
  write.csv(mapping, file.path(d, "SYMBOL_to_ENTREZID.csv"), row.names = FALSE)
  for (nm in names(ans)) write.csv(as.data.frame(ans[[nm]]), file.path(d, paste0(nm, ".csv")), row.names = FALSE)
  ans
}
res <- list(Up = run_one(up, "Up"), Down = run_one(down, "Down"))

top_terms <- function(x, direction, category, n = 5) {
  z <- as.data.frame(x)
  if (!nrow(z)) return(NULL)
  z |> arrange(pvalue) |> slice_head(n = n) |>
    transmute(Direction = .env$direction, Category = .env$category,
              Term = str_wrap(Description, 45), pvalue, p.adjust,
              Score = -log10(pmax(pvalue, .Machine$double.xmin)))
}
plot_df <- bind_rows(lapply(names(res), function(direction)
  bind_rows(lapply(names(res[[direction]]), function(category)
    top_terms(res[[direction]][[category]], direction, category)))))
write.csv(plot_df, file.path(out_dir, "Figure2_top5_enrichment_terms.csv"), row.names = FALSE)
plot_df$Term <- factor(plot_df$Term, levels = rev(unique(plot_df$Term[order(plot_df$pvalue)])))
p <- ggplot(plot_df, aes(Term, Score, fill = Category)) +
  geom_col(width = .68) + coord_flip() + facet_grid(Category ~ Direction, scales = "free_y", space = "free_y") +
  labs(x = NULL, y = expression(-log[10](italic(P))), title = "Functional enrichment of 37 consensus genes") +
  scale_fill_manual(values = c(GO_BP = "#D9B300", GO_MF = "#7A3E9D", GO_CC = "#E68A2E", KEGG = "#5B8FD1")) +
  theme_bw(base_size = 10) + theme(legend.position = "none", axis.text.y = element_text(face = "bold"), strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "Figure2_enrichment.png"), p, width = 14, height = 12, dpi = 300)
ggsave(file.path(out_dir, "Figure2_enrichment.pdf"), p, width = 14, height = 12)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Figure2_enrichment_sessionInfo.txt"))
message("Saved reproducible enrichment tables and panel.")
