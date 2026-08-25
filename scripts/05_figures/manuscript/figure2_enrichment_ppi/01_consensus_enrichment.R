# HISTORICAL PROVENANCE COPY — DO NOT RUN. Use run_figure2.R.
# ============================================================
# ClusterProfiler enrichment (UP/DOWN) + ONE figure:
# GO_BP, GO_MF, GO_CC, KEGG (Top 5 UP vs Top 5 DOWN side-by-side)
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

# ----------------------------
# 1) Input gene lists
# ----------------------------
up_genes <- c("ADAM33","AEBP1","CCDC102A","CLDN9","GFAP","HSPB1","HSPB7","KANK2","KLF15",
              "MRGPRF","NUPR1","PIK3R5","PRELP","PRX","TCEA3","TMPRSS5")

down_genes <- c("CHML","COL25A1","ELOVL4","GAD1","GAD2","HPRT1","ITFG1","MAS1","NAP1L5","NCALD",
                "NEUROD6","NRN1","OPN3","RAB3B","RAB3C","RGS4","RPH3A","SCG2","SERPINI1","STAT4","TRIM36")

# ----------------------------
# 2) SYMBOL -> ENTREZ mapping
# ----------------------------
to_entrez <- function(symbols) {
  bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>%
    distinct(SYMBOL, ENTREZID)
}

up_map   <- to_entrez(up_genes)
down_map <- to_entrez(down_genes)

up_entrez   <- unique(up_map$ENTREZID)
down_entrez <- unique(down_map$ENTREZID)

message("UP mapped:   ", nrow(up_map),   "/", length(up_genes))
message("DOWN mapped: ", nrow(down_map), "/", length(down_genes))

# ----------------------------
# 3) Enrichment runner (NO filtering at this stage)
#    Important: pvalueCutoff=1, qvalueCutoff=1 so results don’t go empty
# ----------------------------
run_enrich <- function(entrez_ids, label, outdir="pathway_clusterprofiler") {
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  od <- file.path(outdir, label)
  dir.create(od, showWarnings=FALSE, recursive=TRUE)
  
  ego_bp <- enrichGO(entrez_ids, OrgDb=org.Hs.eg.db, ont="BP",
                     pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
                     readable=TRUE)
  
  ego_mf <- enrichGO(entrez_ids, OrgDb=org.Hs.eg.db, ont="MF",
                     pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
                     readable=TRUE)
  
  ego_cc <- enrichGO(entrez_ids, OrgDb=org.Hs.eg.db, ont="CC",
                     pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
                     readable=TRUE)
  
  ekegg <- enrichKEGG(entrez_ids, organism="hsa",
                      pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1)
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(ekegg, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
  }
  
  # Save tables
  write.csv(as.data.frame(ego_bp), file.path(od, "GO_BP.csv"), row.names=FALSE)
  write.csv(as.data.frame(ego_mf), file.path(od, "GO_MF.csv"), row.names=FALSE)
  write.csv(as.data.frame(ego_cc), file.path(od, "GO_CC.csv"), row.names=FALSE)
  write.csv(as.data.frame(ekegg),  file.path(od, "KEGG.csv"),  row.names=FALSE)
  
  # Save dotplots (these will still show p.adjust colorbar)
  save_plot <- function(p, fn) ggsave(file.path(od, fn), p, width=9, height=6, dpi=300)
  
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0)
    save_plot(dotplot(ego_bp, showCategory=20) + ggtitle(paste0(label," | GO BP")), "dotplot_GO_BP.png")
  if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0)
    save_plot(dotplot(ego_mf, showCategory=20) + ggtitle(paste0(label," | GO MF")), "dotplot_GO_MF.png")
  if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0)
    save_plot(dotplot(ego_cc, showCategory=20) + ggtitle(paste0(label," | GO CC")), "dotplot_GO_CC.png")
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0)
    save_plot(dotplot(ekegg, showCategory=20) + ggtitle(paste0(label," | KEGG")), "dotplot_KEGG.png")
  
  invisible(list(GO_BP=ego_bp, GO_MF=ego_mf, GO_CC=ego_cc, KEGG=ekegg))
}

res_up   <- run_enrich(up_entrez,   "UP")
res_down <- run_enrich(down_entrez, "DOWN")

# ----------------------------
# 4) ONE combined barplot figure (Top 5) for GO_BP/MF/CC/KEGG
#    We choose Top 5 AFTER enrichment (so UP won’t be empty)
# ----------------------------
topN <- 5
panels <- c("GO_BP","GO_MF","GO_CC","KEGG")

# ============================================================
# TWO SEPARATE FIGURES:
# 1) UPREGULATED: GO_BP, GO_MF, GO_CC, KEGG (Top 5 each)
# 2) DOWNREGULATED: GO_BP, GO_MF, GO_CC, KEGG (Top 5 each)
# + Different colors for GO vs KEGG
# Requires: res_up and res_down (from run_enrich())
# Outputs:
#   Top5_UP_GO_BP_GO_MF_GO_CC_KEGG_barplot.png
#   Top5_DOWN_GO_BP_GO_MF_GO_CC_KEGG_barplot.png
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

topN <- 5
panels <- c("GO_BP","GO_MF","GO_CC","KEGG")

# Color map: GO panels one color, KEGG another color
panel_fill <- c(
  GO_BP = "yellow3",
  GO_MF = "purple4",
  GO_CC = "orange2",
  KEGG  = "cornflowerblue"
)

# ---- Styling tweaks requested:
# 1) Bold main title + panel titles (GO_BP/MF/CC/KEGG)
# 2) Increase gap between bars
# 3) Make left-side pathway terms bold
# 4) Reduce plot width but keep readable

library(ggplot2)
library(stringr)
library(patchwork)
library(dplyr)

get_top_df <- function(enrich_obj, topN = 5) {
  if (is.null(enrich_obj)) return(NULL)
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>%
    mutate(Score = -log10(pvalue),
           Term  = str_wrap(Description, width = 55)) %>%
    arrange(pvalue)
  
  df[seq_len(min(topN, nrow(df))), , drop = FALSE]
}

make_bar <- function(df, title_txt, fill_col) {
  if (is.null(df) || nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(title_txt) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05)))+
        theme(
          plot.title = element_text(face = "bold", size = 11, hjust = 0.5)
        )
    )
  }
  
  df$Term <- factor(df$Term, levels = rev(df$Term))
  
  ggplot(df, aes(x = Term, y = Score)) +
    # gap between bars: reduce bar width (smaller width => more gap)
    geom_col(fill = fill_col, width = 0.6) +
    coord_flip() +
    labs(title = title_txt, x = NULL, y = "-log10(p-value)") +
    theme_bw(base_size = 10) +
    
    theme(
      # 1) bold panel title
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      
      # 3) left-side pathway terms bold
      axis.text.y = element_text(face = "bold", size = 12, lineheight = 0.9),
      
      # make y-axis numbers readable
      axis.text.x = element_text(size = 9),
      
      # a bit of spacing
      plot.margin = margin(6, 6, 6, 6)
    )
}

# ----------------------------
# UPREGULATED FIGURE (4 rows, 1 column)
# ----------------------------
up_plots <- lapply(panels, function(pn) {
  df_up <- get_top_df(res_up[[pn]], topN = topN)
  make_bar(df_up, pn, panel_fill[[pn]])   # pn title is bold via theme()
})

fig_up <- wrap_plots(up_plots, ncol = 1) +
  plot_annotation(title = "UPREGULATED GENES") &
  theme(
    # 1) bold main title
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# 4) reduced width, still readable
ggsave(
  filename = paste0("Top", topN, "_UP_GO_BP_GO_MF_GO_CC_KEGG_barplot.png"),
  plot = fig_up,
  width = 7.5, height = 8, dpi = 600
)

# ----------------------------
# DOWNREGULATED FIGURE (4 rows, 1 column)
# ----------------------------
down_plots <- lapply(panels, function(pn) {
  df_down <- get_top_df(res_down[[pn]], topN = topN)
  make_bar(df_down, pn, panel_fill[[pn]])
})

fig_down <- wrap_plots(down_plots, ncol = 1) +
  plot_annotation(title = "DOWNREGULATED GENES") &
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(
  filename = paste0("Top", topN, "_DOWN_GO_BP_GO_MF_GO_CC_KEGG_barplot.png"),
  plot = fig_down,
  width = 7.5, height = 8, dpi = 600
)

fig_up
fig_down


# Optional sanity check:
nrow(as.data.frame(res_up$GO_BP)); head(as.data.frame(res_up$GO_BP)[,c("Description","pvalue","p.adjust")])
