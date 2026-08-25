library(here)
out_dir <- here("results", "figures", "intermediate", "figure2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# The 37-gene signature is the three-method intersection prepared for Figure 1,
# not every gene called by any one meta-analysis output. Direction is assigned
# from the random-effects meta_logFC.
effect_file <- here("results", "figures", "intermediate", "figure1", "meta_effects_37_genes.csv")
if (!file.exists(effect_file)) stop("Run Figure 1 input preparation first: ", effect_file)
effect <- read.csv(effect_file, stringsAsFactors = FALSE)
stopifnot(all(c("SYMBOL", "meta_logFC") %in% names(effect)), nrow(effect) == 37L)
up <- effect$SYMBOL[effect$meta_logFC > 0]
down <- effect$SYMBOL[effect$meta_logFC < 0]
stopifnot(length(up) == 16L, length(down) == 21L, length(unique(c(up, down))) == 37L)
writeLines(up, file.path(out_dir, "consensus_up_16.txt"))
writeLines(down, file.path(out_dir, "consensus_down_21.txt"))
write.csv(data.frame(SYMBOL = c(up, down), Direction = c(rep("Up", 16), rep("Down", 21))),
          file.path(out_dir, "consensus_37_direction.csv"), row.names = FALSE)
message("Prepared 16 upregulated and 21 downregulated consensus genes.")
