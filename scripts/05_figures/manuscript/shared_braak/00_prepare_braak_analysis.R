suppressPackageStartupMessages({
  library(here); library(readxl); library(dplyr); library(tidyr)
  library(purrr); library(broom); library(sva)
})

input_file <- here("data", "figure_inputs", "braak", "BRAAK_DATA_4studies.xlsx")
out_dir <- here("results", "figures", "intermediate", "braak")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))

df <- read_excel(input_file)
meta_cols <- c("SampleID", "SEX", "AGE", "BRAAK1", "PATHDX", "REGION", "BATCH")
stopifnot(all(meta_cols %in% names(df)), !anyDuplicated(df$SampleID))

# Accept numeric or Roman Braak labels. The archived workbook currently uses 0-6.
braak_map <- c("0"=0,"I"=1,"1"=1,"II"=2,"2"=2,"III"=3,"3"=3,
               "IV"=4,"4"=4,"V"=5,"5"=5,"VI"=6,"6"=6)
braak_raw <- toupper(trimws(as.character(df$BRAAK1)))
df$BRAAK_numeric <- unname(braak_map[braak_raw])
if (anyNA(df$BRAAK_numeric)) stop("Unrecognized Braak labels: ", paste(unique(braak_raw[is.na(df$BRAAK_numeric)]), collapse=", "))

df <- df |>
  mutate(
    SampleID = as.character(SampleID),
    SEX = recode(toupper(trimws(SEX)), "MALE"="M", "FEMALE"="F"),
    PATHDX = recode(toupper(trimws(PATHDX)), "ALZHEIMERS"="AD", "HEALTHY"="CON"),
    REGION = toupper(trimws(REGION)), BATCH = trimws(BATCH),
    Braak_bin = factor(ifelse(BRAAK_numeric <= 3, "Low(<=3)", "High(>3)"),
                       levels=c("Low(<=3)","High(>3)"))
  )

# All staged participants are required: AD-only filtering leaves no low-Braak
# observations in three regions. Diagnosis is preserved during ComBat but is
# not added to the Braak association model because it is confounded with stage.
gene_cols <- setdiff(names(df), c(meta_cols, "BRAAK_numeric", "Braak_bin"))
stopifnot(length(gene_cols)==37L, !anyNA(df[,c(meta_cols,"BRAAK_numeric")]),
          all(vapply(df[,gene_cols], is.numeric, logical(1))))
group_counts <- df |> count(BATCH, REGION, Braak_bin, name="n") |> complete(BATCH, REGION, Braak_bin, fill=list(n=0))
region_counts <- df |> count(REGION, Braak_bin, name="n") |> complete(REGION, Braak_bin, fill=list(n=0))
if (any(region_counts$n==0)) stop("At least one region lacks a low/high Braak group.")

expr <- as.matrix(df[,gene_cols]); rownames(expr) <- df$SampleID
mod <- model.matrix(~ PATHDX + AGE + SEX, data=df)
expr_z <- scale(expr)
expr_combat_z <- t(ComBat(dat=t(expr_z), batch=df$BATCH, mod=mod,
                          par.prior=TRUE, prior.plots=FALSE))
expr_combat <- t(ComBat(dat=t(expr), batch=df$BATCH, mod=mod,
                        par.prior=TRUE, prior.plots=FALSE))
rownames(expr_combat_z) <- rownames(expr_combat) <- df$SampleID

overall_ancova <- map_dfr(gene_cols, function(g) {
  dat <- df |> mutate(Expression=expr_combat_z[SampleID,g])
  tidy(aov(Expression ~ Braak_bin + REGION + AGE + SEX, data=dat)) |>
    filter(term=="Braak_bin") |> transmute(Gene=g, statistic, p_value=p.value)
}) |> mutate(FDR=p.adjust(p_value,"BH")) |> arrange(FDR)

regional_ancova <- split(df, df$REGION) |> imap_dfr(function(dat, region) {
  map_dfr(gene_cols, function(g) {
    model_dat <- dat |> mutate(Expression=expr_combat_z[SampleID,g])
    fit <- lm(Expression ~ Braak_bin + AGE + SEX, data=model_dat)
    co <- tidy(fit) |> filter(term=="Braak_binHigh(>3)")
    if (!nrow(co)) return(tibble(Gene=g, Region=region, statistic=NA_real_, p_value=NA_real_))
    transmute(co, Gene=g, Region=region, statistic, p_value=p.value)
  })
}) |> group_by(Region) |> mutate(FDR=p.adjust(p_value,"BH")) |> ungroup()

long_raw <- as.data.frame(expr_combat) |> tibble::rownames_to_column("SampleID") |>
  pivot_longer(all_of(gene_cols), names_to="Gene", values_to="Expression") |>
  left_join(df |> select(SampleID, Braak_bin, Region=REGION), by="SampleID")

overall_fc <- long_raw |> group_by(Gene,Braak_bin) |> summarise(mean=mean(Expression),.groups="drop") |>
  pivot_wider(names_from=Braak_bin,values_from=mean) |>
  mutate(log2FC=`High(>3)`-`Low(<=3)`, signed_FC=ifelse(log2FC>=0,2^log2FC,-2^abs(log2FC))) |>
  left_join(overall_ancova |> select(Gene,p_value,FDR),by="Gene")
regional_fc <- long_raw |> group_by(Gene,Region,Braak_bin) |> summarise(mean=mean(Expression),.groups="drop") |>
  pivot_wider(names_from=Braak_bin,values_from=mean) |>
  mutate(log2FC=`High(>3)`-`Low(<=3)`, signed_FC=ifelse(log2FC>=0,2^log2FC,-2^abs(log2FC))) |>
  left_join(regional_ancova |> select(Gene,Region,p_value,FDR),by=c("Gene","Region"))

saveRDS(list(df=df, expr_combat=expr_combat, expr_combat_z=expr_combat_z,
             overall_ancova=overall_ancova, regional_ancova=regional_ancova,
             overall_fc=overall_fc, regional_fc=regional_fc),
        file.path(out_dir,"braak_analysis_objects.rds"))
write.csv(group_counts,file.path(out_dir,"sample_counts_by_dataset_region_braak.csv"),row.names=FALSE)
write.csv(region_counts,file.path(out_dir,"sample_counts_by_region_braak.csv"),row.names=FALSE)
write.csv(overall_ancova,file.path(out_dir,"overall_Braak_ANCOVA.csv"),row.names=FALSE)
write.csv(regional_ancova,file.path(out_dir,"regional_Braak_ANCOVA.csv"),row.names=FALSE)
write.csv(overall_fc,file.path(out_dir,"overall_Braak_fold_changes.csv"),row.names=FALSE)
write.csv(regional_fc,file.path(out_dir,"regional_Braak_fold_changes.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"Braak_sessionInfo.txt"))
message("Prepared Braak analysis: ",nrow(df)," staged participants, 37 genes, 4 regions.")
