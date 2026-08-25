# MetaVolcanoR 1.16.0 compatibility wrapper for ggplot2 >= 4.0.
# It uses MetaVolcanoR's plot function but avoids two upstream issues:
# metap::sumlog() silently returning NA without multtest, and the legacy S4
# slot requiring plot class "gg" rather than ggplot2's current class.
combining_mv_compatible <- function(diffexp, pcriteria = "p",
                                    foldchangecol = "effectsize",
                                    genenamecol = "SYMBOL", metafc = "Mean",
                                    metathr = 0.05, jobname = "MetaVolcano",
                                    outputfolder = ".", draw = "HTML") {
  renamed <- Map(function(x, nm) {
    x <- x[, c(genenamecol, pcriteria, foldchangecol), drop = FALSE]
    names(x)[names(x) == pcriteria] <- paste0(pcriteria, ".", nm)
    names(x)[names(x) == foldchangecol] <- paste0(foldchangecol, ".", nm)
    x
  }, diffexp, names(diffexp))

  merged <- Reduce(function(x, y) merge(x, y, by = genenamecol, all = TRUE), renamed)
  pcols <- grep(paste0("^", pcriteria, "\\."), names(merged), value = TRUE)
  fcols <- grep(paste0("^", foldchangecol, "\\."), names(merged), value = TRUE)

  # Fisher's method, identical to metap::sumlog(), without its optional dependency.
  merged$metap <- apply(merged[pcols], 1, function(z) {
    z <- as.numeric(z)
    z <- z[is.finite(z) & z > 0 & z <= 1]
    if (!length(z)) NA_real_ else
      stats::pchisq(-2 * sum(log(z)), df = 2 * length(z), lower.tail = FALSE)
  })
  merged$metafc <- apply(merged[fcols], 1, function(z) {
    z <- as.numeric(z)
    if (metafc == "Mean") mean(z, na.rm = TRUE) else median(z, na.rm = TRUE)
  })
  merged <- merged[is.finite(merged$metap) & is.finite(merged$metafc), , drop = FALSE]
  merged$idx <- merged$metafc * -log10(merged$metap)
  lo <- stats::quantile(merged$idx, metathr / 2, na.rm = TRUE)
  hi <- stats::quantile(merged$idx, 1 - metathr / 2, na.rm = TRUE)
  merged$degcomb <- ifelse(merged$idx < lo, "0.Down-regulated",
                           ifelse(merged$idx > hi, "2.Up-regulated", "1.Unperturbed"))
  merged <- merged[order(-abs(merged$idx)), ]

  gg <- MetaVolcanoR:::plot_mv(merged, NULL, genenamecol, TRUE, metafc)
  dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(outputfolder,
                       paste0("combining_method_MetaVolcano_", jobname,
                              if (draw == "HTML") ".html" else ".pdf"))
  if (draw == "HTML") {
    htmlwidgets::saveWidget(plotly::as_widget(plotly::ggplotly(gg)), outfile,
                            selfcontained = FALSE)
  } else {
    ggplot2::ggsave(outfile, gg, width = 4, height = 5)
  }

  list(metaresult = merged[, c(genenamecol, "metap", "metafc", "idx")],
       MetaVolcano = gg)
}
