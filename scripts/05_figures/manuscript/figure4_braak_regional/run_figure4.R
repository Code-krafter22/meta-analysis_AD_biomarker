root<-file.path(here::here(),"scripts","05_figures","manuscript")
for(f in c(file.path("shared_braak","00_prepare_braak_analysis.R"),file.path("figure4_braak_regional","01_regional_waterfall.R"))){message("Running: ",f);sys.source(file.path(root,f),envir=globalenv())}
