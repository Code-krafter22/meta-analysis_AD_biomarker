suppressPackageStartupMessages({library(here);library(dplyr);library(ggplot2)})
x <- readRDS(here("results","figures","intermediate","braak","braak_analysis_objects.rds"))
out <- here("results","figures","intermediate","figure3");dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- x$overall_fc |> mutate(Direction=ifelse(log2FC>=0,"Up","Down"),
  Star=case_when(FDR<=.001~"***",FDR<=.01~"**",FDR<=.05~"*",TRUE~""),
  label_y=signed_FC+ifelse(signed_FC>=0,.12,-.12))
d$Gene <- factor(d$Gene, levels=d$Gene[order(d$log2FC)])
p <- ggplot(d,aes(Gene,signed_FC,fill=Direction))+geom_col(width=.82)+geom_hline(yintercept=0)+
  geom_text(aes(y=label_y,label=Star),fontface="bold",size=4)+coord_flip()+
  scale_fill_manual(values=c(Up="tomato4",Down="yellowgreen"))+
  labs(x=NULL,y="Signed fold change: Braak High vs Low",fill=NULL)+theme_minimal(base_size=11)+
  theme(axis.text.y=element_text(face="bold"),legend.position="top")
ggsave(file.path(out,"Figure3A_overall_Braak_waterfall.png"),p,width=8,height=7,dpi=300)
ggsave(file.path(out,"Figure3A_overall_Braak_waterfall.pdf"),p,width=8,height=7)
message("Saved Figure 3A; stars indicate BH-adjusted ANCOVA significance.")
