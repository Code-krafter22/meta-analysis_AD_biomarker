suppressPackageStartupMessages({library(here);library(dplyr);library(ggplot2)})
x<-readRDS(here("results","figures","intermediate","braak","braak_analysis_objects.rds"))
out<-here("results","figures","intermediate","figure4");dir.create(out,recursive=TRUE,showWarnings=FALSE)
main<-here("results","figures","main");dir.create(main,recursive=TRUE,showWarnings=FALSE)
d<-x$regional_fc |> mutate(Direction=ifelse(log2FC>=0,"Up","Down"),Star=case_when(FDR<=.001~"***",FDR<=.01~"**",FDR<=.05~"*",TRUE~""),label_y=signed_FC+ifelse(signed_FC>=0,.12,-.12))
order<-d |> group_by(Gene) |> summarise(mean_log2FC=mean(log2FC),.groups="drop") |> arrange(mean_log2FC) |> pull(Gene)
d<-d |> mutate(Gene=factor(Gene,levels=order))
p<-ggplot(d,aes(Gene,signed_FC,fill=Direction))+geom_col(width=.72)+geom_hline(yintercept=0)+
  geom_text(aes(y=label_y,label=Star),fontface="bold",size=2.7)+facet_wrap(~Region,ncol=2,scales="free_y")+
  scale_fill_manual(values=c(Up="skyblue",Down="yellowgreen"))+
  labs(x=NULL,y="Signed fold change: Braak High vs Low",fill=NULL)+theme_minimal(base_size=9)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6,face="bold"),strip.text=element_text(face="bold"),legend.position="top")
ggsave(file.path(out,"Figure4_regional_Braak_waterfall.png"),p,width=12.5,height=7,dpi=300)
ggsave(file.path(out,"Figure4_regional_Braak_waterfall.pdf"),p,width=12.5,height=7)
ggsave(file.path(main,"Figure4_regional_Braak_waterfall.png"),p,width=12.5,height=7,dpi=300)
ggsave(file.path(main,"Figure4_regional_Braak_waterfall.pdf"),p,width=12.5,height=7)
writeLines(capture.output(sessionInfo()),file.path(out,"Figure4_sessionInfo.txt"))
message("Saved Figure 4; stars indicate within-region BH-adjusted ANCOVA significance.")
