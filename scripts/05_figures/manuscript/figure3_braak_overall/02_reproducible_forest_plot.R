suppressPackageStartupMessages({library(here);library(dplyr);library(purrr);library(broom);library(ggplot2)})
x <- readRDS(here("results","figures","intermediate","braak","braak_analysis_objects.rds"))
out <- here("results","figures","intermediate","figure3");dir.create(out,recursive=TRUE,showWarnings=FALSE)
genes <- colnames(x$expr_combat_z)
res <- map_dfr(genes,function(g){
  dat <- x$df |> mutate(Expression=x$expr_combat_z[SampleID,g])
  z <- tidy(glm(Braak_bin~Expression+REGION+AGE+SEX,data=dat,family=binomial()),conf.int=TRUE,exponentiate=TRUE) |> filter(term=="Expression")
  transmute(z,Gene=g,OR=estimate,CI_low=conf.low,CI_high=conf.high,p_value=p.value)
}) |> mutate(FDR=p.adjust(p_value,"BH"),Direction=ifelse(OR>=1,"Higher odds","Lower odds")) |> arrange(FDR)
write.csv(res,file.path(out,"Figure3B_Braak_logistic_regression.csv"),row.names=FALSE)
d <- res |> filter(is.finite(OR),is.finite(CI_low),is.finite(CI_high))
d$Gene <- factor(d$Gene, levels=d$Gene[order(d$OR)])
p <- ggplot(d,aes(OR,Gene,color=Direction))+geom_vline(xintercept=1,linetype=2,color="grey40")+
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.2)+geom_point(size=2.5)+
  scale_x_log10()+scale_color_manual(values=c("Higher odds"="tomato3","Lower odds"="steelblue3"))+
  labs(x="Odds ratio per 1-SD expression (95% CI)",y=NULL,color=NULL)+theme_minimal(base_size=10)+
  theme(axis.text.y=element_text(face="bold"),legend.position="top")
# No significance stars here by design; stars are reserved for waterfall plots.
ggsave(file.path(out,"Figure3B_Braak_forest.png"),p,width=8,height=7,dpi=300)
ggsave(file.path(out,"Figure3B_Braak_forest.pdf"),p,width=8,height=7)
message("Saved Figure 3B forest plot without significance stars.")
