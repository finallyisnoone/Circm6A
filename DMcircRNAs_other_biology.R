options(stringsAsFactors=F)

setwd("/data1/yeying/m6a_circ/pancreatic/DM")

###########m6A circRNA expression level##########
dim(circ_input_RPM_feature)
#allsample_circRNA_edgeR not with specific
#
circ_input_RPM_feature.express.df = left_join(final_all_expres_res,circ_input_RPM_feature,by="id")%>%
  dplyr::select(-c(colnames(circ_input_RPM_feature)[2:7]))%>%
  mutate(DM_status = case_when(id%in%hyper_circ ~ "hyper",
                               id%in%hypo_circ ~ "hypo" ,
                               id%in%m6Acirc.lst ~ "m6A" ,
                               TRUE ~ "non-m6A"))%>%filter(DM_status%in%c("hyper","hypo"))

xtabs(~ DM_status + expreStatus , circ_input_RPM_feature.express.df)
write_delim(dplyr::select(circ_input_RPM_feature.express.df,-contains("S_")),path = "/data1/yeying/m6a_circ/pancreatic/DM/biogenesis/DMcirc_feature.express.df",delim = "\t" )
###########m6A promote circRNA formation##########
circ_ratio.mat <- read_delim("circ_ratio/circ_ratio.filter.final.mat",delim = "\t")
circ_ratio.mat <- read_delim("circ_ratio/circ_ratio_dev8.final.mat",delim = "\t")
m6Acirc.lst[grep("chrM",m6Acirc.lst)]
dim(circ_ratio.mat)

circ_ratio.mat <- circ_ratio.mat%>%mutate(DM_status = case_when(id%in%hyper_circ ~ "hyper",
                                                                id%in%hypo_circ ~ "hypo" ,
                                                                id%in%m6Acirc.lst ~ "m6A" ,
                                                                TRUE ~ "non-m6A"))
table(circ_ratio.mat$DM_status)
circ_ratio.mat[1:5,1:5]
#circ_ratio.mat[circ_ratio.mat==0] = NA
length(intersect(circ_ratio.mat$id,c(hyper_circ,hypo_circ)))
circ.circRatio.mean.df = data.frame(id = circ_ratio.mat$id ,DM_status = circ_ratio.mat$DM_status, 
                                      tumor_ratio = rowMeans(circ_ratio.mat[,2:54],na.rm = T),
                                      normal_ratio = rowMeans(circ_ratio.mat[,55:78],na.rm = T))

###extract DM circRNAs
DMcirc.circRatio.mat = filter(circ_ratio.mat,id%in%c(hyper_circ,hypo_circ))%>%mutate(DM_status = ifelse(id%in%hyper_circ,"hyper","hypo"))
#DMcirc.circRatio.mean.df
dim(DMcirc.circRatio.mat)

#hyper circRNA : tumor vs normal

hypercirc.circRatio.mean.df = melt(filter(circ.circRatio.mean.df,DM_status=="hyper"),
                                 id.vars = colnames(circ.circRatio.mean.df)[1:2])
head(hypercirc.circRatio.mean.df)

ggecdf(data = hypercirc.circRatio.mean.df,x="value",color = "variable",palette = ggcolors[1:2],size = 1,
       ggtheme = theme(legend.spacing.x =unit(2,'mm'),legend.position = c(0.7,0.2)))+
  scale_y_continuous(expand = c(0.001,0.001),limits = c(0,1))+ylab("CDF")+
  scale_x_continuous(expand = c(0.001,0.001))+gg_theme


#hyper vs non-m6A circRNAs

hyper_nonm6A.circRatio.mean.df = melt(filter(circ.circRatio.mean.df,DM_status%in%c("hyper","non-m6A")),
                                   id.vars = colnames(circ.circRatio.mean.df)[1:2])

wilcox.test(filter(hyper_nonm6A.circRatio.mean.df,DM_status=="hyper")$value,
            filter(hyper_nonm6A.circRatio.mean.df,DM_status=="non-m6A")$value,alternative = "great")
table(hyper_nonm6A.circRatio.mean.df$DM_status)
ggecdf(data = hyper_nonm6A.circRatio.mean.df,x="value",color = "DM_status",palette = ggcolors[1:2],size = 1,
       ggtheme = theme(legend.spacing.x =unit(2,'mm'),legend.position = c(0.7,0.2)))+
  scale_y_continuous(expand = c(0.001,0.001),limits = c(0,1))+ylab("CDF")+
  scale_x_continuous(expand = c(0.000,0.000),limits = c(0,0.3))+gg_theme


