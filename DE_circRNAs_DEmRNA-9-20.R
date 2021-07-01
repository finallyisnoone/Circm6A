
#####DE circRNA and DE mRNA result from 8-28######

# model and Circm6A pos
circ.feature.df = read_delim("/data1/yeying/m6a_circ/final_plot/figure3/totalCirc.modelPred.autocirc.txt",delim = "\t")
double.pos.m6Acirc = filter(circ.feature.df,Circm6A=="positive_in_model")$id

mRNA_res.df
final_all_expres_res
# select DE CircRNAs
DEcirc.m6a.expres.df = left_join(filter(final_all_expres_res,expreStatus!="nonSig"),
                                 final_all_meth_res,by="id",suffix = c("_express","_meth"))%>%
  mutate(methStatus2 = case_when( methStatus=="hyper" & id%in%double.pos.m6Acirc ~  "hyper" ,
                                  methStatus=="hypo" & id%in%double.pos.m6Acirc ~  "hypo",
                                  !id%in%double.pos.m6Acirc ~  "non-m6A" , 
                                  TRUE ~ "nonSig" ))
table(final_all_expres_res$expreStatus)
DEcirc.m6a.expres.tmp = dplyr::select(DEcirc.m6a.expres.df,id,methStatus,methStatus2,type_express)
DEcirc.m6a.expres.tmp$methStatus[is.na(DEcirc.m6a.expres.tmp$methStatus)] = "non-m6A"
table(DEcirc.m6a.expres.tmp$methStatus2)
xtabs(~methStatus2+expreStatus,data = DEcirc.m6a.expres.df)

####################








##load cor data
setwd("/data1/yeying/m6a_circ/pancreatic/co-network/DEcircRNA_DEmRNA")

cor_cutoff= 0.5

mRNA_m6Astatus.lst = read_delim("mRNA_m6Astatus.type.lst",delim = "\t",col_names = c("id","mRNA_m6A"))
DEcircRNA_DEmRNA.cor.df = read_delim("DEcircRNA_DEmRNAs.cor.df",delim = "\t")
DEcircRNA_DEmRNA.cor.df = left_join(DEcircRNA_DEmRNA.cor.df,DEcirc.m6a.expres.tmp,by=c("circRNA_id"="id"))%>%
  left_join(mRNA_m6Astatus.lst,by=c("mRNA_id"="id"))%>%
  mutate(cor_status = case_when(cor_value >= cor_cutoff ~ "positive",
                                cor_value <= -cor_cutoff ~ "negative",
                                TRUE ~ "nonSig"),
         mRNA_expres = case_when(mRNA_id%in%up_genes$ensemble_id ~ "up",
                                 mRNA_id%in%down_genes$ensemble_id ~ "down",
                                 TRUE ~ "nonSig"),
         circRNA_m6A = ifelse(methStatus2=="non-m6A","non-m6A","m6A"),
         circRNA_m6A_raw = ifelse(methStatus=="non-m6A","non-m6A","m6A"),
         co_meth_status = case_when(
           circRNA_m6A == "m6A" & mRNA_m6A == "m6A" ~ "co-meth",
           circRNA_m6A == "non-m6A" & mRNA_m6A == "m6A" ~ "linear_meth_uniq",
           circRNA_m6A == "m6A" & mRNA_m6A == "non-m6A" ~ "circ_meth_uniq",
           circRNA_m6A == "non-m6A" & mRNA_m6A == "non-m6A" ~ "co-nonMeth"
         ))

DEcircRNA_DEmRNA.cor.df$co_meth_status = factor(DEcircRNA_DEmRNA.cor.df$co_meth_status,levels = c("co-meth","circ_meth_uniq","linear_meth_uniq","co-nonMeth"))
dim(DEcircRNA_DEmRNA.cor.df)
summary(DEcircRNA_DEmRNA.cor.df$cor_value)

DEcircRNA_DEmRNAs.cor05.df = filter(DEcircRNA_DEmRNA.cor.df,!grepl("nonSig",cor_status))

#write_delim(DEcircRNA_DEmRNAs.cor05.df,path = "DEcircRNA_DEmRNAs.cor05.df",delim = "\t")

xtabs(~circRNA_m6A+mRNA_m6A,DEcircRNA_DEmRNAs.cor05.df)
xtabs(~circRNA_m6A_raw+mRNA_m6A,DEcircRNA_DEmRNAs.cor05.df)
fisher.test(xtabs(~circRNA_m6A_raw+mRNA_m6A,DEcircRNA_DEmRNAs.cor05.df))

####stats of cor de circRNAs an de mRNAs####
#cor de mRNA with m6A-circRNA
table(DEcircRNA_DEmRNAs.cor05.df$co_meth_status)
cor_de_mRNAs = unique(filter(DEcircRNA_DEmRNAs.cor05.df,co_meth_status=="co-meth")$mRNA_id)
cor_de_mRNAs.df = filter(mRNA_res.df,ensemble_id%in%cor_de_mRNAs)
table(cor_de_mRNAs.df$expre_status) # up:623 ; down:182 
#cor de circRNAs
cor_de_circRNAs = unique(filter(DEcircRNA_DEmRNAs.cor05.df,co_meth_status=="co-meth")$circRNA_id)
cor_de_circRNAs.feature.df = filter(circ_feature.TIS,id%in%cor_de_circRNAs)%>%inner_join(DEcirc.m6a.expres.df)
xtabs(~methStatus2+expreStatus,cor_de_circRNAs.feature.df) # expres :up 107 ; down 5
length(cor_de_mRNAs);length(cor_de_circRNAs) # 874 ; 112

##############m6A affect circRNA-demRNA network#############
######all cor coefficient
DEcirc_DEmRNA_cor.hist <- ggplot(DEcircRNA_DEmRNA.cor.df,aes(x=cor_value, y=..density..))+theme_classic()+
  geom_histogram(fill = ggcolors[2],color = "white")+
  geom_density(size=1)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+gg_theme
######group by circ methylation status
wilcox_res <- wilcox.test(filter(DEcircRNA_DEmRNA.cor.df , circRNA_m6A != "non-m6A")$cor_value, 
                  filter(DEcircRNA_DEmRNA.cor.df , circRNA_m6A == "non-m6A")$cor_value,alternative = "greater")
wilcox_res$p.value
DEcirc_DEmRNA_cor.ecdf = ECDF_plot(df = DEcircRNA_DEmRNA.cor.df,value_var = "cor_value" , group_var = "circRNA_m6A" ,test_result = wilcox_res ,
                                  plot_title = "")+
  scale_color_manual(values=c("#e16428","#2980b9"))#c("#DE8E43","#364D54")) ### color config manually

DEcirc_DEmRNA_cor.ecdf

ggplot(DEcircRNA_DEmRNA.cor.df,aes(x=cor_value,color=co_meth_status))+theme_test()+
  stat_ecdf(size = 1)+theme(legend.position=c(0.75,0.25))+
  #annotate("text",x=-Inf,y=Inf,vjust=1.5,hjust=-.12,label=test_anno)+ 
  scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+ #,limits = c(-0.5,0.5)
  labs(y="Cumulative fraction")+scale_color_manual(values = c("#a8026f",ggcolors[4],ggcolors[3],"#05004e","#974949"))+
  theme(plot.title = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"))

ggdensity(DEcircRNA_DEmRNA.cor.df,x = "cor_value",color = "co_meth_status")

DEcirc_DEmRNA_cor.ecdf


##########enrichment : cor de mRNAs #############
#gsea
#cor de mRNA list 
cor_de_mRNAs.fc.lst = cor_de_mRNAs.df$log2FoldChange
names(cor_de_mRNAs.fc.lst) = cor_de_mRNAs.df$symbol
#calcu gsea
gsea.kegg.res = gsea_R(db_name = "/data1/yeying/database/human/KEGG_2019_Human.txt.new",total_genes.lst = cor_de_mRNAs.fc.lst)
grid.draw(gsea.kegg.res$top_pathway)
grid.draw(gsea.kegg.res$plot_table)
gsea.panther.res = gsea_R(db_name = "/data1/yeying/database/human/Panther_2016.txt.new",total_genes.lst = cor_de_mRNAs.fc.lst)
gsea.reactome.res = gsea_R(db_name = "/data1/yeying/database/human/Reactome_2016.txt.new",total_genes.lst = cor_de_mRNAs.fc.lst)

#enrichment
GoKegg(gene_list = as.character(filter(cor_de_mRNAs.df,expre_status=="up")$symbol),outdir = "./",prefix ="cor_up_mRNAs_9-23")
GoKegg(gene_list = as.character(filter(cor_de_mRNAs.df,expre_status=="down")$symbol),outdir = "./",prefix ="cor_down_mRNAs_9-23")



###############RBP ###################





##########8-28 : wer prefer to bind co-Meth DM-DE pairs############
#load data
circRNA_wer_bind.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/WER/circRNA_wer_bind.info.txt",delim = "\t")

mRNA_wer_bind.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/WER/mRNA_wer_bind.info.txt",delim = "\t")

##add wer bind data
dim(DMcirc_DEcircRNA.allsample.cor.df)
DMcirc_DEcircRNA.cor.allsample.wer_bind.df = left_join(DMcirc_DEcircRNA.allsample.cor.df,
                                                       dplyr::select(circRNA_wer_bind.df,-(chr),-(start),-(end),-(m6Astatus),-(strand)),by = c("circRNA_id" = "id"))%>%
  left_join(dplyr::select(mRNA_wer_bind.df,-(chr),-(start),-(end),-(m6Astatus),-(strand)),by = c("mRNA_id" = "id"),suffix = c("_circ","_mRNA"))%>%
  mutate(co_meth_status = case_when(mRNA_m6A=="m6A" & circRNA_m6A=="m6A" ~ "co-Meth",
                                    circRNA_m6A=="m6A" & mRNA_m6A=="non-m6A" ~ "circRNA-meth-only"))%>%dplyr::select(-(log2FoldChange))

dim(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)#1226544
length(colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)[grep("_circ",colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df))])
colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)[grep("_mRNA",colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df))]
table(DMcirc_DEcircRNA.cor.allsample.wer_bind.df$co_meth_status)
# write_delim(DMcirc_DEcircRNA.cor.allsample.wer_bind.df,
#             path = "/data1/yeying/m6a_circ/pancreatic/co-network/WER/DMcirc_DEcircRNA.cor.allsample.wer_bind.df",delim = "\t")

##add co bind status

DMcirc_DEcircRNA.co_bind.df = dplyr::select(DMcirc_DEcircRNA.cor.allsample.wer_bind.df,-contains("_circ"),-contains("_mRNA")) 
colnames(DMcirc_DEcircRNA.co_bind.df)
#loop
for (i in 1:17){
  #i=1
  tmp.wer = gsub("_circ","",colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)[(i+10)])
  tmp.wer.bind.df = DMcirc_DEcircRNA.cor.allsample.wer_bind.df[,c(1,4,(i+10),(i+10+17))]
  colnames(tmp.wer.bind.df) = c("mRNA_id","circRNA_id","circ_bind","mRNA_bind")
  tmp.wer.bind.df =  mutate(tmp.wer.bind.df,co_bind_status = case_when( 
    circ_bind > 0 & mRNA_bind > 0 ~ "co-bind",
    circ_bind > 0 & mRNA_bind == 0 ~ "circRNA-bind-only",
    circ_bind == 0 & mRNA_bind > 0 ~ "mRNA-bind-only",
    circ_bind == 0 & mRNA_bind == 0 ~ "co-non-bind"))%>% dplyr::select(-(circ_bind),-(mRNA_bind))
  tmp.wer.bind.df$co_bind_status = factor(tmp.wer.bind.df$co_bind_status,levels = c("co-bind","circRNA-bind-only","mRNA-bind-only","co-non-bind"))
  DMcirc_DEcircRNA.co_bind.df = inner_join(DMcirc_DEcircRNA.co_bind.df,tmp.wer.bind.df)
  colnames(DMcirc_DEcircRNA.co_bind.df)[i+11] = tmp.wer
  print(tmp.wer)
}

colnames(DMcirc_DEcircRNA.co_bind.df)






#save.image()





