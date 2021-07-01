###DM-DE circRNAs
setwd("/data1/yeying/m6a_circ/pancreatic/DE_DM/9-23")

#####hyper promote up-regulated####

final_m6acirc_meth_expres_res.df = left_join(final_all_meth_res,final_all_expres_res,
          by="id",suffix = c("_meth","_express"))%>%
  mutate(log2meth = log2((Tumor_mean+0.1)/(Normal_mean+0.1)),
    #log2meth = log2((Tumor_mean+1)/(Normal_mean+1)),
         diff_meth = Tumor_mean - Normal_mean,
         label = paste0(methStatus,"_",expreStatus))%>%
  filter(!grepl("chrM",id))
table(final_m6acirc_meth_expres_res.df$methStatus)
##extract circ
hyper.circ = filter(final_m6acirc_meth_expres_res.df,methStatus=="hyper")$id
hypo.circ = filter(final_m6acirc_meth_expres_res.df,methStatus=="hypo")$id

summary(abs(color.df$log2meth))
test = log2((color.df$Tumor_mean+0.1)/(color.df$Normal_mean+0.1))
summary(abs(test))
summary(abs(color.df$diff_meth))
dim(final_m6acirc_meth_expres_res.df)
table(final_m6acirc_meth_expres_res.df$methStatus)
xtabs(~expreStatus+methStatus,data = filter(final_m6acirc_meth_expres_res.df,methStatus!="non-m6A"))

#recalcu log2FC of expression                                                              
expression_count.mat <- as.data.frame(read_delim("/data1/yeying/m6a_circ/pancreatic/matrix/11_7/filter_input_count.mat",delim = "\t"))%>%filter(!grepl("chrM",id))
rownames(expression_count.mat) = as.character(expression_count.mat$id)
dim(expression_count.mat[,-1])
Circ_norm_edgeR=cpm(expression_count.mat[,-1])
m6Acirc_cpm.mat = Circ_norm_edgeR[final_m6acirc_meth_expres_res.df$id,]

m6Acirc.expre.res.df <- data.frame(id = rownames(m6Acirc_cpm.mat), tumor_cpm_mean=rowMeans(m6Acirc_cpm.mat[,grep("T",colnames(m6Acirc_cpm.mat))]),
                                   normal_cpm_mean=rowMeans(m6Acirc_cpm.mat[,grep("N",colnames(m6Acirc_cpm.mat))]),
                                   freq=rowSums(m6Acirc_cpm.mat > 0))%>%
  mutate(log2FC_express = log2((tumor_cpm_mean+0.1)/(normal_cpm_mean+0.1)))

#####merge log2fc of expression and m6A####
final_m6acirc_meth_expres_res.df2 = left_join(final_m6acirc_meth_expres_res.df,m6Acirc.expre.res.df)

###scatter plot##

#unique(final_m6acirc_meth_expres_res.df2$label)
color.df = filter(final_m6acirc_meth_expres_res.df2,!grepl("nonSig",label))
color.df$label = factor(color.df$label,levels = c("hyper_upRegu","hyper_downRegu","hypo_downRegu"))

ggplot()+geom_point(data = filter(final_m6acirc_meth_expres_res.df2,grepl("nonSig",label)),
                    aes_string(y= "log2FC_express" ,x="log2meth"),size = 1,color="grey50")+
  geom_point(data =color.df ,aes_string(y= "log2FC_express" ,x="log2meth",color = "label"),size = 1)+
  xlim(-max(abs(final_m6acirc_meth_expres_res.df2$log2meth)),max(abs(final_m6acirc_meth_expres_res.df2$log2meth)))+
  ylim(-max(abs(final_m6acirc_meth_expres_res.df2$log2FC_express)),max(abs(final_m6acirc_meth_expres_res.df2$log2FC_express)))+
  geom_vline(xintercept = lfc,linetype="dashed")+geom_vline(xintercept = -lfc,linetype="dashed")+
  geom_hline(yintercept = -lfc,linetype="dashed")+geom_hline(yintercept  = lfc,linetype="dashed")+
  scale_color_manual(values = ggcolors)+
  theme_test()+gg_theme


###extact specific type circRNAs
#hyper-up
hyper_up.circ = filter(final_m6acirc_meth_expres_res.df,methStatus=="hyper",expreStatus=="upRegu")
#hyper nonDE
hyper_noDE.circ = filter(final_m6acirc_meth_expres_res.df,methStatus=="hyper",expreStatus=="nonSig")
write_delim(hyper_noDE.circ,"/data1/yeying/m6a_circ/pancreatic/DM/hyper_noDE.circ.df",delim = "\t")
dim(hyper_up.circ) #659
dim(hyper_noDE.circ) # 391
write_delim(hyper_up.circ,path = "hyper_up.circ.txt",delim = "\t")

######RBP pos. cor with hyper circRNAs in Tumor####
reader.lst
wers_DE_res = filter(mRNA_res.df,symbol%in%final.wers)
wers_DE_res
#write_delim(wers_DE_res,"/data1/yeying/m6a_circ/pancreatic/co-network/only_tumor/wers_DE_res.df",delim = "\t")
#load total_circ
circRNA_wers.cor.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/only_tumor/circRNA_wers.cor.df",delim = "\t")%>%filter(!is.na(cor_value))
dim(circRNA_wers.cor.df)
ggecdf(circRNA_wers.cor.df,x="cor_value")
##hyper up cor wer
cor_cutoff = 0.5
hyper_upcirc.reader.cor.anno.df = filter(circRNA_wers.cor.df,circRNA_id%in%hyper_up.circ$id)%>%
  left_join(gene.info,res,by= c("mRNA_id" = "ensemble_id"))%>%filter(symbol%in%reader.lst$symbol)

max(hyper_upcirc.reader.cor.anno.df$cor_value)

hyper_up.reader.ecdf.plot.df = rbind(dplyr::select(circRNA_wers.cor.df,cor_value)%>%mutate(label="total circRNAs"),
                                     data.frame(cor_value = hyper_upcirc.reader.cor.anno.df$cor_value,
                                                label = "hyper-up CircRNAs"))#hyper_upcirc.reader.cor.anno.df$symbol))
#hyper_up.reader.ecdf.plot.df$label = factor(hyper_up.reader.ecdf.plot.df$label,levels = c("totalCirc"))


hyper_up.reader.density.plot <- ggplot(data = hyper_up.reader.ecdf.plot.df,aes(x=cor_value,color = label))+stat_density(geom = "line",position = "identity",
                      size = 1.5,adjust = 0.8)+
  theme_test()+scale_color_manual(values = c("grey50",ggcolors[1]))+
  scale_y_continuous(expand = c(0.001,0.001))+scale_x_continuous(expand = c(0,0))+
  labs( y="Cumulative fraction")+
  theme(
        axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
        legend.position=c(0.2,.8),
        legend.text = element_text(size = 15, face = NULL, color = "black"),
        legend.title= element_blank())+geom_vline(xintercept = 0.5,linetype="dashed")
hyper_up.reader.density.plot

##########DM translation################
setwd("/data1/yeying/m6a_circ/pancreatic/DM/translation")

circRNA_TIS_eIFs.anno.df = read_delim("/data1/yeying/m6a_circ/pancreatic/DM/translation/circRNA_TIS_eIFs.anno.xls",delim = "\t")
colnames(circRNA_TIS_eIFs.anno.df)
dim(circRNA_TIS_eIFs.anno.df)
#hyper
dim(filter(final_m6acirc_meth_expres_res.df,methStatus=="hyper"))
hyper.circ = filter(final_m6acirc_meth_expres_res.df,methStatus=="hyper")$id
hypercircRNA_TIS_eIFs.anno.df = filter(circRNA_TIS_eIFs.anno.df,id%in%hyper.circ)
dim(hypercircRNA_TIS_eIFs.anno.df)
hypercircRNA_TIS_eIFs.pos.df = hypercircRNA_TIS_eIFs.anno.df[rowSums(hypercircRNA_TIS_eIFs.anno.df[,c(3,6)]) > 0,]
dim(hypercircRNA_TIS_eIFs.pos.df)

##hyper non-DE
hyper_nonDE_circRNA_TIS_eIFs.anno.df = filter(circRNA_TIS_eIFs.anno.df,id%in%hyper_noDE.circ$id)
hyper_nonDE_circRNA_TIS_eIFs.pos.df = hyper_nonDE_circRNA_TIS_eIFs.anno.df[rowSums(hyper_nonDE_circRNA_TIS_eIFs.anno.df[,c(3,6)]) > 0,]
dim(hyper_nonDE_circRNA_TIS_eIFs.pos.df)

###p-value = 0.002927
fisher.test(matrix(c(dim(hyper_nonDE_circRNA_TIS_eIFs.pos.df),length(hyper_noDE.circ$id),
                     dim(hypercircRNA_TIS_eIFs.pos.df),length(hyper.circ)),nrow = 2,byrow = T),alternative = "greater")

tmp = hyper_nonDE_circRNA_TIS_eIFs.anno.df%>%filter(symbol%in%oncogene.df$symbol)
unique(tmp$symbol)
dim(tmp)

