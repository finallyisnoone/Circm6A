options(stringsAsFactors=F)
setwd("/data1/yeying/m6a_circ/pancreatic/8_12/DE_mRNA/")

library(DESeq2)

#####DE mRNAs#####
raw.total_RNA.mat = read.table("gencode.rsem.count.txt",sep = "\t")
raw.total_RNA.mat$mRNA_ensemble_id = as.character(rownames(raw.total_RNA.mat))
colnames(raw.total_RNA.mat) = gsub("X","S_",colnames(raw.total_RNA.mat))
mRNA_fpkm.mat = read_delim("mRNA_FPKM.3samples.mat",delim = "\t")
colnames(mRNA_fpkm.mat)=gsub("_linear","",colnames(mRNA_fpkm.mat))

##filter mRNA
mRNA.count.mat = as.data.frame(raw.total_RNA.mat[,colnames(mRNA_fpkm.mat)[-2]]%>%filter(mRNA_ensemble_id%in%mRNA_fpkm.mat$mRNA_ensemble_id))
rownames(mRNA.count.mat) = mRNA.count.mat$mRNA_ensemble_id
mRNA.count.mat2 = mRNA.count.mat[,-1]
dim(mRNA.count.mat2)
#####DESeq2
fc = 1.5
lfc = log2(fc)

keep <- rowSums(mRNA.count.mat2 > 0) >= 3 #a Count>0 in at least 3 samples
countData <- mRNA.count.mat2[keep,]
dim(countData)
pheno_data = get_pheno(x = colnames(countData),label1 = "T",label2 = "N",group1 = "Tumor",group2 = "Normal")

dds = DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~ Type)
dds = DESeq(dds)
res = results(dds, contrast=c("Type","Tumor","Normal"))
res = as.data.frame(res)
res$mRNA_ensemble_id = rownames(res)

##anno
gene.info = read_delim("/data1/yeying/database/human/hg38/Gencode/hg38_gencode.v25_region_ensembleTosymble_map.info",
                       delim = "\t")
head(gene.info)
mRNA_res.df =  right_join(gene.info,res,by= c( "ensemble_id" = "mRNA_ensemble_id" ))%>%
  mutate(expre_status = case_when( log2FoldChange >= lfc &  padj <= p_value ~ "up" ,
                                   log2FoldChange <= -lfc &  padj <= p_value ~ "down",
                                   TRUE ~ "non_diff"))%>%arrange(desc(log2FoldChange))
# 
# mutate(expre_status = case_when( log2FoldChange >= lfc &  pvalue <= p_value ~ "up" ,
#                                                  log2FoldChange <= -lfc &  pvalue <= p_value ~ "down",
#                                                  TRUE ~ "non_diff"))%>%arrange(desc(log2FoldChange))


table(mRNA_res.df$expre_status)

up_genes =  filter(mRNA_res.df,expre_status=="up")%>%dplyr::select(symbol,ensemble_id,log2FoldChange)%>%arrange(desc(log2FoldChange))

circ.factor = c("ADAR1","QKI","DHX9")

circ.factor.de_res.df = filter(mRNA_res.df,symbol%in%circ.factor)
filter(mRNA_res.df,grepl("ADAR",symbol))
circ.factor.de_res.df
down_genes =  filter(mRNA_res.df,expre_status=="down")%>%dplyr::select(symbol,ensemble_id,log2FoldChange)%>%arrange(desc(log2FoldChange))


tmp = dplyr::select(mRNA_res.df,ensemble_id,padj,log2FoldChange)
colnames(tmp) = c("id","FDR","logFC")
volcano_plot(tmp,contrast_factor=rev(c("tumor","normal")),ylab_variable = "FDR",pval = 0.05,fc = 1.5)

#####sup : volcano plot for DEmRNAs####
source("/data1/yeying/m6a_circ/script/R_function.R")

volcano.plot( DE_res.df = dplyr::select(mRNA_res.df,ensemble_id,padj,log2FoldChange))

#write_delim(mRNA_res.df,path = "/data1/yeying/m6a_circ/pancreatic/8_12/DE_mRNA/DE_mRNA_result.xls",delim = "\t")
#write_delim(filter(mRNA_res.df,expre_status%in%c("up","down")),path = "/data1/yeying/m6a_circ/pancreatic/8_12/DE_mRNA/DE_mRNA_list.xls",delim = "\t")


###########pearson : new correct correlation (only in tumor samples)##############
setwd("/data1/yeying/m6a_circ/pancreatic/co-network/only_tumor")
cor_cutoff = 0.5
DMcircRNA_DEmRNAs.cor.df = read_delim("DMcircRNA_DEmRNAs.cor.nonNA.df",delim = "\t")
#all_circ_mRNA.m6Astatus.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/allcircRNA_mRNA.m6Astatus.df",delim = "\t")
#length(unique(DMcircRNA_DEmRNAs.cor.df$circRNA_id))
head(DMcircRNA_DEmRNAs.cor.df)
summary(DMcircRNA_DEmRNAs.cor.df$cor_value)

ggplot(DMcircRNA_DEmRNAs.cor.df,aes(x=cor_value, y=..density..))+theme_classic()+
  geom_histogram(fill = ggcolors[2],color = "white")+
  geom_density(size=1)+xlab("Spearman Rank Correlation")+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+gg_theme

### filter correlated circRNA-mRNAs pairs ####
#filter with correlation cut-off 
DMcircRNA_DEmRNAs.tumor.cor.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/only_tumor/DMcircRNA_DEmRNAs.cor.nonNA.df",delim = "\t")%>%
  mutate(cor_status = case_when(cor_value >= cor_cutoff ~ "positive",
                                cor_value <= -cor_cutoff ~ "negative",
                                TRUE ~ "nonSig"),
         pairs_id = paste0(mRNA_id,"_",circRNA_id))

table(DMcircRNA_DEmRNAs.tumor.cor.df$cor_status)
dim(DMcircRNA_DEmRNAs.tumor.cor.df)
##DMcircRNA_DEmRNAs.cor.df : DMcircRNA_DEmRNAs.tumor.cor.df > cor_cutoff
DMcircRNA_DEmRNAs.cor.df = filter(DMcircRNA_DEmRNAs.tumor.cor.df,!grepl("nonSig",cor_status))
dim(DMcircRNA_DEmRNAs.cor.df)
table(DMcircRNA_DEmRNAs.cor.df$cor_status)

#cor de mRNAs
cor_de_mRNAs = unique(filter(DMcircRNA_DEmRNAs.cor.df,!grepl("nonSig",cor_status))$mRNA_id)
cor_de_mRNAs.df = filter(mRNA_res.df,ensemble_id%in%cor_de_mRNAs)
table(cor_de_mRNAs.df$expre_status) # up:551 ; down:325 
#cor dm circRNAs
cor_dm_circRNAs = unique(filter(DMcircRNA_DEmRNAs.cor.df,!grepl("nonSig",cor_status))$circRNA_id)
cor_dm_circRNAs.feature.df = filter(circ_feature.TIS,id%in%cor_dm_circRNAs)%>%inner_join(final_all_meth_res)
table(cor_dm_circRNAs.feature.df$methStatus)
length(cor_de_mRNAs);length(cor_dm_circRNAs) # 876 ; 127
length(unique(cor_dm_circRNAs.feature.df$symbol)) # 92

##########enrichment : cor de mRNAs #############
#gsea
#cor de mRNA list 
cor_de_mRNAs.fc.lst = cor_de_mRNAs.df$log2FoldChange
names(cor_de_mRNAs.fc.lst) = cor_de_mRNAs.df$symbol
#calcu gsea
# gsea.kegg.res = gsea_R(db_name = "/data1/yeying/database/human/KEGG_2019_Human.txt.new",total_genes.lst = cor_de_mRNAs.fc.lst)
# grid.draw(gsea.kegg.res$top_pathway)
# grid.draw(gsea.kegg.res$plot_table)
# gsea.panther.res = gsea_R(db_name = "/data1/yeying/database/human/Panther_2016.txt.new",total_genes.lst = cor_de_mRNAs.fc.lst)
# gsea.reactome.res = gsea_R(db_name = "/data1/yeying/database/human/Reactome_2016.txt.new",total_genes.lst = cor_de_mRNAs.fc.lst)

#enrichment
#GoKegg(gene_list = as.character(filter(cor_de_mRNAs.df,expre_status=="up")$symbol),outdir = "./",prefix ="_cor_up_mRNAs")
#GoKegg(gene_list = as.character(filter(cor_de_mRNAs.df,expre_status=="down")$symbol),outdir = "./",prefix ="_cor_down_mRNAs")
#write_delim(cor_de_mRNAs.df,path = "DMcircRNA_cor_de_mRNAs.df",delim = "\t")

####enrichr#####

enrichr_barplot <- function(res.df,top_number=10,bar_color = ggcolors[5]){ 
  ggbarplot(data = res.df[1:top_number,],x = "Description",y = "log10qvalue",
            color = "white",fill = bar_color,sort.val = "asc",orientation = "horiz")+
    xlab("")+ylab("-log10(qValue)")+scale_y_continuous(expand = c(0,0))+theme(
      axis.title = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.spacing.x = unit(0.2,"cm"),
      legend.title= element_blank(),
      legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
      legend.background = element_rect(fill = NA))
}

pathway_dir = "up_mRNA_enrichr/"
pathway_dir = "down_mRNA_enrichr/"
pathway.plot.lst = list()
list.files(pathway_dir)
for (i in 1:4){
  #i=list.files(pathway_dir)[1]
  path.name = strsplit(list.files(pathway_dir)[i],"_")[[1]][1]
  
  pathway.res.df = read_delim(paste0(pathway_dir,list.files(pathway_dir)[i]),delim = "\t")%>%
    mutate(Description=Term,log10qvalue=-log10(`Adjusted P-value`) )%>%filter(log10qvalue>= -log10(0.05))
  
  if (nrow(pathway.res.df) > 10) {
    pathway.res.df = pathway.res.df[1:10,]
  }
  p <- enrichr_barplot(pathway.res.df,nrow(pathway.res.df),bar_color = ggcolors[i])+ggtitle(path.name)
  pathway.plot.lst[[i]] = p
}

#pathway.plot.lst[[3]]

plot_grid(plotlist = pathway.plot.lst,ncol = 2,scale = c(1,1,1,1))

###########dm : cor dm circRNAs ###############
#cor dm circRNA list
# GoKegg(gene_list = as.character(unique(filter(cor_dm_circRNAs.feature.df,type=="tumor_uniq")$symbol)),outdir = "./",prefix ="_cor_hyper_circRNAs")
# GoKegg(gene_list = as.character(unique(filter(cor_dm_circRNAs.feature.df,type=="normal_uniq")$symbol)),outdir = "./",prefix ="_cor_hypo_circRNAs")

DEmRNAs_DMcircRNAs.cor.mat = circRNA_mRNA.simple.cor.mat[DE_mRNAs.lst,intersect(DMcircRNAs.lst,colnames(circRNA_mRNA.simple.cor.mat))]

########### new correct correlation ( only in normal samples)##############

DMcircRNA_DEmRNAs.normal.cor.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/only_normal/DEmRNAs_DMcircRNAs.cor.anno.df",delim = "\t")
DMcircRNA_DEmRNAs.normal.cor.df = filter(DMcircRNA_DEmRNAs.normal.cor.df,!is.na(cor_value))%>%
  mutate(cor_status = case_when(cor_value >= cor_cutoff ~ "positive",
                                cor_value <= -cor_cutoff ~ "negative",
                                TRUE ~ "nonSig"),
         pairs_id = paste0(mRNA_id,"_",circRNA_id))

dim(DMcircRNA_DEmRNAs.normal.cor.df)
table(DMcircRNA_DEmRNAs.normal.cor.df$cor_status)


##############merge tumor and normal cor.df#############
cor_cutoff = 0.5
intersect(colnames(DMcircRNA_DEmRNAs.cor.df),colnames(DMcircRNA_DEmRNAs.normal.cor.df))

totalCirc_DEmRNA.allsample.cor.df=read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/DEmRNA/share_MRE/circRNA_DEmRNAs.allsample.cor.non_shareMER.df",
                                             delim = "\t")%>%
  left_join(dplyr::select(circ_feature,id,m6Astatus),by = c("circRNA_id" = "id"))%>%
  mutate(cor_status_tumor = case_when(cor_value_tumor >= cor_cutoff ~ "positive",
                                      cor_value_tumor <= -cor_cutoff ~ "negative",
                                      TRUE ~ "nonSig"),
         cor_status_normal= case_when(cor_value_normal >= cor_cutoff ~ "positive",
                                      cor_value_normal <= -cor_cutoff ~ "negative",
                                      TRUE ~ "nonSig"),
         pairs_id = paste0(mRNA_id,"_",circRNA_id))%>%
  filter(!circRNA_id%in%ambiguou.circRNA.lst)

dim(totalCirc_DEmRNA.allsample.cor.df)
head(totalCirc_DEmRNA.allsample.cor.df$m6Astatus)
head(totalCirc_DEmRNA.allsample.cor.df$pairs_id)


# in tumor
m6ACirc_DEmRNA.tumor.stat = table(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A")$cor_status_tumor)
nonm6ACirc_DEmRNA.tumor.stat = table(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="non-m6A")$cor_status_tumor)

# in normal
m6ACirc_DEmRNA.normal.stat = table(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A")$cor_status_normal)
nonm6ACirc_DEmRNA.normal.stat = table(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="non-m6A")$cor_status_normal)

####### random sample 100k pairs########
set.seed(2019)
sample100k.totalCirc_DEmRNA.allsample.cor.df = sample_n(totalCirc_DEmRNA.allsample.cor.df,size = 100000)%>%
  left_join(mRNA_m6Astatus.lst,by = c("mRNA_id"="id"))%>%
  mutate(meth_type = case_when(
    m6Astatus == "m6A" & mRNA_m6A == "m6A" ~ "co-meth",
    m6Astatus == "non-m6A" & mRNA_m6A == "m6A" ~ "linear_meth_uniq",
    m6Astatus == "m6A" & mRNA_m6A == "non-m6A" ~ "circ_meth_uniq",
    m6Astatus == "non-m6A" & mRNA_m6A == "non-m6A" ~ "co-nonMeth"))

########comparison of cor value in tumor and normal###########
setwd("/data1/yeying/m6a_circ/pancreatic/co-network/DMcircRNA_DEmRNA")
# m6A vs non-m6A

m6A_nonm6A.wilcox.res = wilcox.test(filter(sample100k.totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A")$cor_value_tumor,
            filter(sample100k.totalCirc_DEmRNA.allsample.cor.df,m6Astatus!="m6A")$cor_value_tumor,alternative = "greater")
m6A_nonm6A.wilcox.res
# hyper.m6A.wilcox.res = wilcox.test(filter(totalCirc_DEmRNA.allsample.cor.df,circRNA_id%in%hyper.circ)$cor_value_tumor,
#                                     filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A")$cor_value_tumor,alternative = "greater")

#hyper.m6A.wilcox.res
table(sample100k.totalCirc_DEmRNA.allsample.cor.df$m6Astatus)
m6a_nonm6a_cor.plot <- 
  ggplot(sample100k.totalCirc_DEmRNA.allsample.cor.df)+
  geom_density(aes(x=cor_value_tumor, y=..density..,color=m6Astatus),size=1,adjust = 2.5)+xlab("Spearman Rank Correlation")+
  # geom_density(data = filter(totalCirc_DEmRNA.allsample.cor.df,circRNA_id%in%hyper.circ),
  #              aes(x= cor_value_tumor),size=1,color=ggcolors[3],adjust = 2.5)+
  theme_classic()+
  scale_color_manual(values = ggcolors)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  geom_vline(xintercept = cor_cutoff)+
  ggtitle("Random co-expression pairs(n = 100,000)")+
  annotate("rect",xmin = cor_cutoff,xmax = 1,alpha=0.2,fill=ggcolors[3], ymin = 0,ymax = 3)+
  geom_vline(xintercept = -cor_cutoff)+
  annotate("rect",xmin = -cor_cutoff,xmax = -1,alpha=0.2,fill=ggcolors[3], ymin = 0,ymax = 3)+
  annotate("text",x= -0.75,y=2,
           label=paste0("co-expression \npairs\nwilcox test ,\np = ",signif(m6A_nonm6A.wilcox.res$p.value,digits = 3)))+
  annotate("text",x= 0.75,y=2,label=paste0("N = ",m6ACirc_DEmRNA.tumor.stat[1]+m6ACirc_DEmRNA.tumor.stat[3]),color = ggcolors[1])+
  annotate("text",x= 0.75,y=1.8,label=paste0("N = ",nonm6ACirc_DEmRNA.tumor.stat[1]+nonm6ACirc_DEmRNA.tumor.stat[3]),color = ggcolors[2])+gg_theme

m6a_nonm6a_cor.plot


#rm(m6A_nonm6A.wilcox.res)

#  tumor vs normal 
# tmp = data.frame(cor_value = c(sample100k.totalCirc_DEmRNA.allsample.cor.df$cor_value_tumor,sample100k.totalCirc_DEmRNA.allsample.cor.df$cor_value_normal),
#                  type=c(rep("tumor",nrow(sample100k.totalCirc_DEmRNA.allsample.cor.df)),rep("normal",nrow(sample100k.totalCirc_DEmRNA.allsample.cor.df))),
#                  cor_status = c(sample100k.totalCirc_DEmRNA.allsample.cor.df$cor_status_tumor,sample100k.totalCirc_DEmRNA.allsample.cor.df$cor_status_normal))
# 
# tmp$type = factor(tmp$type,levels = c("tumor","normal"))
# 
# 
# 
# tumor_normal_cor.plot <- ggplot(tmp,aes(x=cor_value, y=..density..,color=type))+theme_classic()+
#   geom_density(size=1,adjust = 2)+xlab("Spearman Rank Correlation")+
#   scale_color_manual(values = ggcolors)+
#   scale_y_continuous(expand = c(0, 0))+
#   scale_x_continuous(expand = c(0, 0))+
#   geom_vline(xintercept = cor_cutoff)+
#   annotate("rect",xmin = cor_cutoff,xmax = 1,alpha=0.2,fill=ggcolors[3], ymin = 0,ymax = 3)+
#   geom_vline(xintercept = -cor_cutoff)+
#   annotate("rect",xmin = -cor_cutoff,xmax = -1,alpha=0.2,fill=ggcolors[3], ymin = 0,ymax = 3)+
#   annotate("text",x= 0.75,y=2,label="co-expression \npairs")+gg_theme
# 
# tumor_normal_cor.plot

######### m6A vs non-m6A circRNA pairs ########

library("VennDiagram")
library(grDevices)

dim(totalCirc_DEmRNA.allsample.cor.df)

#total circRNA-DEmRNA
total.TvsN.venn = venn.diagram(x = list(tumor_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,abs(cor_value_tumor) >= cor_cutoff)$pairs_id),
                                      normal_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,abs(cor_value_normal) >= cor_cutoff)$pairs_id)),
                             filename = NULL,
                             main= "total circRNA co-expression pairs in Tumor and normal",
                             col= ggcolors[1:2],
                             fill = c(alpha(ggcolors[1],0.3), alpha(ggcolors[2],0.3)),
                             cex = 1.5,
                             category.names = c("Tumor" , "Normal"),
                             cat.pos = c(-180, 180),
                             cat.col = ggcolors[1:2])
grid.newpage()
grid.draw(total.TvsN.venn)


m6A.TvsN.lst  = list(tumor_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A",abs(cor_value_tumor) >= cor_cutoff)$pairs_id),
                normal_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A",abs(cor_value_normal) >= cor_cutoff)$pairs_id))

m6A.TvsN.venn = venn.diagram(x = m6A.TvsN.lst,filename = NULL,
                             main= "m6A circRNA co-expression pairs in Tumor and normal",
                             col= ggcolors[1:2],
                             fill = c(alpha(ggcolors[1],0.3), alpha(ggcolors[2],0.3)),
                             cex = 2,
                             category.names = c("Tumor" , "Normal"),
                             cat.pos = c(-180, 180),
                             cat.col = ggcolors[1:2])
grid.newpage()
grid.draw(m6A.TvsN.venn)

nonm6A.TvsN.lst  = list(tumor_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="non-m6A",abs(cor_value_tumor) >= cor_cutoff)$pairs_id),
                    normal_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="non-m6A",abs(cor_value_normal) >= cor_cutoff)$pairs_id))
nonm6A.TvsN.venn = venn.diagram(x = nonm6A.TvsN.lst,filename = NULL,
                                main= "non-m6A circRNA co-expression pairs in Tumor and normal",
                                col= ggcolors[1:2],
                                fill = c(alpha(ggcolors[1],0.3), alpha(ggcolors[2],0.3)),
                                cex = 1.5,
                                category.names = c("Tumor" , "Normal"),
                                cat.pos = c(-180, 180),
                                cat.col = ggcolors[1:2])


# pdf("/data1/yeying/m6a_circ/final_plot/figure6/co-expression_venn-10-11.pdf",height = 4,width = 4)
# grid.newpage()
# grid.draw(total.TvsN.venn)
# grid.newpage()
# grid.draw(m6A.TvsN.venn)
# grid.newpage()
# grid.draw(nonm6A.TvsN.venn)
# dev.off()


#########DM RNA pairs in tumor and normal : Venn ######
###########DMcircRNA - DEmRNA network#######

DMcirc_DEmRNA.allsample.cor.df0 = full_join(DMcircRNA_DEmRNAs.tumor.cor.df,DMcircRNA_DEmRNAs.normal.cor.df,
                                           by =c("mRNA_id","mRNA_express","circRNA_id","meth_status","pairs_id"),
                                           suffix=c("_tumor","_normal"))
head(totalCirc_DEmRNA.allsample.cor.df)
DMcirc_DEmRNA.allsample.cor.df = filter(totalCirc_DEmRNA.allsample.cor.df,circRNA_id%in%DM_list)%>%
  left_join(dplyr::select(DMcirc_DEmRNA.allsample.cor.df0,c("mRNA_id","mRNA_express","circRNA_id","meth_status")))

rm(DMcirc_DEmRNA.allsample.cor.df0)
dim(DMcirc_DEmRNA.allsample.cor.df)
head(DMcirc_DEmRNA.allsample.cor.df)
#hyper circRNAs
hypercirc_DEmRNA.allsample.cor.df = filter(DMcirc_DEmRNA.allsample.cor.df,circRNA_id%in%hyper.circ)

hyper.pairs.lst  = list(tumor_pairs = as.character(filter(DMcircRNA_DEmRNAs.tumor.cor.df,meth_status=="hyper",abs(cor_value)>=cor_cutoff)$pairs_id),
            normal_pairs = as.character(filter(DMcircRNA_DEmRNAs.normal.cor.df,meth_status=="hyper",abs(cor_value)>=cor_cutoff)$pairs_id))
 
hypo.pairs.lst  = list(tumor_pairs = as.character(filter(DMcircRNA_DEmRNAs.tumor.cor.df,meth_status=="hypo",abs(cor_value)>=cor_cutoff)$pairs_id),
                normal_pairs = as.character(filter(DMcircRNA_DEmRNAs.normal.cor.df,meth_status=="hypo",abs(cor_value)>=cor_cutoff)$pairs_id))

#venn.diagram(x = tmp.lst,filename = NULL,fill = ggcolors[1:2])
hyper.pairs.venn = venn.diagram(x = hyper.pairs.lst,filename = NULL,
               main= "hypo circRNA co-expression pairs in tumor and normal",
               col= ggcolors[1:2],
               fill = c(alpha(ggcolors[1],0.3), alpha(ggcolors[2],0.3)),
               cex = 2,
               category.names = c("Tumor" , "Normal"),
               cat.pos = c(-180, 180),
               cat.col = ggcolors[1:2])
grid.newpage()
grid.draw(hyper.pairs.venn)

xtabs(~ cor_status_normal + cor_status_tumor,data = hypercirc_DEmRNA.allsample.cor.df)

#hypo
hypocirc_DEmRNA.allsample.cor.df = filter(DMcirc_DEmRNA.allsample.cor.df,meth_status=="hypo")
xtabs(~ cor_status_tumor+ cor_status_normal ,data = hypocirc_DEmRNA.allsample.cor.df)


###########hyper/hypo circRNA contribute to co-expression network#########

normlzs.venn <- function(raw.lst){
  tumor.pairs = raw.lst$tumor_pairs
  normal.pairs = raw.lst$normal_pairs
  share.pairs = length(intersect(tumor.pairs,normal.pairs))
  tumor.pairs.gain = length(setdiff(tumor.pairs,normal.pairs))
  tumor.pairs.loss = length(setdiff(normal.pairs,tumor.pairs))
  venn.vec = c(share.pairs = share.pairs,pairs_gain = tumor.pairs.gain,pairs_loss = tumor.pairs.loss,
    pairs_gain_normlz =  tumor.pairs.gain/share.pairs , pairs_loss_normlz =  tumor.pairs.loss/share.pairs )
  return(venn.vec)
}

##total circRNA co-expression alter

total_circ.alter = normlzs.venn(list(tumor_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,abs(cor_value_tumor) >= cor_cutoff)$pairs_id),
                                     normal_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,abs(cor_value_normal) >= cor_cutoff)$pairs_id)))
round(total_circ.alter,digits  = 2)
dim(circ_feature)
sample1k.Circ.lst = sample(circ_feature$id,size = 1052)
sample1k.Circ_DEmRNA.allsample.cor.df = filter(totalCirc_DEmRNA.allsample.cor.df,circRNA_id%in%sample1k.Circ.lst)
sample1k_circ.alter = normlzs.venn(list(tumor_pairs = as.character(filter(sample1k.Circ_DEmRNA.allsample.cor.df,abs(cor_value_tumor) >= cor_cutoff)$pairs_id),
                     normal_pairs = as.character(filter(sample1k.Circ_DEmRNA.allsample.cor.df,abs(cor_value_normal) >= cor_cutoff)$pairs_id)))

m6acirc.alter = normlzs.venn(m6A.TvsN.lst)

nonDMcirc.pairs.lst = list(tumor_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,circRNA_id%in%nonDM.lst$id,abs(cor_value_tumor) >= cor_cutoff)$pairs_id),
                           normal_pairs = as.character(filter(totalCirc_DEmRNA.allsample.cor.df,circRNA_id%in%nonDM.lst$id,abs(cor_value_normal) >= cor_cutoff)$pairs_id))

nonDMcirc.alter = normlzs.venn(nonDMcirc.pairs.lst)

hypercirc.alter = normlzs.venn(hyper.pairs.lst)
hypocirc.alter = normlzs.venn(hypo.pairs.lst)

###hyper vs total
chisq.test(matrix(c(hypercirc.alter[4:5],total_circ.alter[4:5]),nrow = 2,byrow = T,dimnames = list(c("m6A","total"),c("tumor gain", "normal gain"))),correct = T)
###hyper vs nonDM
chisq.test(matrix(round(c(hypercirc.alter[4:5],nonDMcirc.alter[4:5])),nrow = 2,byrow = T,
                  dimnames = list(c("hyper","nonDM"),c("tumor gain", "normal gain"))),correct = T)
c(nonDMcirc.alter)

matrix(c(hypercirc.alter[4:5],nonDMcirc.alter[4:5]),nrow = 2,byrow = T,dimnames = list(c("hyper","nonDM"),c("tumor gain", "normal gain")))

##final
tmp = matrix(c(hypercirc.alter[2:3],nonDMcirc.alter[2:3]),nrow = 2,byrow = T,dimnames = list(c("hyper","nonDM"),c("tumor gain", "normal gain")))
tmp
chisq.test(tmp)

# hyper circRNAs vs non sig m6A circRNAs, permutation ########

# nonDM.lst = filter(final_all_meth_res,methStatus=="nonSig",!grepl("chrM",id))
# table(final_all_meth_res$methStatus)
# dim(nonDM.lst)
#write_delim(nonDM.lst,path = "/data1/yeying/m6a_circ/pancreatic/co-network/DEmRNA/nonDM_m6AcircRNA.df",delim = "\t")

###sample hyer circRNAs(100)  for tumor gain####
hypercircRNA.sample100.res.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/DMcircRNA_DEmRNA/hypercircRNA.sample100.res.df",delim = "\t")%>%
  filter(share.pairs>0)%>%mutate(label = "hyper")
dim(hypercircRNA.sample100.res.df)
summary(hypercircRNA.sample100.res.df)
hyper.quant0.95 = quantile(hypercircRNA.sample100.res.df$pairs_gain_normlz,0.25)
hyper.quant0.95

###non DM circRNAs
# sample1000_network_res.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/DEmRNA/sample10k_network_res.df",delim = "\t")%>%
#   filter(share.pairs>0)%>%mutate(label = "non_diff_meth")

nonDMcircRNA.sample802.10k.res.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/DMcircRNA_DEmRNA/nonDMcircRNA.sample802.10k.res.df",delim = "\t")%>%
  filter(share.pairs>0)%>%mutate(label = "non_diff_meth")

summary(sample802_network_res.df)

nonDM.quant0.975 = quantile(nonDMcircRNA.sample802.10k.res.df$pairs_gain_normlz,0.975)
nonDM.quant0.975

##########hyper circRNA with more pairs gain plot : permutation test##########

library(exactRankTests)
library("coin")

###calcu p value for permutation test
permut.p.val = sum(nonDMcircRNA.sample802.10k.res.df$pairs_gain_normlz >= hypercirc.alter[4])/nrow(nonDMcircRNA.sample802.10k.res.df)

#plot#
if (permut.p.val == 0 ){
  pval.text = "p < 1 X 10e4"
}

#normalz
permut.plot <- ggplot(nonDMcircRNA.sample802.10k.res.df,aes(x=pairs_gain_normlz,y=..density..))+
  geom_histogram(color ="white",fill = ggcolors[3],bins = 80)+#scale_fill_manual(values = ggcolors[3:4])+
  stat_density(geom = "line",position = "identity",size=1,adjust = 1)+
  xlab("Co-expression pairs gain in Tumor(normalized)")+theme_classic()+
  geom_vline(xintercept = hypercirc.alter[4],color=ggcolors[1],size=1)+
  geom_vline(xintercept = nonDM.quant0.975,linetype="dashed",size=1)+
  annotate("text",x = c(nonDM.quant0.975+3),y = 0.6,label="alpha=0.05")+
  annotate("text",x = c(hypercirc.alter[4]-4),y = 0.6,label=paste0("Permutation test,\niteration = 10k\n",pval.text))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,(hypercirc.alter[4]+5)))+gg_theme


##############Readers############

reader.DEmRNA.allsamples.cor.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/DEmRNA/DEmRNA_WER.fpkm.allsample.cor.df",delim = "\t")

reader.DEmRNA.allsamples.cor.df = reader.DEmRNA.allsamples.cor.df%>%
  filter(!is.na(cor_value_tumor),wer%in%reader.lst$ensemble_id)%>%
  left_join(dplyr::select(mRNA_res.df,"ensemble_id","symbol"),by=c("wer" = "ensemble_id"))%>%
  mutate(cor_status_tumor = case_when(cor_value_tumor >= cor_cutoff ~ "positive",
                                      cor_value_tumor <= -cor_cutoff ~ "negative",
                                      TRUE ~ "nonSig"),
         cor_status_normal= case_when(cor_value_normal >= cor_cutoff ~ "positive",
                                      cor_value_normal <= -cor_cutoff ~ "negative",
                                      TRUE ~ "nonSig"),
         DEmRNA_circRNAs = ifelse(DEmRNA%in%cor_de_mRNAs,"cor_with_hyperCirc","non_cor_hyper"))
  
dim(reader.DEmRNA.allsamples.cor.df) # 23228
table(reader.DEmRNA.allsamples.cor.df$cor_status_tumor)

##DE mRNA cor with readers
length(unique(filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor!="nonSig")$DEmRNA)) #547

length(unique(filter(reader.DEmRNA.allsamples.cor.df,
                     DEmRNA_circRNAs=="cor_with_hyperCirc",cor_status_tumor!="nonSig")$DEmRNA)) #326

#reader.DEmRNA.df = filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor!="nonSig")
plot.df = filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor!="nonSig")

#order by cor_with_hyperCirc
order.wers.df = as.data.frame(prop.table(xtabs(~symbol+DEmRNA_circRNAs,
                               data = filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor!="nonSig")),margin = 1))%>%
  filter(DEmRNA_circRNAs=="cor_with_hyperCirc")%>%
  dplyr::arrange(desc(Freq))
order.wers.df

# order by count
order.wers.df = as.data.frame(prop.table(xtabs(~symbol,
                                               data = filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor!="nonSig")),margin = 1))%>%
  #filter(DEmRNA_circRNAs=="cor_with_hyperCirc")%>%
  dplyr::arrange(desc(Freq))

order.wers.df = as.data.frame(xtabs(~symbol,data = filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor!="nonSig")))%>%
  dplyr::arrange(desc(Freq))

###plot
#order
order.wers = as.character(order.wers.df$symbol)

plot.df$DEmRNA_circRNAs = factor(plot.df$DEmRNA_circRNAs , levels = c("non_cor_hyper","cor_with_hyperCirc"))
plot.df$symbol = factor(plot.df$symbol,levels = order.wers)

levels(plot.df$symbol)
str(plot.df)

DEmRNA.cor.reader.plot <- ggplot(plot.df,aes(x=symbol))+
  #geom_bar(aes(fill = DEmRNA_circRNAs),stat = "count",position="stack",width = 0.7)+
  stat_count(geom="bar",position = "stack",aes(fill = DEmRNA_circRNAs),width = 0.8) +
  ylab("Count")+xlab("")+
  theme_classic()+scale_fill_manual(values = rev(ggcolors[5:6]))+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12 , color = "black",angle = 60,hjust = 1.1,vjust = 1.1),
        axis.text.y = element_text(size = 12 , color = "black"),
        legend.title= element_blank())
DEmRNA.cor.reader.plot

## 2020-1-8
#width=521&height=357
pdf(file = "/data1/yeying/m6a_circ/final_plot/2020-1-8/figure4/cor_readers.DEmRNA.pdf",width = 5.2,height = 3.57)
DEmRNA.cor.reader.plot
dev.off()

###cor de mRNA###

# reader.cor_de_mRNA.allsamples.cor.df = filter(reader.DEmRNA.allsamples.cor.df,DEmRNA%in%cor_de_mRNAs)
# dim(reader.cor_de_mRNA.allsamples.cor.df) #8759
# 
# plot.df = as.data.frame(xtabs(~symbol+cor_status_tumor,data = reader.cor_de_mRNA.allsamples.cor.df))%>%
#   filter(cor_status_tumor=="positive")%>%arrange(desc(Freq))
# plot.df$symbol = factor(plot.df$symbol,levels = plot.df$symbol)
# 
# ggplot(plot.df,aes(x=symbol,y=Freq))+
#   geom_bar(stat = "identity",width = 0.7,fill = ggcolors[4])+
#   ylab("Count of Correlated DEmRNA")+xlab("")+
#   theme_classic()+
#   scale_y_continuous(expand = c(0, 0))+
#   theme(axis.title = element_text(size = size_axistitle, color = "black"),
#         axis.text.x = element_text(size = 15 , color = "black",angle = 60,hjust = 1.1,vjust = 1.1),
#         axis.text.y = element_text(size = 15 , color = "black"),
#         legend.position = "none")
#strip.text = element_text(size = rel(1)),
#strip.background = element_blank(),strip.placement="outside",
#legend.title= element_blank(),
#legend.background = element_rect(fill = NA),
#legend.justification=c(1,0))

##########8-28 : wer prefer to bind co-Meth DM-DE pairs############
#load data
circRNA_wer_bind.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/WER/circRNA_wer_bind.info.txt",delim = "\t")%>%
  dplyr::select(colnames(circRNA_wer_bind.df)[1:6],reader.lst$symbol)

mRNA_wer_bind.df = read_delim("/data1/yeying/m6a_circ/pancreatic/co-network/WER/mRNA_wer_bind.info.txt",delim = "\t")%>%
  dplyr::select(colnames(mRNA_wer_bind.df)[1:6],reader.lst$symbol)

##add wer bind data
dim(DMcirc_DEcircRNA.allsample.cor.df)
DMcirc_DEcircRNA.cor.allsample.wer_bind.df = left_join(DMcirc_DEcircRNA.allsample.cor.df,
                                                       dplyr::select(circRNA_wer_bind.df,-(chr),-(start),-(end),-(m6Astatus),-(strand)),
                                                       by = c("circRNA_id" = "id"))%>%
  left_join(dplyr::select(mRNA_wer_bind.df,-(chr),-(start),-(end),-(m6Astatus),-(strand)),by = c("mRNA_id" = "id"),suffix = c("_circ","_mRNA"))%>%
  mutate(co_meth_status = case_when(mRNA_m6A=="m6A" & circRNA_m6A=="m6A" ~ "co-Meth",
                                    circRNA_m6A=="m6A" & mRNA_m6A=="non-m6A" ~ "circRNA-meth-only"))%>%dplyr::select(-(log2FoldChange))

hypercirc_DEcircRNA.cor.allsample.wer_bind.df  = DMcirc_DEcircRNA.cor.allsample.wer_bind.df%>%
  filter(meth_status=="hyper")

dim(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)#1226544
dim(hypercirc_DEcircRNA.cor.allsample.wer_bind.df)#938492

length(colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)[grep("_circ",colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df))])
colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df)[grep("_mRNA",colnames(DMcirc_DEcircRNA.cor.allsample.wer_bind.df))]
table(DMcirc_DEcircRNA.cor.allsample.wer_bind.df$co_meth_status)
# write_delim(DMcirc_DEcircRNA.cor.allsample.wer_bind.df,
#             path = "/data1/yeying/m6a_circ/pancreatic/co-network/WER/DMcirc_DEcircRNA.cor.allsample.wer_bind.df",delim = "\t")

##add co bind status

hypercirc_DEcircRNA.co_bind.df = dplyr::select(hypercirc_DEcircRNA.cor.allsample.wer_bind.df,-contains("_circ"),-contains("_mRNA")) 

colnames(hypercirc_DEcircRNA.cor.allsample.wer_bind.df)
colnames(DMcirc_DEcircRNA.co_bind.df)
reader.lst$symbol

#loop
for (i in 1:10){
  #i=1
  tmp.wer = gsub("_circ","",colnames(hypercirc_DEcircRNA.cor.allsample.wer_bind.df)[(i+10)])
  tmp.wer.bind.df = hypercirc_DEcircRNA.cor.allsample.wer_bind.df[,c(1,4,(i+10),(i+10+10))]
  colnames(tmp.wer.bind.df) = c("mRNA_id","circRNA_id","circ_bind","mRNA_bind")
  tmp.wer.bind.df =  mutate(tmp.wer.bind.df,co_bind_status = case_when( 
    circ_bind > 0 & mRNA_bind > 0 ~ "co-bind",
    circ_bind > 0 & mRNA_bind == 0 ~ "circRNA-bind-only",
    circ_bind == 0 & mRNA_bind > 0 ~ "mRNA-bind-only",
    circ_bind == 0 & mRNA_bind == 0 ~ "co-non-bind"))%>% dplyr::select(-(circ_bind),-(mRNA_bind))
  tmp.wer.bind.df$co_bind_status = factor(tmp.wer.bind.df$co_bind_status,levels = c("co-bind","circRNA-bind-only","mRNA-bind-only","co-non-bind"))
  
  hypercirc_DEcircRNA.co_bind.df = inner_join(hypercirc_DEcircRNA.co_bind.df,tmp.wer.bind.df)
  colnames(hypercirc_DEcircRNA.co_bind.df)[i+11] = tmp.wer
  print(tmp.wer)
}

colnames(hypercirc_DEcircRNA.co_bind.df)
#write_delim(hypercirc_DEcircRNA.co_bind.df,path = "hypercirc_DEcircRNA.co_bind.df",delim = "\t")
#write_delim(DMcirc_DEcircRNA.co_bind.df,path = "DMcirc_DEcircRNA.wer_bindstatus.df",delim = "\t")

##### plot df #####
#hyper_cor_pairs_wer_bind.df = filter(DMcirc_DEcircRNA.cor.allsample.wer_bind.df,meth_status=="hyper")
###dcast bind.df
co_bind.plot.df = melt(DMcirc_DEcircRNA.co_bind.df,id.vars = colnames(DMcirc_DEcircRNA.co_bind.df)[c(1:11)])%>%
  mutate(bind_status = ifelse(value=="co-bind","co-bind pairs","others"),
         pairs_meth = ifelse(co_meth_status=="co-Meth","co-Meth pairs","others"))

co_bind.plot.df$value = factor(co_bind.plot.df$value,levels = c("co-bind","circRNA-bind-only","mRNA-bind-only","co-non-bind"))
final.wers = as.character(unique(co_bind.plot.df$variable))
length(final.wers)
co_bind.plot.df$variable = factor(co_bind.plot.df$variable , 
                                  levels = c(sort(final.wers[c(grep("YTH",final.wers),grep("IGF",final.wers),grep("HNR",final.wers))])))
unique(co_bind.plot.df$variable)
head(co_bind.plot.df$co_meth_status)
colnames(co_bind.plot.df)

#hyper and hypo circRNAs co pairs 
colnames(co_bind.plot.df)[which(colnames(co_bind.plot.df)=="meth_status")] = "circRNA_DM"
hyper_cor_pairs_wer_bind.df = filter(co_bind.plot.df,circRNA_DM=="hyper")
#hypo_cor_pairs_wer_bind.df = filter(co_bind.plot.df,circRNA_DM=="hypo")
head(hyper_cor_pairs_wer_bind.df)

#########total_wer.plot####
hyper_wer.plot.df = rbind(dplyr::select(hyper_cor_pairs_wer_bind.df,variable,bind_status,pairs_meth)%>%mutate(meth_status = "total_pairs"),
                          dplyr::select(filter(hyper_cor_pairs_wer_bind.df,co_meth_status=="co-Meth"),variable,bind_status,pairs_meth)%>%
                            mutate(meth_status = "co-meth_pairs"))
head(hyper_wer.plot.df$meth_status)
hypo_wer.plot.df = rbind(dplyr::select(hypo_cor_pairs_wer_bind.df,variable,bind_status)%>%mutate(meth_status = "total_pairs"),
                         dplyr::select(filter(hypo_cor_pairs_wer_bind.df,co_meth_status=="co-Meth"),variable,bind_status)%>%
                           mutate(meth_status = "co-meth_pairs"))
head(hypo_wer.plot.df$meth_status)
# save.image()

##########calcu fisher test p-value for co-bind####
###function
calcu_fisher = function(wer.df){
  #wer.df=hypo_cor_pairs_wer_bind.df
  wer_co_bind.fisher.df = data.frame()
  for (i in 1:length(levels(wer.df$variable))){
    #i=2
    wer.name = levels(wer.df$variable)[i]
    wer.plot.stats.tb = xtabs( ~ pairs_meth+bind_status,
                               data = filter(wer.df,variable==wer.name))
    #print(wer.plot.stats.tb)
    tmp.fisher = fisher.test(wer.plot.stats.tb,alternative = "greater")
    tmp.vec = c(id=wer.name,odd = tmp.fisher$estimate ,pvalue = as.numeric(tmp.fisher$p.value))
    print(tmp.vec)
    wer_co_bind.fisher.df = rbind(wer_co_bind.fisher.df,tmp.vec)
  }
  
  colnames(wer_co_bind.fisher.df) = c("wer","odds","pvalue")
  wer_co_bind.fisher.df = mutate(wer_co_bind.fisher.df,group1="total_pairs",
                                 group2 = "co-meth_pairs", p_value  =signif(as.numeric(pvalue),digits = 3))
  return(wer_co_bind.fisher.df)
}

###fisher test
#levels(hyper_cor_pairs_wer_bind.df$variable)
hyper_fisher.df = calcu_fisher(hyper_cor_pairs_wer_bind.df)
str(hyper_fisher.df)

#hypo_fisher.df = calcu_fisher(hypo_cor_pairs_wer_bind.df)
###add pvalue to plot
min(hyper_fisher.df$p_value[hyper_fisher.df$p_value>0])
hyper_fisher.df = mutate(hyper_fisher.df,pvalue_revised = ifelse(p_value==0,min(hyper_fisher.df$p_value[hyper_fisher.df$p_value>0])/1000,p_value),log10p=-log10(pvalue_revised))


reader.co_bind.pvalue.plot <- ggbarplot(data = hyper_fisher.df,x = "wer",y = "log10p",color = "white",
                                        fill = ggcolors[3],sort.val ="asc", orientation = c("horizontal") )+
  ylab("-log10(pvalue)")+scale_y_continuous(expand = c(0, 0))+gg_theme

reader.co_bind.pvalue.plot
