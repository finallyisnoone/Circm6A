library(survival)
library(dplyr)
library(survminer)
library(ggplot2)
library(readr)
library(cowplot)
#setwd("/data/xingyang/m6A_zhengjian")
options(stringsAsFactors=F)
setwd("/data1/yeying/m6a_circ/pancreatic/survival/3-24-2020/")
#
#DM_circ_mat0 = read_delim("/data1/yeying/m6a_circ/analysis/survival/DM_DE_input_correct_RPM_m6Astatus_9_14.matrix",delim = "\t")
#str(DM_circ_mat)
########warning: has been changed to DM

#meth.mat$id = rownames(meth.mat)

DM_circ_mat0 = meth.mat[DM_list,]
dim(DM_circ_mat0)

DM_circ_mat=as.data.frame(t(as.matrix(DM_circ_mat0)))
DM_circ_mat[,1:10]

DM_circ_mat=cbind(t_id=gsub(pattern = "S_",replacement = "",as.character(rownames(DM_circ_mat))),DM_circ_mat)
rownames(DM_circ_mat)
str(DM_circ_mat)
#####load clinic data
#clinc_data=read.csv("survival/m6A_clinc_analysis_2.csv")
#clinc_data=read_delim("cdclinical_data-3-24.txt",delim = "\t")
clinc_data=read_delim("/data1/yeying/m6a_circ/pancreatic/survival/3-24-2020/tmp.txt",delim = "\t")
os_data0=inner_join(clinc_data,as.data.frame(DM_circ_mat),by="t_id")
dim(os_data0);os_data0[,c(34:40)]
os_data = as.data.frame(os_data0[which(!is.na(os_data0$sample_id_num)),])
os_data$t_id=as.character(os_data$t_id)

str(os_data)


# i=888
# k=5
PFS.res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
head(PFS.res)
PFS.res$FDR = p.adjust(PFS.res$p.value,method = "fdr")
PFS.res = data.frame(id = rownames(PFS.res),PFS.res)%>%arrange(FDR)
write.xlsx(PFS.res,file = "PFS_cox_model.result.xlsx")

#########significance
# sig_survival_res=as.data.frame(t(apply(survival_res,1,function(x){p_min <-min(as.numeric(x[3:6]))
# return(c(x,p_min=p_min))
# })))%>%filter(p_min <= 0.05 )%>%arrange(p_min,zero_num)

survival_res$p_logtest =as.numeric(survival_res$p_logtest)
survival_res$p_waldtest =as.numeric(survival_res$p_waldtest)
survival_res$p_logranktest = as.numeric(survival_res$p_logranktest)
survival_res$padj_logrank = p.adjust(survival_res$p_logranktest)
sig_survival_res = filter(survival_res,as.numeric(p_logranktest) <= 0.05)
unique(sig_survival_res$gene)

names(sig_survival_res)[1] <- "id"

sig_survival_DM_feature = merge(sig_survival_res,circ_feature,by = "id")

#write_delim(sig_survival_DM_feature,path = "survival/sig_survival_DM_feature_11_7.txt",delim = "\t",col_names = T)

#write.table(sig_survival_res,"/data1/yeying/m6a_circ/pancreatic/survival/sig_survival_res.txt",col.names = T,row.names = F,quote = F,sep="\t")

##########plot

plot_survial <- function(data_df,title_sur){
  library(survminer)
  
  p = ggsurvplot(survfit(Surv(time,event)~label,data = data_df),data = data_df,title=title_sur,
                 palette = c("#00468B", "#F24747"),legend.labs=c("low","high"),pval = TRUE,pval.method = TRUE, 
                 risk.table=T)
  p =plot_grid(p$plot+scale_x_continuous(expand = c(0,0)),p$table,ncol = 1,rel_heights = c(3, 1))
  
  return(p)
}

#################plot all sig circRNA###############
#pdf(file = "/data1/yeying/m6a_circ/final_plot/DM_survival_res_12_6.pdf")
pdf(file = "/data1/yeying/m6a_circ/pancreatic/RBP/11_7/2_21/COL1A1_survival_res.pdf")

for (c in 1:dim(sig_survival_res)[1]){
  #####i is col_num in sur
  # i=1162
  # k=13
  # data_df=os_data
  gene_name=paste(sig_survival_DM_feature$symbol[c],
                  sig_survival_DM_feature$id[c],sep = ":")
  title_sur=paste("survival analysis of",gene_name)
  col_num=as.numeric(sig_survival_res$col_number[c])
  cutoff=as.numeric(sig_survival_res$n_low[c])
  cat(gene_name,"P value = ",sig_survival_DM_feature$p_logranktest[c],"\n",
      "HR =",sig_survival_DM_feature$HR[c],"\n")
  
  temp=data.frame(t_id=os_data$t_id,RPM=as.numeric(os_data[,col_num]))
  temp=temp%>%arrange(RPM)%>% mutate(label= c(rep(1,cutoff),rep(2,nrow(temp)-cutoff)))
  sur_data=data.frame(t_id=os_data$t_id,time = os_data$month,event = os_data$status)
  sur_data = inner_join(sur_data,temp,by="t_id")
  ###df : time ,event , label(high or low)
  p = plot_survial(data_df = sur_data,title_sur=title_sur)
  
  print(p)
  
}

dev.off()

####################other risk factor
str(os_data)

###########multi risk factor vocano plot##############
#colnames(result)=c("Peak.id","Clinic.factor","P.value","log2FC")
clinc_data_circ <- os_data[,1:34]
clinc_data_circ = mutate(clinc_data_circ,Differentiation_new = case_when(Differentiation==1 ~ 1,
                                                                         Differentiation==2 ~ 1,
                                                                         Differentiation==3 ~ 2))
rownames(clinc_data_circ) = paste0("S_",clinc_data_circ$t_id)


peak2gene= filter(circ_feature, id%in%DM_list)
peak2gene = data.frame(Peak.id = peak2gene$id,gene = peak2gene$symbol)
peak2gene_nonNFPB = filter(peak2gene,!grepl("NBPF",gene)) 
peak2gene_nonNFPB = filter(peak2gene_nonNFPB,!grepl("RP11-458D21.5",gene))
tmp = as.data.frame(table(peak2gene_nonNFPB$gene))
filter(tmp,Freq > 10)
dim(peak2gene_nonNFPB)
tumor_meth.mat = meth.mat[peak2gene_nonNFPB$Peak.id,rownames(clinc_data_circ)]
dim(tumor_meth.mat)
clinc_data_circ = clinc_data_circ[colnames(tumor_meth.mat),]
####PNI perineural invasion

risk_factors = c("gender","smoke","drink","PNI","Angiogenesis",
                 "node","TNM_3","Relapse","status") #"Metastasis.Relapse"

clinc_data_circ$TNM_3[which(is.na(clinc_data_circ$TNM_3))] <- 1
all_risk_factor = data.frame()
attach(clinc_data_circ)
i=1
for ( i in 1:length(risk_factors) ) {
  pheno.df <- data.frame(Type= as.factor(get(risk_factors[i])))
  rownames(pheno.df) = rownames(clinc_data_circ)
  MethCirc_tmp<- mat_wilcox(expre_mat = tumor_meth.mat,test_mode = "",
                            group_mat = pheno.df)
  #MethCirc_edgeR <-edgeR_test(expre_mat = tumor_meth.mat,group_mat = pheno.df)
  factor_res.df = data.frame(id = rownames(MethCirc_tmp),
                             Clinic.factor =rep(risk_factors[i],nrow(MethCirc_tmp)),
                             P.value = MethCirc_tmp$PValue ,FDR=MethCirc_tmp$FDR, log2FC = MethCirc_tmp$logFC)
  all_risk_factor <- rbind(all_risk_factor,factor_res.df)
} 

detach(clinc_data_circ)
head(all_risk_factor)

write.xlsx(all_risk_factor,file = "other_risk_factor.res.xlsx")

clinic_heatmap_for_common_hyper=function(peak_associate_with_clinic,top=10,peak2gene){
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$P.value<0.05),]
  peak_associate_with_clinic=merge(peak_associate_with_clinic,peak2gene,by="Peak.id",all.x = T)
  peak_associate_with_clinic=peak_associate_with_clinic[which(!is.na(peak_associate_with_clinic$gene)),]
  peak_associate_with_clinic=as.data.frame(cbind(peak_associate_with_clinic,type="non"))
  peak_associate_with_clinic=peak_associate_with_clinic[order(peak_associate_with_clinic$P.value),]
  i=1
  Clinic.factor=unique(peak_associate_with_clinic$Clinic.factor)
  n=length(Clinic.factor)
  while(i<=n){
    if(length(peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'Peak.id'])<top){
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i]),'type']='sig'
    } else {
      peak_associate_with_clinic[which(peak_associate_with_clinic$Clinic.factor==Clinic.factor[i])[1:top],'type']='sig'
    }
    i=i+1
  }
  peak_associate_with_clinic=peak_associate_with_clinic[which(peak_associate_with_clinic$type=="sig"),]
  peak_associate_with_clinic$P.value=as.numeric(peak_associate_with_clinic$P.value)
  peak_associate_with_clinic$Peak.id=paste(peak_associate_with_clinic$Peak.id,":",peak_associate_with_clinic$gene,sep="")
  # peak_associate_with_clinic$Peak.id=ordered(peak_associate_with_clinic$Peak.id,
  #                                            levels=names(table(peak_associate_with_clinic$Peak.id)[
  #                                              order(table(peak_associate_with_clinic$Peak.id),decreasing = T)]))
  peak_associate_with_clinic=data.frame(cbind(peak_associate_with_clinic,Difftype="Up regulation"))
  peak_associate_with_clinic[which(peak_associate_with_clinic$log2FC<0),"Difftype"]="Down regulation"
  #peak_associate_with_clinic$Clinic.factor=ordered(peak_associate_with_clinic$Clinic.factor,
  #                                                 levels=rev(c("age","gender","drink","smoke","TNM_3","xueguan","node","nerve",
  #                                                              "fenhua_class : 1  vs Rest","fenhua_class : 2  vs Rest",
  #                                                              "fenhua_class : 3  vs Rest")))
  write.table(sort(peak_associate_with_clinic$Peak.id),"survival/peak_ass.temp.txt",row.names = F,
              quote = F)
  p=ggplot(peak_associate_with_clinic)+
    geom_point(aes(y=peak_associate_with_clinic$Clinic.factor,
                   x=peak_associate_with_clinic$Peak.id,size=(-log10(peak_associate_with_clinic$P.value)),
                   color=as.factor(peak_associate_with_clinic$Difftype)))+
    scale_color_manual(values=c("#5707d5","#d50709"),name="log2(Flod change)")+
    scale_size_continuous(name="-log10(P value)")+
    xlab("Clinic factor")+ylab("Tumor specific peaks top10 in each clinic factor")+
    theme_classic()+theme(axis.text.x = element_text(angle=90))
  return(p)
}

#ggsave(risk_factor_DM_dotplot,filename = "../final_plot//all_risk_factor_DM.pdf",width = 20,height = 8)

risk_factor_DM_dotplot_pvalue <- clinic_heatmap_for_common_hyper(peak_associate_with_clinic = all_risk_factor,
                                                                 top = 10,peak2gene = peak2gene) ##### peak2gene : peakID genesymbol

all_risk_factor_FDR = all_risk_factor 
all_risk_factor_FDR$P.value = all_risk_factor_FDR$FDR
head(all_risk_factor_FDR);head(peak2gene)
filter(all_risk_factor_FDR,FDR <= 0.05)%>%left_join(peak2gene,by="Peak.id")


risk_factor_DM_dotplot_FDR <- clinic_heatmap_for_common_hyper(peak_associate_with_clinic = all_risk_factor_FDR,
                                                              top = 10,peak2gene = peak2gene) ##### peak2gene : peakID genesymbol




#####3-27#####
###all circ
circ2symbol = dplyr::select(circ_feature,id,symbol)
#
colnames(os_data)[1:34]
os_data$gender = as.factor(os_data$gender)
os_data$age
os_data$smoke = as.factor(os_data$smoke)
os_data$drink = as.factor(os_data$drink)
DM.circ.lst = colnames(os_data)[35:1269] ###all the circRNA
#DM.circ.lst = hyper.circ
#+PNI+Angiogenesis+node+TNM_3+zhuanyi+Relapse
# univ_formulas <- sapply(covariates,
#                         function(x) as.formula(paste0('Surv(PFS_month)~ `', x,
#                                                      "`+gender+age+smoke+drink")))
cox_model <- function(clinical.data,covariates,responer,circ2symbol){
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste0('Surv(',responer,')~ `', x,
                                                        "`+gender+age+smoke+drink")))
  
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = clinical.data)})
  univ_models
  
  #result extraction
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           ######
                           p.value<-signif(x$coefficients[1,5], digits=2)
                           #wald.test<-signif(x$wald["test"], digits=2)
                           #beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coefficients[1,2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[1,"lower .95"],2)
                           HR.confint.upper <- signif(x$conf.int[1,"upper .95"],2)
                           HR.ci <- paste0(HR, " (", 
                                           HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(HR,HR.ci, p.value)
                           names(res)<-c("HR","CI_95", "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  class(univ_results)
  
  res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  res$FDR = p.adjust(res$p.value,method = "fdr")
  res = data.frame(id = rownames(res),res)%>%arrange(FDR)%>%left_join(circ2symbol)
  return(res)
}

PFS.cox.df = cox_model(clinical.data = os_data,covariates = DM.circ.lst,responer = "PFS_month",circ2symbol = circ2symbol)
head(PFS.cox.df);dim(PFS.cox.df)
#write.xlsx(PFS.res,file = "PFS_cox_model.result.xlsx")

#OS

os_data$status = as.factor(os_data$status)
os.cox.df = cox_model(clinical.data = os_data,covariates = DM.circ.lst,responer = "month",circ2symbol)

# write.xlsx(os.cox.df,file = "OS_cox_model.result.xlsx")

####plot####
fdr.cutoff = 0.15

#vocano plot
volcano.cox.model <- function(cox.df,plot.title = "",fdr.cutoff= 0.15,x_max = ""){
  #
  tmp = cox.df%>%
    mutate(log2HR = log2(as.numeric(HR)),
           label = case_when(log2HR > 0 & FDR <= fdr.cutoff ~ "up",
                             log2HR < 0 & FDR <= fdr.cutoff  ~ "down",
                             TRUE ~ "nonSig"),
           logFDR = -log10(FDR))
  tmp$label = factor(tmp$label , levels = c("up","nonSig","down"))
  if (x_max == ""){
    x_max = max(abs(tmp$log2HR),na.rm = T) 
  }
  summary(tmp$log2HR)
  sig.stats = table(tmp$label)
  #plot
  p = ggplot(data = tmp,aes(x= log2HR ,y=logFDR,color=label))+geom_point(size = 2,alpha = 0.8)+
    geom_vline(xintercept = 0,linetype="dashed")+
    geom_hline(yintercept = -log10(fdr.cutoff),linetype="dashed")+
    xlim(-x_max,x_max)+ylab("-log10(FDR)")+ggtitle(plot.title)+
    scale_color_manual(values = c(ggcolors[1],"grey50",ggcolors[2]))+
    annotate("text",label = paste("Up = ",sig.stats[1]),color = ggcolors[1],x = x_max/2 , y= max(tmp$logFDR,na.rm = T))+ 
    annotate("text",label = paste("Down =",sig.stats[3]),color = ggcolors[2],x = -x_max/2,y= max(tmp$logFDR,na.rm = T))+
    annotate("text",label = paste0("FDR = ",fdr.cutoff),x = -x_max*0.8,y=-log10(p_value)*0.5)+
    theme_test()+
    gg_theme+theme(legend.position = "top",legend.spacing.x = unit(0.2,"cm"))
  
  return(p)
}

summary(as.numeric(PFS.cox.df$HR))
PFS.volcano.plot = volcano.cox.model(cox.df = PFS.cox.df,plot.title = "PFS",x_max = 2) 
PFS.volcano.plot
## x.max = 1.5 -> Removed 7 rows containing missing values
#width=354&height=370
pdf(file = "PFS.cox.volcano.pdf",width = 3.54,height = 3.7)
PFS.volcano.plot
dev.off()

#OS plot
summary(as.numeric(os.cox.df$HR))
OS.volcano.plot = volcano.cox.model(cox.df = os.cox.df,plot.title = "OS",x_max = 2) 
OS.volcano.plot
## x.max = 2 -> Removed 1 rows containing missing values
#width=334&height=366
pdf(file = "OS.cox.volcano.pdf",width = 3.34,height = 3.66)
OS.volcano.plot
dev.off()

##
#enrichr_barplot(res.df = PFS.cox.df,top_number = 10,bar_color = )
PFS.cox.df = PFS.cox.df%>%mutate(log10FDR=-log10(FDR),
                                 label = paste(symbol,id,sep=":"))

#width=585&height=374
pdf("PFS.cox.barplot.pdf",width = 5.85,height = 3.74)
ggbarplot(data = filter(PFS.cox.df,FDR <= fdr.cutoff),x = "label",y = "log10FDR",
          color = "white",fill = "#f00a36",sort.val = "asc",orientation = "horiz")+
  geom_hline(yintercept = -log10(fdr.cutoff),linetype = "dashed")+
  ggtitle("PFS")+
  xlab("")+ylab("-log10(FDR)")+scale_y_continuous(expand = c(0,0))+theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    legend.spacing.x = unit(0.2,"cm"),
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA))

dev.off()

###########os barplot#####
os.cox.df = os.cox.df%>%mutate(log10FDR=-log10(FDR),
                                 label = paste(symbol,id,sep=":"))

#width=595&height=340
pdf("OS.cox.barplot.pdf",width = 6,height = 3.4)
ggbarplot(data = filter(os.cox.df,FDR <= fdr.cutoff),x = "label",y = "log10FDR",
          color = "white",fill = "#007cc0",sort.val = "asc",orientation = "horiz")+
  geom_hline(yintercept = -log10(fdr.cutoff),linetype = "dashed")+
  ggtitle("OS")+
  xlab("")+ylab("-log10(FDR)")+scale_y_continuous(expand = c(0,0))+theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    legend.spacing.x = unit(0.2,"cm"),
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA))
dev.off()
