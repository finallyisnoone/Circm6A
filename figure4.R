#########
#m6A effect on circRNA
##########plot init############
##font size
size_axisblacktext <- 18 #24
size_axisgreytext <- 12 #18
size_legendtext <- 10  #22
size_axistitle <- 15

##colors for ggplot 
ggcolors <- c("#BC3C29","#0072B5","#E18727","#20854E","#80397B")
show_col(ggcolors)

## basic theme config for all ggplot
gg_theme <- theme(
  axis.title = element_text(size = size_axistitle, color = "black"),
  axis.text.x = element_text(size = size_axistitle, color = "black"),
  axis.text.y = element_text(size = size_axistitle, color = "black"),
  legend.key.size = unit(0.3,"cm"),
  legend.title= element_blank(),
  legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
  legend.background = element_rect(fill = NA))

#usual for color config : +scale_y_continuous(expand = c(0, 0))+scale_color_manual(values = ggcolors[1:2])+gg_theme


################m6A & expression correlation ##########
###hist of each m6A circRNA
# circMeth_RPM <- ggplot(circ_meth_RPM_cor_df2,aes(x=cor_value, y=..density..))+theme_classic()+
#   geom_histogram(fill = ggcolors[1],color = "white",breaks = seq(-0.2, 1.0, 0.05))+
#   geom_density(size=1)+xlab("Spearman Rank Correlation")+
#   scale_y_continuous(expand = c(0, 0))+#scale_fill_manual(values = ggcolors[2])+
#   scale_x_continuous(expand = c(0, 0),limits  = c(-0.2,1.0),breaks = seq(-0.2, 1.0, 0.2))+gg_theme

circMeth_RPM <- ggplot(circ_meth_RPM_cor_df2,aes(x=cor_value))+theme_classic()+
  geom_histogram(fill = ggcolors[5],color = "white",breaks = seq(-0.2, 1.0, 0.1))+
  xlab("Spearman Rank Correlation (each circRNA) ")+
  scale_y_continuous(expand = c(0, 0))+#scale_fill_manual(values = ggcolors[2])+
  scale_x_continuous(expand = c(0.01, 0.01),limits  = c(-0.2,1.0),breaks = seq(-0.2, 1.0, 0.1))+gg_theme

circMeth_RPM

summary(circ_meth_RPM_cor_df.samples)
circMeth_RPM.samples <- ggplot(circ_meth_RPM_cor_df.samples,aes(x=cor_value))+theme_classic()+
  geom_histogram(fill = ggcolors[4],color = "white",breaks = seq(0.1, 0.9, 0.05))+
  xlab("Spearman Rank Correlation (each sample)")+
  scale_y_continuous(expand = c(0, 0))+#scale_fill_manual(values = ggcolors[2])+
  scale_x_continuous(expand = c(0.01, 0.01),limits  = c(0.1,0.9),breaks = seq(0.1, 0.8, 0.1))+gg_theme



pdf("/data1/yeying/m6a_circ/final_plot/figure4/meth_expres.pdf",height = 5.3,width = 12.5)
plot_grid(circMeth_RPM,circMeth_RPM.samples)
dev.off()
###point plot
m6A_expre_cor.plot

#################m6A & circularization ratio##########

ratio_ks_res <- ks.test(filter(circ_ratio_m6A,m6Astatus=="m6A")$mean,
                        filter(circ_ratio_m6A,m6Astatus=="non-m6A")$mean)
ratio_wilx_res <- wilcox.test(filter(circ_ratio_m6A,m6Astatus=="m6A")$mean,
                              filter(circ_ratio_m6A,m6Astatus=="non-m6A")$mean,alternative = "greater")

circRatio_ecdf <- ECDF_plot(df = circ_ratio_m6A,value_var = "mean",group_var = "m6Astatus",test_result = ratio_wilx_res)+
  theme(legend.position=c(0.9,.1))+
  xlab("circular ratio")+scale_color_manual(values = ggcolors[1:2])+gg_theme

circRatio_ecdf

###############total_cov density##################################
total_cov.density.plot = ggplot(total_circ.cov.stats.df,
                                aes_string(x="bin_final",y="total_cov_density",linetype="m6Astatus",color="tag"))+facet_grid(tag ~ .  ,scales = "free_y")+
  geom_line(size = 1.5)+scale_color_manual(values = ggcolors )+
  xlab("")+ylab(label ="binding sites density")+scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
        axis.text.x = element_text(size = size_axisblacktext, color = "black",angle = 60,hjust = 1.1,vjust = 1.1),
        axis.text.y = element_text(size = size_axisblacktext, color = "black"),
        strip.text = element_text(size = rel(1)),
        strip.background = element_blank(),strip.placement="outside",
        legend.title= element_blank(),
        legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
        legend.background = element_rect(fill = NA),
        legend.justification=c(1,0),
        legend.key.height = unit(1.5,"line"))

pdf(paste0(outdir,"/total_cov.density.pdf"),width = 7,height = 10)
total_cov.density.plot
dev.off()


###################### m6A circRNAs with TIS#############################

circm6a_TIS_stats.df = as.data.frame(prop.table(xtabs(~m6Astatus+TIS,circ_feature.TIS),margin = 1))
circm6a_TIS.plot = ggplot(circm6a_TIS_stats.df,aes(x=m6Astatus,fill=TIS,y=Freq))+
  geom_bar(stat = "identity",width = 0.5)+geom_text(aes(label=rev(signif(Freq,3))), color = "white", vjust = 0.6,size = 5)+
  xlab("")+ylab("Percentage")+ggtitle("Fisher's Exact Test , p-value = 0.038")+
  theme_classic()+scale_fill_manual(values=ggcolors[3:4])+ theme_classic ()+ scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
        axis.text.x = element_text(size = size_axisblacktext, color = "black"),
        axis.text.y = element_text(size = size_axisblacktext, color = "black"),
        legend.title= element_blank(),
        legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
        legend.background = element_rect(fill = NA),
        legend.justification=c(1,0),
        legend.key.height = unit(1.5,"line"))#

circm6a_TIS.plot

######################m6A and eIFs###############################
#####m6A circRNA eIFs and YTH##########
eIFs_total.df = rbind(mRNA_eif.df,circRNA_eif.df,circRNA.m6A.TIS.df)%>%filter(RBP!="eIF4A3")
eIFs_total.df$RBP = factor(eIFs_total.df$RBP,levels = c("eIF3","eIF4G2","eIF4G1",paste0("YTHDF",1:3)))
eIFs_total.df$type = factor(eIFs_total.df$type,levels = c( "m6AcircRNA_withTIS","circRNA","mRNA"))

eIFs.plot <- ggplot(filter(eIFs_total.df,RBP!="eIF4G1"),aes(x=RBP,fill=type,y=fraction))+
  geom_bar(stat = "identity",position = "dodge",width = 0.5)+
  theme_classic()+scale_fill_manual(values=ggcolors)+ theme_classic ()+ scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
        axis.text.x = element_text(size = size_axisblacktext, color = "black"),
        axis.text.y = element_text(size = size_axisblacktext, color = "black"),
        legend.title= element_blank(),
        legend.position = "top",
        legend.text = element_text(size = size_legendtext), 
        legend.background = element_blank(),
        legend.key.height = unit(1.5,"line"))+xlab("")+ylab("Percentage")

eIFs.plot

#####eIFs choose##########
eIFs_total.df2 = rbind(mRNA_eif.df,circRNA_eif.df,circRNA.m6A.TIS.df)%>%filter(RBP%in%c("eIF4G1","eIF4G2"))
eIFs_total.df2$RBP = factor(eIFs_total.df2$RBP,levels = c("eIF4G2","eIF4G1"))
eIFs_total.df2$type = factor(eIFs_total.df2$type,levels = c( "m6AcircRNA_withTIS","circRNA","mRNA"))

eIFs_circRNA_mRNA.plot = ggplot(eIFs_total.df2,aes(x=RBP,fill=type,y=fraction))+
  geom_bar(stat = "identity",position = "dodge",width = 0.5)+
  theme_classic()+scale_fill_manual(values=ggcolors)+ scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
        axis.text.x = element_text(size = size_axisblacktext, color = "black"),
        axis.text.y = element_text(size = size_axisblacktext, color = "black"),
        legend.title= element_blank(),
        legend.position = "top",
        legend.text = element_text(size = size_legendtext), 
        legend.background = element_blank(),
        legend.key.height = unit(1.5,"line"))+xlab("")+ylab("Percentage")
eIFs_circRNA_mRNA.plot
pdf(paste0(outdir,"/eIFs_circRNA_mRNA.plot.pdf"),width = 4.5,height = 4)
eIFs_circRNA_mRNA.plot
dev.off()

##############merge plot################

first_row = plot_grid(circMeth_RPM,circMeth_RPM.samples,labels = c("A",""),nrow = 1)
first_row
v_col = plot_grid(circRatio_ecdf,circm6a_TIS.plot,labels = c("B","D"),ncol = 1)
median_row = plot_grid(v_col,total_cov.density.plot,labels = c("","C"),ncol = 2,rel_widths = c(1, 1.5))
#output
outdir="/data1/yeying/m6a_circ/final_plot/figure4/"
pdf(paste(outdir,"figure4.pdf",sep = "/"),height = 15,width = 9.7)
plot_grid(first_row,median_row,eIFs.plot,ncol = 1,rel_heights = c(1.2,2,1.2))
dev.off()



###############figure 3 part#################### 
################miRNA##############
allcirc_miRNAtarget.info <-  read_delim("/data1/yeying/m6a_circ/analysis/miRNA/11_7/all_filter_circ_m6Astatus_targetInfo.txt",delim = "\t")
allcirc_miRNAtarget.info= filter(allcirc_miRNAtarget.info,chr!="chrM")%>%mutate(id = paste(chr,start,end,strand,sep = "_"))
dim(allcirc_miRNAtarget.info)
m6ACirc.count = length(which(allcirc_miRNAtarget.info$m6Astatus=="m6A"))
nonm6ACirc.count = length(which(allcirc_miRNAtarget.info$m6Astatus=="non-m6A"))
  
######fisher test
allcirc_stats.tb = fisher.test(xtabs(~m6Astatus+miRNAtarget,data = allcirc_miRNAtarget.info))
allcirc_stat2 = filter(as.data.frame(xtabs(~m6Astatus+miRNAtarget,data = allcirc_miRNAtarget.info)),miRNAtarget=="miRNA_target")
allcirc_stat2$Percent = c(allcirc_stat2$Freq[1]/m6ACirc.count,
                          allcirc_stat2$Freq[2]/nonm6ACirc.count)
allcirc_stats.tb
allcirc_miRNA_barplot <- ggplot(allcirc_stat2,aes(x=m6Astatus,y=Percent,fill=m6Astatus))+geom_bar(stat = "identity",width = 0.5)+
  scale_fill_manual(values = ggcolors[1:2])+xlab("")+scale_y_continuous(expand = c(0,0))+gg_theme
allcirc_miRNA_barplot

allcirc_stats.tb = as.data.frame(prop.table(xtabs(~m6Astatus+miRNAtarget,data = allcirc_miRNAtarget.info),margin = 1))
miRNA.target.plot <- ggplot(allcirc_stats.tb,aes(x=m6Astatus,y=Freq,fill=miRNAtarget))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+
  theme_classic()+scale_fill_manual(values=ggcolors)+ theme_classic ()+ scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
    axis.text.x = element_text(size = size_axisgreytext, color = "black"),
    axis.text.y = element_text(size = size_axisblacktext, color = "black"),
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA),
    legend.justification=c(1,0),
    #legend.position=c(1, 0),
    legend.key.height = unit(1.5,"line")
  )+xlab("")+ylab("Percentage")


pdf(paste0(outdir,"/miRNA_target_percent.pdf"))
miRNA.target.plot
dev.off()

