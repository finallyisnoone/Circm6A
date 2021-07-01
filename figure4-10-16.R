
####plot for fig4
#####co-network in tumor#####

m6a_nonm6a_cor.plot <- 
  ggplot(sample100k.totalCirc_DEmRNA.allsample.cor.df)+
  stat_density(geom = "line",position = "identity",aes(x=cor_value_tumor, y=..density..,color=m6Astatus),size=1,adjust = 2.5)+xlab("Pearson Correlation")+
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
  annotate("text",x= 0.75,y=1.8,label=paste0("N = ",nonm6ACirc_DEmRNA.tumor.stat[1]+nonm6ACirc_DEmRNA.tumor.stat[3]),color = ggcolors[2])+
  theme(
    axis.title = element_text(size = size_axistitle, color = "black"),
    axis.text.x = element_text(size = size_axistitle, color = "black"),
    axis.text.y = element_text(size = size_axistitle, color = "black"),
    legend.position = c(0.35,0.85),
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA))

m6a_nonm6a_cor.plot
# venn for tumor pairs (N of m6A vs N of nonm6A)



#######normalz#####
permut.plot <- ggplot(nonDMcircRNA.sample802.10k.res.df,aes(x=pairs_gain_normlz,y=..density..))+
  geom_histogram(color ="white",fill = ggcolors[3],bins = 60)+#scale_fill_manual(values = ggcolors[3:4])+
  stat_density(geom = "line",position = "identity",size=1,adjust = 1)+
  xlab("Co-expression pairs gain in Tumor(normalized)")+theme_classic()+
  geom_vline(xintercept = hypercirc.alter[4],color=ggcolors[1],size=1)+
  geom_vline(xintercept = nonDM.quant0.975,linetype="dashed",size=1)+
  annotate("text",x = c(nonDM.quant0.975+3),y = 0.6,label="alpha=0.05")+
  annotate("text",x = c(hypercirc.alter[4]-4),y = 0.6,label=paste0("Permutation test,\niteration = 10k\n",pval.text))+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,(hypercirc.alter[4]+5)))+gg_theme


#####readers###

#cor reader DEmRNA count
plot.df = as.data.frame(xtabs(~symbol+cor_status_tumor,data = reader.DEmRNA.allsamples.cor.df))%>%
  filter(cor_status_tumor=="positive")%>%arrange(desc(Freq))
plot.df$symbol = factor(plot.df$symbol,levels = plot.df$symbol)

DEmRNA.count.cor_reader.plot <- ggplot(plot.df,aes(x=symbol,y=Freq))+
  geom_bar(stat = "identity",width = 0.7,fill = ggcolors[4])+
  ylab("Count of Correlated DEmRNA")+xlab("")+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
        axis.text.x = element_text(size = 15 , color = "black",angle = 60,hjust = 1.1,vjust = 1.1),
        axis.text.y = element_text(size = 15 , color = "black"),
        legend.position = "none")


##fraction
pos.df = filter(reader.DEmRNA.allsamples.cor.df,cor_status_tumor =="positive")
plot.df = as.data.frame(prop.table(xtabs(~symbol+DEmRNA_circRNAs,data = pos.df),margin = 1))%>%dplyr::arrange(desc(Freq))

plot.df$symbol = factor(plot.df$symbol,levels = dplyr::arrange(filter(plot.df,DEmRNA_circRNAs=="cor_with_hyperCirc"),desc(Freq))$symbol)

DEmRNA.cor.reader.plot <- ggplot(plot.df,aes(x=symbol,y=Freq,fill = DEmRNA_circRNAs))+
  geom_bar(stat = "identity",position="stack",width = 0.7)+
  ylab("Fraction of Correlated DEmRNA")+xlab("")+
  theme_classic()+scale_fill_manual(values = ggcolors[5:6])+
  scale_y_continuous(expand = c(0, 0))+
  theme(axis.title = element_text(size = size_axistitle, color = "black"),
        axis.text.x = element_text(size = 15 , color = "black",angle = 60,hjust = 1.1,vjust = 1.1),
        axis.text.y = element_text(size = 15 , color = "black"),
        legend.position = "top",
        legend.title= element_blank(),
        legend.spacing.x = unit(0.1,"cm")
        )
DEmRNA.cor.reader.plot

#merge plot
fig4.plot.lst = list(m6a_nonm6a_cor.plot,
                     permut.plot,
                     DEmRNA.count.cor_reader.plot,
                     DEmRNA.cor.reader.plot)
plot_grid(plotlist = fig4.plot.lst,ncol = 2,scale = rep(0.85,4))

pdf("/data1/yeying/m6a_circ/final_plot/figure6/co-network.pdf",width = 9.92,height = 7.8)
plot_grid(plotlist = fig4.plot.lst,ncol = 2,scale = rep(0.85,4))
dev.off()

#save.image(file = "/data1/yeying/m6a_circ/script/10-17.RData")


######sup######
#####sup 2 : volcano plot for DEmRNAs####
source("/data1/yeying/m6a_circ/script/R_function.R")

mRNA.volcano.plot = volcano.plot( DE_res.df = dplyr::select(mRNA_res.df,ensemble_id,padj,log2FoldChange))


mRNA.volcano.plot 
grid.newpage()
grid.draw(total.TvsN.venn)

##venn : tumor pairs : m6A circRNA m6A mRNA
tumor.pairs.cor.df = filter(totalCirc_DEmRNA.allsample.cor.df,abs(cor_value_tumor) >= cor_cutoff)%>%
  left_join(mRNA_m6Astatus.lst,by = c("mRNA_id"="id"))
head(tumor.pairs.cor.df)
tumor.pairs.stats = as.data.frame(xtabs(~m6Astatus + mRNA_m6A ,data = tumor.pairs.cor.df))
# mRNA_m6A
# m6Astatus   m6A non-m6A
# m6A     14704    3901
# non-m6A  3297     882
fisher.test(xtabs(~m6Astatus + mRNA_m6A ,data = tumor.pairs.cor.df))

####
tmp = filter(totalCirc_DEmRNA.allsample.cor.df,m6Astatus=="m6A",!is.na(cor_value_tumor))%>%
  left_join(mRNA_m6Astatus.lst,by = c("mRNA_id"="id"))

m6Acirc_mRNAm6A_cor.wilcox.res = wilcox.test(filter(tmp,mRNA_m6A=="m6A")$cor_value_tumor,
                                             filter(tmp,mRNA_m6A=="non-m6A")$cor_value_tumor,alternative ="greater")
m6Acirc_mRNAm6A_cor.wilcox.res

rm(tmp)
m6Acirc_mRNAm6A_cor.plot <- 
  ggplot(filter(sample100k.totalCirc_DEmRNA.allsample.cor.df, m6Astatus=="m6A"))+
  stat_density(geom = "line",position = "identity",aes(x=cor_value_tumor,y=..density..,color=meth_type),size=0.8,adjust = 0.5)+
  xlab("Pearson Correlation")+ylab("Density")+
  theme_classic()+
  scale_color_manual(values = ggcolors)+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  geom_vline(xintercept = cor_cutoff)+
  ggtitle("Random m6AcircRNA-DEmRNA pairs(n = 100,000)")+
  annotate("rect",xmin = cor_cutoff,xmax = 1,alpha=0.2,fill=ggcolors[3], ymin = 0,ymax = 3)+
  geom_vline(xintercept = -cor_cutoff)+
  annotate("rect",xmin = -cor_cutoff,xmax = -1,alpha=0.2,fill=ggcolors[3], ymin = 0,ymax = 3)+
  annotate("text",x= -0.75,y=2,
           label=paste0("wilcox test ,\np = ",signif(m6Acirc_mRNAm6A_cor.wilcox.res$p.value,digits = 3)))+
  annotate("text",x= 0.75,y=2,label=paste0("N = ",tumor.pairs.stats[1,3]),color = ggcolors[1])+
  annotate("text",x= 0.75,y=1.8,label=paste0("N = ",tumor.pairs.stats[3,3]),color = ggcolors[2])+
  theme(
    axis.title = element_text(size = size_axistitle, color = "black"),
    axis.text.x = element_text(size = size_axistitle, color = "black"),
    axis.text.y = element_text(size = size_axistitle, color = "black"),
    legend.position = c(0.39,0.9),
    legend.title= element_blank(),
    legend.text = element_text(size = size_legendtext, face = NULL, color = "black"), 
    legend.background = element_rect(fill = NA))



tumor.pairs.mRNA = filter(mRNA_m6Astatus.lst,id%in%unique(tumor.pairs.cor.df$mRNA_id))
table(tumor.pairs.mRNA$mRNA_m6A)
tumor.pairs.circRNA = distinct(dplyr::select(tumor.pairs.cor.df,circRNA_id,m6Astatus))
table(tumor.pairs.circRNA$m6Astatus)

chisq.test(matrix(round(c(hypercirc.alter[2:3],nonDMcirc.alter[2:3])),nrow = 2,byrow = T,dimnames = list(c("hyper","nonDM"),c("tumor gain", "normal gain"))))

##cobind#####
reader.co_bind.pvalue.plot <- ggbarplot(data = hyper_fisher.df,x = "wer",y = "log10p",color = "white",
                                        fill = new_ggcolors[1],sort.val ="asc",orientation = c("horizontal"))+
  xlab("")+ylab("-log10(pvalue)")+scale_y_continuous(expand = c(0, 0),limits = c(0,300))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+gg_theme

reader.co_bind.pvalue.plot

sup.plot = list(mRNA.volcano.plot,m6Acirc_mRNAm6A_cor.plot,reader.co_bind.pvalue.plot)

pdf("/data1/yeying/m6a_circ/final_plot/figure6/sup_fig6.pdf",width = 9.4,height = 8.3)
plot_grid(plotlist = sup.plot,ncol = 2,scale = rep(0.9,3))
dev.off()
