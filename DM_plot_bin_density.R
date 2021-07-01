
#outdir = "/data1/yeying/m6a_circ/analysis/m6A_to_circRNA/biogenesis/intron_region/density_plot/region_2k/bin_100bp/"
outdir = "/data1/yeying/m6a_circ/pancreatic/DM/biogenesis"
bin_length = 50

gg_theme2 = theme(
  axis.title = element_text(size = size_axistitle, color = "black"),
  axis.text.x = element_text(size = size_axisgreytext, color = "black"),
  axis.text.y = element_text(size = size_axisblacktext, color = "black"),
  legend.key.size = unit(0.5,"cm"),
  legend.title= element_blank(),
  legend.text = element_text(size = 12, face = NULL, color = "black"), 
  legend.background = element_rect(fill = NA))

#######################calculate density data.frame######################
#####m6A
setwd(paste0(outdir,"/","m6A"))
setwd(paste0(outdir,"/","hyper"))
m6A.cov.file = list.files("./",pattern = "*cov.txt")
m6A.cov.file
m6A.cov.df = data.frame()
for ( i in 1:length(m6A.cov.file) ){
  splice_site = strsplit(m6A.cov.file[i],"_")[[1]][3]
  print(splice_site)
  tmp.df = read_delim(m6A.cov.file[i],delim = "\t")%>%
    mutate(m6Astatus = "m6A",splice_site = splice_site,
           bin_final = ifelse(splice_site=="3ss",-bin_length*(bin-1),bin_length*(bin-1)))
  tmp.df$coverage_density = tmp.df$coverage/sum(tmp.df$coverage)
  m6A.cov.df = rbind(m6A.cov.df,tmp.df)
}


m6A.cov.stats.df = as.data.frame(group_by(m6A.cov.df,tag,bin_final,m6Astatus)%>%
                                   dplyr::summarise(total_cov = sum(coverage),total_cov_density = sum(coverage_density)))

sum(filter(m6A.cov.stats.df,tag=="ADAR1")$total_cov_density)
head(m6A.cov.stats.df)
######non-m6A
setwd(paste0(outdir,"/","non-m6A"))
nonm6A.cov.file = list.files("./",pattern = "*cov.txt")
nonm6A.cov.file
nonm6A.cov.df = data.frame()
for ( i in 1:length(nonm6A.cov.file) ){
  splice_site = strsplit(nonm6A.cov.file[i],"_")[[1]][3]
  print(splice_site)
  tmp.df = read_delim(nonm6A.cov.file[i],delim = "\t")%>%
    mutate(m6Astatus = "non-m6A",splice_site = splice_site,
           bin_final = ifelse(splice_site=="3ss",-bin_length*(bin-1),bin_length*(bin-1)))
  tmp.df$coverage_density = tmp.df$coverage/sum(tmp.df$coverage)
  nonm6A.cov.df = rbind(nonm6A.cov.df,tmp.df)
}

nonm6A.cov.stats.df = as.data.frame(group_by(nonm6A.cov.df,tag,bin_final,m6Astatus)%>%
                                      dplyr::summarise(total_cov = sum(coverage),total_cov_density = sum(coverage_density)))
head(nonm6A.cov.stats.df)
#############merge 
tmp = data.frame(tag = rep(unique(nonm6A.cov.df$tag),2),bin_final = 0,
                 m6Astatus=c(rep("m6A",4),rep("non-m6A",4)),total_cov = 0,total_cov_density = 0)
tmp
total_circ.cov.stats.df = rbind(m6A.cov.stats.df,nonm6A.cov.stats.df)%>%mutate(group_tag = factor(paste0(m6Astatus,"_",tag)))
head(total_circ.cov.stats.df)
tmp_sum = sum(total_circ.cov.stats.df$total_cov)
total_circ.cov.stats.df$group_tag = factor(total_circ.cov.stats.df$group_tag,
                                           levels = levels(total_circ.cov.stats.df$group_tag))
str(total_circ.cov.stats.df)




#####################plot###########################

#####function

gg_line <- function(line.df,tag_vec,y_var="total_cov_density"){
  
  p = ggplot(filter(line.df,tag%in%tag_vec),
             aes_string(x="bin_final",y=y_var,color = "tag",linetype="m6Astatus"))+
               geom_line(size = 1.5)+scale_color_manual(values = ggcolors[1:4] )+xlab("")+ylab("RBP binding density")#+gg_theme2
  
  return(p)
}


plot_bin_density <- function(density.df, tag_str,y_var ="total_cov_density",lineColor=ggcolors[1],y_tag ="density"){
  p = ggplot(filter(density.df,tag%in%tag_str),
         aes_string(x="bin_final",y=y_var,linetype="m6Astatus"))+
    geom_line(size = 1.5,color = lineColor)+#scale_color_manual(values = lineColor )+
    xlab("")+ylab(paste0(tag_str," binding ",y_tag))+theme(
      axis.title = element_text(size = size_axistitle, color = "black"),
      axis.text.x = element_text(size = size_axisgreytext, color = "black"),
      axis.text.y = element_text(size = size_axisblacktext, color = "black"),
      legend.position = "none")
  return(p)
}

#############final coverage#################
# #########merge RBP&TE###########
# tmp.lst = list(gg_line(total_circ.cov.stats.df,unique(total_circ.cov.stats.df$tag),y_var="total_cov")
#    ,gg_line(total_circ.cov.stats.df,unique(total_circ.cov.stats.df$tag),y_var="total_cov_density"))
# 
# plot_grid(plotlist = tmp.lst,nrow = 1)


########final coverage plot###########
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
total_cov.density.plot
total_circ.cov.stats.df

wilcox.rbps.res.df = data.frame()

for (i in unique(total_circ.cov.stats.df$tag)){
  wilcox.rbps = wilcox.test(dplyr::filter(total_circ.cov.stats.df,tag==i, m6Astatus=="m6A")$total_cov ,
                            dplyr::filter(total_circ.cov.stats.df,tag==i, m6Astatus=="non-m6A")$total_cov)
  tmp.vec = c(elements=i,w = wilcox.rbps$statistic,pvalue=wilcox.rbps$p.value)
  wilcox.rbps.res.df = rbind(wilcox.rbps.res.df,tmp.vec)
}

wilcox.rbps.res.df
outdir
# 1       TE     2    1.76563312724843e-14
# 2      QKI    16    4.70481369400077e-14
# 3     DHX9     0    1.48652279718363e-14
# 4    ADAR1    28    1.24441891197156e-13


pdf(paste0(outdir,"/hyper/circular_factor_cov.density.pdf"),width = 7,height = 10)
total_cov.density.plot
dev.off()

