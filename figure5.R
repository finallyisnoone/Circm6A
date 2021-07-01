
########data from  9-23 methy bio##########
#######DM circRNA expression status####
#data from 9-23 meth_bio###

meth_expre_scatter.plot <- ggplot()+geom_point(data = filter(final_m6acirc_meth_expres_res.df2,grepl("nonSig",label)),
                    aes_string(y= "log2FC_express" ,x="log2meth"),size = 2,color="grey50")+
  geom_point(data =color.df ,aes_string(y= "log2FC_express" ,x="log2meth",color = "label"),size = 2)+
  xlim(-max(abs(final_m6acirc_meth_expres_res.df2$log2meth)),max(abs(final_m6acirc_meth_expres_res.df2$log2meth)))+
  ylim(-max(abs(final_m6acirc_meth_expres_res.df2$log2FC_express)),max(abs(final_m6acirc_meth_expres_res.df2$log2FC_express)))+
  geom_vline(xintercept = lfc,linetype="dashed")+geom_vline(xintercept = -lfc,linetype="dashed")+
  geom_hline(yintercept = -lfc,linetype="dashed")+geom_hline(yintercept  = lfc,linetype="dashed")+
  ylab("log2[(SPRBM_T+0.1)/(SPRBM_N+0.1)]")+xlab("log2[(m6A_T+0.1)/(m6A_N+0.1)]")+
  scale_color_manual(values = ggcolors)+
  theme_test()+gg_theme

meth_expre_scatter.plot
########hyper_up.reader#######
hyper_up.reader.ecdf.plot.df = rbind(dplyr::select(circRNA_wers.cor.df,cor_value)%>%mutate(label="total circRNAs"),
                                     data.frame(cor_value = hyper_upcirc.reader.cor.anno.df$cor_value,
                                                label = "hyper-up CircRNAs"))#hyper_upcirc.reader.cor.anno.df$symbol))

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

#width=536&height=418
pdf(file =  "/data1/yeying/m6a_circ/final_plot/figure5/meth_expres_scatter.pdf",width = 5.36,height = 4.18)
meth_expre_scatter.plot
dev.off()

#492&height=426&scale=1
pdf("/data1/yeying/m6a_circ/final_plot/figure5/figure5_add-10-15.pdf",width = 5,height = 4.3)
meth_expre_scatter.plot
hyper_up.reader.density.plot
#plot_grid(meth_expre_scatter.plot,hyper_up.reader.density.plot)
dev.off()

