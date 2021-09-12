monocle_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
setwd("/wtmp/user170/project/singlecell_Pal/25.cell-cycle/01.phloem_tra")
subsetscd<- monocle_cds[,rownames(subset(pData(monocle_cds),seurat_clusters %in% c("6","13","19")))]
HSMM<-subsetscd
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3 )


ordering_genes <- subset(dispersionTable(HSMM), mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
HSMM <- setOrderingFilter(HSMM, ordering_genes$gene_id)



HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM,reverse=TRUE)

write.table(pData(HSMM),file="HSMM_state.txt",sep="\t",quote=F)
write.table(HSMM@reducedDimS,file="HSMM_xy.txt",sep="\t",quote=F)

pdf("02.trajectory01.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters",)
dev.off()

pdf("02.trajectory02.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")+facet_wrap(~seurat_clusters, nrow = 1)
dev.off()

pdf("02.trajectory03.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

pdf("02.trajectory04.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "orig.ident")
dev.off()

##寻找分支基因

BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 6)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

branch_gene<-row.names(subset(BEAM_res, qval < 1e-4))
write.table(BEAM_res,file="BEAM_res.txt",sep="\t",quote=F)

##把分支基因画在图上
sink("branch_gene_cluster.txt")
pdf("02.branch_gene.pdf",width=4,height=4)
plot_genes_branched_heatmap(HSMM[branch_gene], branch_point = 1, num_clusters = 4,cores = 8,use_gene_short_name = F,show_rownames = F,return_heatmap=TRUE)
dev.off()
sink()

pdf("02.branch_gene.pdf",width=4,height=4)
plot_genes_branched_heatmap(HSMM[branch_gene], branch_point = 1, num_clusters = 4,cores = 8,use_gene_short_name = F)
dev.off()

#提出其中的基因分类，手动调整一下
grep -A 3000 "annotation_row" branch_gene_cluster.txt>branch_gene_cluster.txt2
mv branch_gene_cluster.txt2 branch_gene_cluster.txt
 perl ../../03.Pseudotime/get_function.pl ../../03.Pseudotime/pal_ptr_ath.txt branch_gene_cluster.txt>branch_gene_cluster_funct.txt


##展示某几个基因在主图上

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17659", use_color_gradient = TRUE)
p1$data$value[p1$data$value>4]<-4
pp1=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=2)+labs(title ="pou17659 (AUX1)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10))+theme(panel.grid=element_blank())


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou35162", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-6
pp2=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=3)+labs(title ="pou35162 (PIN1-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10))+theme(panel.grid=element_blank())


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17833", use_color_gradient = TRUE)
p1$data$value[p1$data$value>15]<-15
pp3=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2))+geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=7.5)+labs(title ="pou17833 (APL-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10))+theme(panel.grid=element_blank())



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05235", use_color_gradient = TRUE)
pp4=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=3)+labs(title ="pou05235 (FTIP1)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10))+theme(panel.grid=element_blank())



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou26862", use_color_gradient = TRUE)
p1$data$value[p1$data$value>50]<-50
pp5=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2))+geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=25)+labs(title ="pou26862 (SEOR1-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10))+theme(panel.grid=element_blank())

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou11679", use_color_gradient = TRUE)
p1$data$value[p1$data$value>20]<-20
pp6=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=10)+labs(title ="pou11679 (CalS7)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10))+theme(panel.grid=element_blank())

pdf("figure2B.pdf",width=8.27,height=11.7)
plot_grid(pp1,pp2,pp3,NA,pp4,pp5,pp6,ncol=4,nrow=5,hjust=1.5)
dev.off()

######韧皮发育附图
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05702", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp1=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou05702 OPS")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp1
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17169", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp2=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou17169 (SMXL3-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp2
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou32689", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp3=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou32689 (SMXL3-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp3
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou06647", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp4=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou06647 (SMXL3-c)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp4
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou24401", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp5=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou24401 (SMXL3-d)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp5
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05427", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp6=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou05427 (SMXL5-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp6
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou28037", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp7=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou28037 (SMXL5-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp7
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou19837", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp8=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou19837 (SUC2-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp8
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou31414", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp9=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou31414 (SUC2-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp9
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou01208", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp10=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou01208 (SEOR1-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp10
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou37783", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp11=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou37783 (SEOR1-c)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp11
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou01212", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp12=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou01212 (SEOR1-d)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp12
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou01209", use_color_gradient = TRUE)
p1$data$value[p1$data$value>50]<-50
pp13=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=25)+labs(title ="pou01209 (SEOR1-e)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp13
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou01211", use_color_gradient = TRUE)
p1$data$value[p1$data$value>50]<-50
pp14=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=25)+labs(title ="pou01211 (SEOR1-f)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp14
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou18979", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp15=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou18979 (SEOR1-g)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp15
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou21255", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp16=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou21255 CHER1")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp16
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou36012", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp17=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou36012 (SMXL6)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp17
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17683", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp18=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou17683 (SMXL7-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp18
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou02128", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp19=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou02128 (SMXL7-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp19
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou11099", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp20=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou11099 (D14)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp20
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou37525", use_color_gradient = TRUE)
p1$data$value[p1$data$value>20]<-20
pp21=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=10)+labs(title ="pou37525 (PP2-A1)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp21
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou37526", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp22=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou37526 (PP2-A10-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp22
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou34250", use_color_gradient = TRUE)
p1$data$value[p1$data$value>7]<-7
pp23=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=3.5)+labs(title ="pou34250 (PP2-A10-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp23
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou37527", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp24=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou37527 (PP2-A10-c)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp24
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou34248", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp25=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou34248 (PP2-A4-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp25
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou34293", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp26=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou34293 (PP2-A4-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp26
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou38918", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp27=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou38918 (ENODL1-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp27
dev.off()



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou26650", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp28=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou26650 (ENODL17-a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp28
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou04196", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp29=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou04196 (ENODL17-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp29
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou07937", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp30=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou07937 (ENODL17-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp31
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou26651", use_color_gradient = TRUE)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp31=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title ="pou26651 (ENODL20)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp31
dev.off()






p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou00500", use_color_gradient = TRUE)

p1$data$value[p1$data$value>25]<-25
pp32=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=12.5)+labs(title ="pou00500 (ENODL9-b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.2,0.3),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-10),plot.margin = margin(-0.7, 0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pp32
dev.off()

pdf("supp15-1.pdf",width=8.27,height=11.7)
plot_grid(pp1,pp6,pp7,pp8,pp9,pp10,pp11,pp12,ncol=4,nrow=5,hjust=1.5)
dev.off()

pdf("supp15-2.pdf",width=8.27,height=11.7)
plot_grid(pp13,pp14,pp15,pp16,pp17,pp18,pp19,pp20,ncol=4,nrow=5,hjust=1.5)
dev.off()

pdf("supp15-3.pdf",width=8.27,height=11.7)
plot_grid(pp21,pp22,pp23,pp24,pp25,pp26,pp27,pp28,ncol=4,nrow=5,hjust=1.5)
dev.off()

pdf("supp15-4.pdf",width=8.27,height=11.7)
plot_grid(pp29,pp30,pp31,pp32,ncol=4,nrow=5,hjust=1.5)
dev.off()



