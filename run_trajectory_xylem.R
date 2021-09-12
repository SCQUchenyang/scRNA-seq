##先看看HDZIP-III的表达
setwd("/wtmp/user170/project/singlecell_Pal/25.cell-cycle/02.xylem_tra/01.HDZIP")
test<-read.table('pal_HDzipIII.list',header=F)

pdf("pal_HDzipIII.pdf")
DotPlot(combined, features = test$V1,col.min=0)+RotatedAxis()+scale_x_discrete("")+scale_y_discrete("")
dev.off()

####开始拟时分析
data <- as(as.matrix(combined@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = combined@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
setwd("/wtmp/user170/project/singlecell_Pal/25.cell-cycle/02.xylem_tra")

subsetscd<- monocle_cds[,rownames(subset(pData(monocle_cds),seurat_clusters %in% c("5","9","8","4","2")))]
HSMM<-subsetscd
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 3 )


ordering_genes <- subset(dispersionTable(HSMM), mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
HSMM <- setOrderingFilter(HSMM, ordering_genes$gene_id)



HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM)
HSMM <- orderCells(HSMM,root_state=2)
write.table(pData(HSMM),file="HSMM_state.txt",sep="\t",quote=F)
write.table(HSMM@reducedDimS,file="HSMM_xy.txt",sep="\t",quote=F)

pdf("02.trajectory01.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters",theta=160)
dev.off()

pdf("02.trajectory02.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "seurat_clusters",theta=160)+facet_wrap(~seurat_clusters, nrow = 1)
dev.off()

pdf("02.trajectory03.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",theta=160)
dev.off()

pdf("02.trajectory04.pdf",width=4,height=4)
plot_cell_trajectory(HSMM, color_by = "orig.ident",theta=160)
dev.off()


##展示ACL5 HB家族基因表达

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou06855", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>15]<-15
pp1=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=7.5)+labs(title ="pou06855 (ACL5a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou18623", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>50]<-50
pp2=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=25)+labs(title ="pou18623 (ACL5b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou23052", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>6]<-6
pp3=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=3)+labs(title ="pou23052 (ACL5c)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou39865", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>20]<-20
pp4=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=10)+labs(title ="pou39865 (ACL5d)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou07035", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>3]<-3
pp5=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=1.5)+labs(title ="pou07035 (PalHB8)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou28995", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>8]<-8
pp6=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=4)+labs(title ="pou28995 (PalHB7)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou00919", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>10]<-10
pp7=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=5)+labs(title ="pou00919 (PalHB4a)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pdf("test.pdf")
pp7
dev.off()


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou00937", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>10]<-10
pp8=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=5)+labs(title ="pou00937 (PalHB4b)")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())


pdf("figure3B.pdf",width=8.27,height=11.7)
plot_grid(pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,ncol=4,nrow=5,hjust=1.5)
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
perl ../../03.Pseudotime/get_function.pl branch_gene_cluster.txt2 BEAM_res.txt>branch_gene_cluster.txt

perl ../../03.Pseudotime/get_function.pl branch_gene_cluster.txt ../../03.Pseudotime/pal_ptr_ath.txt >branch_gene_cluster_funct.txt



###fiber和vessel的maker
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou14792", use_color_gradient = TRUE,theta=160)
sort(p1$data$value)
p1$data$value[which(p1$data$value>8)]=8
pp1=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=4)+labs(title ="pal-pou14792 PdDUF579")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())
pp1

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou20622", use_color_gradient = TRUE,theta=160)
sort(p1$data$value)
p1$data$value[which(p1$data$value>70)]=70
pp2=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=35)+labs(title ="pal-pou20622 PtAP66")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())
pp2


p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou12358", use_color_gradient = TRUE,theta=160)
sort(p1$data$value)
p1$data$value[which(p1$data$value>30)]=30
pp3=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=15)+labs(title ="pal-pou12358 PtAP17")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())
pp3



p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou14427", use_color_gradient = TRUE,theta=160)
sort(p1$data$value)
p1$data$value[which(p1$data$value>10)]=10
pp4=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=5)+labs(title ="pal-pou14427 XCP1")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())
pp4




p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou30825", use_color_gradient = TRUE,theta=160)
p1$data$value[which(p1$data$value>10)]=10

pp5=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=5)+labs(title ="pal-pou30825 VND1")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())
pp5

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05808", use_color_gradient = TRUE,theta=160)
p1$data$value[which(p1$data$value>8)]=8
pp6=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=4)+labs(title ="pal-pou05808 MAN6")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())
pp6
dev.off()

pdf("figure3E.pdf",width=8.27,height=11.7)
plot_grid(pp1,pp2,pp3,pp4,pp5,pp6,ncol=4,nrow=5,hjust=1.5)
dev.off()



##两个分支各选一个基因来验证
pou23901
pou17689
pou15262

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou23901", use_color_gradient = TRUE,theta=160)
p1$data$value[which(p1$data$value>90)]=90
pp1=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=45)+labs(title ="pal-pou23901 ESK1")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17689", use_color_gradient = TRUE,theta=160)
p1$data$value[which(p1$data$value>16)]=16
pp2=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=8)+labs(title ="pal-pou17689 ESK1")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou15262", use_color_gradient = TRUE,theta=160)
p1$data$value[which(p1$data$value>110)]=110
pp3=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value))+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=55)+labs(title ="pal-pou15262 ABR1")+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 9,hjust =0.70,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank())

pdf("ESK_ABR1.pdf",width=8.27,height=11.7)
plot_grid(pp1,pp2,pp3,ncol=4,nrow=5,hjust=1.5)
dev.off()



####SCW相关基因
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou01999", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp1=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou01999 ("MYB83a/PtrMYB3")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp1
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou36168", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp2=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou36168 ("MYB83b/PtrMYB20")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp2
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou02068", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp3=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou02068 ("MYB46a/PtrMYB2")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp3
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou36096", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp4=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou36096 ("MYB46b/PtrMYB021")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp4
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou06182", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp5=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou06182 ("NAC075a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp5
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou28740", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp6=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou28740 ("NAC075b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp6
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou12752", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp7=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou12752 ("SND2a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp7
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou21824", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp8=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou21824 ("SND2b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp8
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou27315", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp9=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou27315 ("PtrSND2/PtrNAC154")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp9
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou30579", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp10=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou30579 ("SND2d/PtrNAC156")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp10
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou03491", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp11=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou03491 ("KNAT7")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp11
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou06391", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp12=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou06391 ("CesA7a/PtrCesA7a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp12
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou28337", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp13=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou28337 ("CesA7b/PtrCesA7b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp13
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou28338", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp14=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou28338 ("CesA7c/PtrCesA7c")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp14
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou12874", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp15=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou12874 ("CesA8b/PtrCesA8b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp15
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou21698", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp16=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou21698 ("CesA8c/PtrCesA8c")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp16
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou21759", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp17=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou21759 ("CesA8a/PtrCesA8a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp17
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou02166", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp18=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou02166 ("LAC4a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp18
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05706", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp19=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou05706 ("LAC4b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp19
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17636", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp20=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou17636 ("LAC4c")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp20
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou33798", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp21=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou33798 ("LAC4d")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp21
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou33800", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp22=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou33800 ("LAC4e")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp22
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou35970", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp23=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou35970 ("LAC4f")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp23
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05594", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp24=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou05594 ("LAC17a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp24
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05601", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp25=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou05601 ("LAC17b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp25
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou21262", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp26=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou21262 ("LAC17c")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp26
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou03934", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp27=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou03934 ("IRX10a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp27
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou08134", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp28=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou08134 ("IRX10b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp28
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou00477", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp29=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou00477 ("IRX8a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp29
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou21116", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp30=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou21116 ("IRX8b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp30
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou14790", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp31=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou14790 ("IRX14-La")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp31
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou29629", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp32=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou29629 ("IRX14-Lb")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp32
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou29630", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp33=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou29630 ("IRX14-Lc")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp33
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou05996", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp34=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou05996 ("IRX9a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp34
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou32547", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp35=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou32547 ("IRX9b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp35
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou17689", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp36=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou17689 ("ESK1a")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp36
dev.off()
p1=plot_cell_trajectory(HSMM,color_by = "State",markers = "pal-pou23901", use_color_gradient = TRUE,theta=160)
p1$data$value[p1$data$value>quantile(p1$data$value,0.99)]<-quantile(p1$data$value,0.99)
pp37=ggplot(p1$data,aes(x=data_dim_1,y=data_dim_2)) +geom_point(aes(color=value),size=0.5)+scale_color_gradient2(low="lightgrey",mid="orange",high="#4B0082",midpoint=(quantile(p1$data$value,0.99))/2)+labs(title =~italic(pou23901 ("ESK1b")))+theme_bw()+xlab("")+ ylab("")+theme(axis.text=element_blank(),axis.ticks= element_blank(),legend.position=c(0.8,0.7),legend.key.width = unit(0.15, "cm"),legend.key.height = unit(0.3, "cm"),legend.title=element_blank(),legend.text=element_text(size=7),plot.title = element_text(size = 8,hjust =0.60,vjust=-60),plot.margin = margin(-0.4,0.2, -0.2,-0.4, "cm"))+theme(panel.grid=element_blank(),legend.position = 'none')

pp37
dev.off()

pdf("supp14-SCW.pdf",width=8.27,height=11.7)
plot_grid(pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp9,pp10,pp11,pp12,pp13,pp14,pp15,pp16,pp17,pp18,pp19,pp20,pp21,ncol=4,nrow=5,hjust=1.5)
plot_grid(pp22,pp23,pp24,pp25,pp26,pp27,pp28,pp29,pp30,pp32,pp34,pp35,pp36,pp37,ncol=4,nrow=5,hjust=1.5)
dev.off()


