rm(list = ls())
library(ggplot2)
library(pheatmap)
require(reshape2)
library(MetaCycle)
source("../../materials//scTFA_calculation.R")

themo_demo=theme(
  text = element_text(size=22),
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  axis.text=element_text(color='black'),
  plot.title = element_text(hjust = 0.5),
  title=element_text(size = 30),
  axis.title.x = element_text(size=24),
  axis.title.y = element_text(size=24),
  axis.text.x = element_text(size=22,angle = 45,hjust = 1),
  axis.text.y = element_text(size=22))

data_dir="../raw data/"
result_dir="../results/"

#import RNA data
RNA_data=read.table(paste(data_dir,"GSE73554_WT_AL_Intron_Exon_RFP.txt",sep=""),header=T)
RNA_RPKM=RNA_data
RNA_RPKM[,3:ncol(RNA_RPKM)]=2^RNA_RPKM[,3:ncol(RNA_RPKM)]#original data is log2(RPKM)

#aggregate rows with the same gene symbol by summing
RNA_RPKM_aggregated=aggregate(RNA_RPKM[,3:ncol(RNA_RPKM)],by=list(Gene_Symbol=RNA_RPKM$Gene_Symbol),FUN = sum)
rownames(RNA_RPKM_aggregated)=RNA_RPKM_aggregated$Gene_Symbol
RNA_RPKM_aggregated=RNA_RPKM_aggregated[,-1]
# colnames(RNA_RPKM_aggregated)

#split by time points and intron/exon
condition="WT_AL"
type=c("Intron","Exon")
time=c("00","02","04","06","08","10","12","14","16","18","20","22")

RNA_RPKM_split=matrix(,nrow(RNA_RPKM_aggregated),length(type)*length(time),
                      dimnames = list(rownames(RNA_RPKM_aggregated),paste(rep(type,each=length(time)),time,sep="_")))
for(type_i in type){
  for(time_i in time){
    RNA_RPKM_split[,paste(type_i,time_i,sep="_")]=apply(RNA_RPKM_aggregated[,grep(paste(condition,type_i,time_i,sep="_"),
                                                                                  colnames(RNA_RPKM_aggregated),value=T)],
                                                        1,mean)
  }
}
write.csv(RNA_RPKM_split,paste(result_dir,"RNA_RPKM_split",date_analysis,".csv",sep = ""))

# import GRN data
dorothea_logic_matrix=read.csv("../../materials/mouse_GRN_classAB_activation_only/dorothea_mouse_AB.csv",
                               row.names = 1, check.names = F)
dorothea_logic_matrix=as.matrix(dorothea_logic_matrix)

# TFA of RNA
GA_RNA=RNA_RPKM_split
GA_RNA=GA_RNA[rowSums(GA_RNA)>0,]

TFA_RNA=scTFA_calculation(as.matrix(GA_RNA),dorothea_logic_matrix,zscore=F,
                                     gene_weight_method=NULL)
write.csv(TFA_RNA,paste(result_dir,"TFA_RNA data",date_analysis,".csv",sep = ""))

TFA_RNA_Intron=TFA_RNA[,grep("Intron",colnames(TFA_RNA),value=T)]
write.csv(TFA_RNA_Intron,paste(result_dir,"TFA_RNA_Intron",date_analysis,".csv",sep = ""))

TFA_RNA_Exon=TFA_RNA[,grep("Exon",colnames(TFA_RNA),value=T)]
write.csv(TFA_RNA_Exon,paste(result_dir,"TFA_RNA_Exon",date_analysis,".csv",sep = ""))

#calculate circadian related parameters with MetaCycle
Intron_meta2d_result=meta2d(infile=paste(result_dir,"TFA_RNA_Intron",date_analysis,".csv",sep = ""), filestyle="csv", 
                            outdir=result_dir,
                            timepoints=as.numeric(sub("Intron_","",colnames(TFA_RNA_Intron))),
                            minper = 24,maxper = 24,ARSdefaultPer=24,
                            outputFile=F)
Intron_meta2d_result=Intron_meta2d_result$meta

Exon_meta2d_result=meta2d(infile=paste(result_dir,"TFA_RNA_Exon",date_analysis,".csv",sep = ""), filestyle="csv", result_dir,
                          timepoints=as.numeric(sub("Exon_","",colnames(TFA_RNA_Exon))),
                          minper = 24,maxper = 24,ARSdefaultPer=24,
                          outputFile=F)
Exon_meta2d_result=Exon_meta2d_result$meta

Intron_Exon_meta2d_result=merge(Intron_meta2d_result,Exon_meta2d_result,by.x="CycID",by.y="CycID",all=F,
                                suffixes=c("_Intron","_Exon"))
write.csv(Intron_Exon_meta2d_result,paste(result_dir,"Intron_Exon_meta2d_result",date_analysis,".csv",sep =""))

#choose TFs with circadian TFA from intron and exon (ARS p<0.05)
Intron_Exon_meta2d_result_subset=Intron_Exon_meta2d_result[(Intron_Exon_meta2d_result$meta2d_pvalue_Intron<0.05)&
                                                             (Intron_Exon_meta2d_result$meta2d_pvalue_Exon<0.05),]
Intron_Exon_circadian_TF_vector=Intron_Exon_meta2d_result_subset$CycID

phase_plot_data_frame=Intron_Exon_meta2d_result_subset[,c("CycID","meta2d_phase_Intron","meta2d_phase_Exon"),]
phase_plot_data_frame$Exon_minus_Intron_phase=phase_plot_data_frame$meta2d_phase_Exon-phase_plot_data_frame$meta2d_phase_Intron 

phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase<(-12)]=
  phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase<(-12)]+24
phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase>12]=
  phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase>12]-24

phase_plot_data_frame=phase_plot_data_frame[order(phase_plot_data_frame$Exon_minus_Intron_phase,decreasing = T),]
phase_plot_data_frame$CycID=factor(phase_plot_data_frame$CycID,levels=phase_plot_data_frame$CycID)

p1<-ggplot(data=phase_plot_data_frame,mapping=aes(x=CycID,y=Exon_minus_Intron_phase))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  labs(x = "Exon_minus_Intron_phase", y = "phase difference (h)",
       title = "")+
  theme_bw()+themo_demo

pdf(file =paste(result_dir,"phase_difference",date_analysis,".pdf",sep = ""),width = 10,height = 10)
print(p1)
dev.off()

TFA_RNA_Intron_subset=TFA_RNA[Intron_Exon_circadian_TF_vector,grep("Intron",colnames(TFA_RNA),value=T)]
TFA_RNA_Exon_subset=TFA_RNA[Intron_Exon_circadian_TF_vector,grep("Exon",colnames(TFA_RNA),value=T)]

# scale per TF
TFA_RNA_Intron_subset_scale=t(apply(TFA_RNA_Intron_subset,1,scale))
colnames(TFA_RNA_Intron_subset_scale)=colnames(TFA_RNA_Intron_subset)
TFA_RNA_Exon_subset_scale=t(apply(TFA_RNA_Exon_subset,1,scale))
colnames(TFA_RNA_Exon_subset_scale)=colnames(TFA_RNA_Exon_subset)

TFA_all_scale=as.data.frame(cbind(TFA_RNA_Intron_subset_scale,TFA_RNA_Exon_subset_scale))
TFA_all_scale$TF=rownames(TFA_all_scale)
write.csv(TFA_all_scale,paste(result_dir,"TFA_all_scale",".csv",sep = ""))

TFA_all_scale_melt=melt(TFA_all_scale,id.vars = "TF",variable.name = "index",value.name ="scaled_TFA" )
TFA_all_scale_melt$normalization=factor(sub("_.*$","",TFA_all_scale_melt$index),levels=c("Intron","Exon"))
TFA_all_scale_melt$time=as.numeric(sub("^.*_","",TFA_all_scale_melt$index))
TFA_all_scale_melt$TF=factor(TFA_all_scale_melt$TF,levels=phase_plot_data_frame$CycID)
write.csv(TFA_all_scale_melt,paste(result_dir,"TFA_all_scale_melt",".csv",sep = ""))

p_TFA_all_scale_melt<-ggplot(data=TFA_all_scale_melt,mapping=aes(x=time,y=scaled_TFA,colour=normalization))+
  geom_line(size=1)+
  geom_point(size = 2)+
  labs(x = "Time (h)", y = "scaled_TFA",
       title ="")+
  scale_x_continuous(limits=c(0,22),breaks=seq(0,22,2))+
  facet_grid(TF~.)+
  theme_bw()+themo_demo

pdf(paste(result_dir,"p_TFA_all_scale_melt",date_analysis,".pdf",sep = ""),height=12,width=5.5)
print(p_TFA_all_scale_melt)
dev.off()
