rm(list = ls())

library(Seurat)
require(ggplot2)
require(reshape2)
library("org.Mm.eg.db")
library(EnsDb.Mmusculus.v79)
source("//LINLAB5//LabNAS-5//WY//scTFA//resource//R code//scTFA_related//scTFA_calculation_zscore_v20200709_1.R")

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
  axis.text.y = element_text(size=22)
)

data_dir="//LINLAB5//LabNAS-5//WY//scTFA//code//RAW 264.7//17CS//data//"
result_dir="//LINLAB5//LabNAS-5//WY//scTFA//code//RAW 264.7//17CS//code//check intron//"
date_analysis="analyzed20210201_1(combined)"
dir.create(paste(result_dir,date_analysis,sep=""))
setwd(paste(result_dir,date_analysis,sep=""))

#import RNA data
RNA_fpkm=read.csv(paste(data_dir,"GSE94383_fpms.csv",sep=""),row.names = 1)
colnames(RNA_fpkm)=sub("\\.","-",colnames(RNA_fpkm))

#import metadata
RNA_cell_id=read.csv(paste(data_dir,"GSE94383_cell_ids.csv",sep=""),row.names = 1)
RNA_cell_id=RNA_cell_id[RNA_cell_id$Condition %in% c("NoStim","Stim"),]
rownames(RNA_cell_id)=paste("X",rownames(RNA_cell_id),sep="")

#import movie data
RNA_movie_data=read.csv(paste(data_dir,"GSE94383_single_cell_dynamics.csv",sep=""),row.names = 1)
rownames(RNA_movie_data)=paste("X",rownames(RNA_movie_data),sep="")
RNA_movie_data$movie_time_points=apply(RNA_movie_data,1,function(x) sum(!is.na(x)))

#create seurat object
RAW_seurat=CreateSeuratObject(counts = RNA_fpkm, project = "RAW_fpkm", min.cells = 0, min.features = 0)

#add info of movie and cell id to metadata
#View(A549_seurat@meta.data)
RAW_seurat@meta.data=cbind(RAW_seurat@meta.data,RNA_cell_id[rownames(RAW_seurat@meta.data),],
                           RNA_movie_data[rownames(RAW_seurat@meta.data),])

#import GRN data
mouse_dorothea_logic_matrix=read.csv("//LINLAB5//LabNAS-5//WY//scTFA/resource/mouse/GRN/dorothea/mouse_GRN_classAB_activation_only/mouse_dorothea_logic_matrix classAB_activation_only_analyzed20200521_2.csv",
                                     row.names = 1,check.names=FALSE)

#Convert from gene.symbol to ensembl.gene
geneSymbols <-colnames(mouse_dorothea_logic_matrix)
geneIDs <- ensembldb::select(EnsDb.Mmusculus.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

mouse_dorothea_logic_matrix_ensembl=matrix(0,nrow(mouse_dorothea_logic_matrix),1,
                                           dimnames = list(rownames(mouse_dorothea_logic_matrix),"NA"))
for(col_i in colnames(mouse_dorothea_logic_matrix)){
  if(col_i %in% geneIDs$SYMBOL){
    gene_id_i=geneIDs[geneIDs$SYMBOL==col_i,"GENEID"]
    for(j_gene_id_i in gene_id_i){
      mouse_dorothea_logic_matrix_ensembl=cbind(mouse_dorothea_logic_matrix_ensembl,mouse_dorothea_logic_matrix[,col_i])
      colnames(mouse_dorothea_logic_matrix_ensembl)[ncol(mouse_dorothea_logic_matrix_ensembl)]=j_gene_id_i
    }
  }
}
mouse_dorothea_logic_matrix_ensembl=mouse_dorothea_logic_matrix_ensembl[,-1]

#import intron/exon CPM data and compute TFA
out_file_dir="//LINLAB5//LabNAS-5//WY//scTFA/Code/RAW 264.7/17CS/data/raw data/tpmcalculator_out_file/"
out_file_name_vector=list.files(out_file_dir)
out_file_data=list()
for(out_file_i in out_file_name_vector){
  data_i=read.table(paste(out_file_dir,out_file_i,sep=""),fill = T,header=T)
  data_i=na.omit(data_i)
  rownames(data_i)=sub("\\..*$","",data_i$Gene_Id)
  
  data_i$totalCPM=1e6*(data_i$Reads)/sum(data_i$Reads)
  data_i$ExonCPM=1e6*(data_i$ExonReads)/sum(data_i$Reads)
  data_i$IntronCPM=1e6*(data_i$IntronReads)/sum(data_i$Reads)
  
  TFA_ExonCPM=as.data.frame(scTFA_calculation_zscore_v20200709_1(as.matrix(data_i[,"ExonCPM",drop=F]),mouse_dorothea_logic_matrix_ensembl,zscore=F,
                                                                 gene_weight_method=NULL))
  colnames(TFA_ExonCPM)="TFA_ExonCPM"
  TFA_IntronCPM=as.data.frame(scTFA_calculation_zscore_v20200709_1(as.matrix(data_i[data_i$IntronLength>0,"IntronCPM",drop=F]),mouse_dorothea_logic_matrix_ensembl,zscore=F,
                                                                   gene_weight_method=NULL))
  colnames(TFA_IntronCPM)="TFA_IntronCPM"
  
  TFA_all=cbind(TFA_IntronCPM[rownames(TFA_IntronCPM),,drop=F],
                TFA_ExonCPM[rownames(TFA_IntronCPM),,drop=F])
  TFA_all$TF_name=rownames(TFA_all)
  out_file_data[[out_file_i]]=TFA_all
}

#change names
cell_sra_result=read.csv("//LINLAB5//LabNAS-5//WY//scTFA//code//RAW 264.7//17CS//data//raw data//sra_result.csv")
cell_SraRunInfo=read.csv("//LINLAB5//LabNAS-5//WY//scTFA//code//RAW 264.7//17CS//data//raw data//SraRunInfo.csv")
cell_metadata=merge(cell_sra_result,cell_SraRunInfo,all=F,by.x="Experiment.Accession",by.y="Experiment")
rownames(cell_metadata)=cell_metadata$Run

names_out_file_data_new=sub("_genes.out","",names(out_file_data))
names_out_file_data_new=cell_metadata[names_out_file_data_new,"Experiment.Title"]
names_out_file_data_new=sub("\\).*$","",sub("^.*\\(","",names_out_file_data_new))
names_out_file_data_new=paste("X",names_out_file_data_new,sep="")

names(out_file_data)=names_out_file_data_new
saveRDS(out_file_data,file=paste("out_file_data",".rds"))
# out_file_data=readRDS(paste(result_dir,date_analysis,"//out_file_data .rds",sep=""))

#subset
time_point_vector=c(75,150,300)
max_movie_time_point_vector=c(16,31,61)
TF_chosen=c("Rela")

correlation_matrix=matrix(,length(time_point_vector),3,
                          dimnames = list(paste("time",time_point_vector,sep="_"),c("IntronCPM_vs_movie","ExonCPM_vs_movie","cellnum")))
for(i_time_point in 1:length(time_point_vector)){
  time_point=time_point_vector[i_time_point]
  max_movie_time_point=max_movie_time_point_vector[i_time_point]
  RAW_subset=subset(RAW_seurat, subset = Time.point == time_point)
  
  RAW_subset=subset(RAW_subset, subset = nCount_RNA < 9e6 ) ##according to ???###
  
  out_file_data_subset=out_file_data[which(names(out_file_data) %in% colnames(RAW_subset))]
  RAW_subset=RAW_subset[,intersect(colnames(RAW_subset),names(out_file_data_subset))]
  
  #merge data
  index=0
  for(name_i in names(out_file_data_subset)){
    colnames(out_file_data_subset[[name_i]])=paste(colnames(out_file_data_subset[[name_i]]),name_i,sep="_")
    if(index==0){
      out_file_data_merged=out_file_data_subset[[name_i]]
    }else{
      out_file_data_merged=merge(out_file_data_merged,out_file_data_subset[[name_i]],all=F,
                                 by.x=paste("TF_name",names(out_file_data_subset)[1],sep="_"),
                                 by.y=paste("TF_name",name_i,sep="_"))
    }
    index=index+1
  }
  rownames(out_file_data_merged)=out_file_data_merged[,paste("TF_name",names(out_file_data_subset)[1],sep="_")]
  write.csv(out_file_data_merged,paste("out_file_data_merged",time_point,".csv"))
  
  TFA_ExonCPM_all=out_file_data_merged[,paste("TFA_ExonCPM",colnames(RAW_subset),sep="_")]
  colnames(TFA_ExonCPM_all)=colnames(RAW_subset)
  rownames(TFA_ExonCPM_all)=paste("TFA_ExonCPM",rownames(TFA_ExonCPM_all),sep="_")
  
  TFA_IntronCPM_all=out_file_data_merged[,paste("TFA_IntronCPM",colnames(RAW_subset),sep="_")]
  colnames(TFA_IntronCPM_all)=colnames(RAW_subset)
  rownames(TFA_IntronCPM_all)=paste("TFA_IntronCPM",rownames(TFA_IntronCPM_all),sep="_")
  
  #add metadata
  movie_TFA_data=cbind(RAW_subset@meta.data,
                       t(TFA_ExonCPM_all)[rownames(RAW_subset@meta.data),],
                       t(TFA_IntronCPM_all)[rownames(RAW_subset@meta.data),])
  
  index_vector=c("TFA_IntronCPM","TFA_ExonCPM")
  
  #correlation between TFA and nucl loc of every frame
  row_names_every=paste("corr with nucl of the",1:(max_movie_time_point+1),"frame before")
  corr_matrix_every=matrix(NA,length(row_names_every),length(index_vector),
                           dimnames = list(row_names_every,index_vector))
  for(index_i in index_vector){
    for(time_point_num in 1:(max_movie_time_point+1)){
      corr_matrix_every[paste("corr with nucl of the",time_point_num,"frame before"),index_i]=
        cor(movie_TFA_data[,paste(index_i,TF_chosen,sep="_")],
            movie_TFA_data[,paste("X",max_movie_time_point-time_point_num+1,sep="")],use="pairwise.complete.obs")
    }
  }
  corr_matrix_every=as.data.frame(corr_matrix_every)
  corr_matrix_every$label=as.numeric(sub(" frame before","",sub("corr with nucl of the ","",rownames(corr_matrix_every))))
  write.csv(corr_matrix_every,paste("corr_matrix_every data",time_point,".csv"))
  
  correlation_matrix[paste("time",time_point,sep="_"),c("IntronCPM_vs_movie","ExonCPM_vs_movie","cellnum")]=
    c(as.numeric(corr_matrix_every[corr_matrix_every$label==1,c("TFA_IntronCPM","TFA_ExonCPM")]),ncol(RAW_subset))
  
  #correlation between TFA and nucl loc of last frames
  row_names_last=paste("corr with sum nucl last",1:(max_movie_time_point+1),"frame",sep="_")
  corr_matrix_last=matrix(NA,length(row_names_last),length(index_vector),
                          dimnames = list(row_names_last,index_vector))
  for(index_i in index_vector){
    for(time_point_num in 1:(max_movie_time_point+1)){
      corr_matrix_last[paste("corr with sum nucl last",time_point_num,"frame",sep="_"),index_i]=
        cor(movie_TFA_data[,paste(index_i,TF_chosen,sep="_")],
            apply(movie_TFA_data[,paste("X",(max_movie_time_point-time_point_num+1):max_movie_time_point,sep=""),drop=F],1,sum),
            use="pairwise.complete.obs")
    }
  }
  corr_matrix_last=as.data.frame(corr_matrix_last)
  corr_matrix_last$label=sub("_frame","",sub("corr with sum nucl last_","",rownames(corr_matrix_last)))
  corr_matrix_last$label=as.numeric(corr_matrix_last$label)
  write.csv(corr_matrix_last,paste("corr_matrix_last data",time_point,".csv"))
}
correlation_matrix=as.data.frame(correlation_matrix)
correlation_matrix$time=factor(rownames(correlation_matrix),levels=rownames(correlation_matrix))
write.csv(correlation_matrix,paste("correlation_matrix data",".csv"))

correlation_matrix_plot=melt(correlation_matrix,id.vars = c("time","cellnum"),variable.name = "index",value.name = "corr")
write.csv(correlation_matrix_plot,paste("correlation_matrix_plot",".csv"))

p1<-ggplot(data=correlation_matrix_plot,mapping=aes(x=time,y=corr,fill=index,group=index))+
  geom_bar(stat="identity",position='dodge',width=0.5)+
  labs(x = "time", y = "Pearson correlation",
       title = paste("TF",TF_chosen))+
  geom_text(mapping=aes(x=time,y=corr+0.05,label=round(corr,3)),colour="black",position = position_dodge(0.5),size=6)+
  geom_text(mapping=aes(x=time,y=0,label=paste(cellnum,"\n","cell",sep="")),colour="black",position = position_dodge(0.5),size=6)+
  scale_y_continuous(limits=c(-0.03,0.53),breaks = seq(0,0.5,0.1))+
  theme_bw()+themo_demo

pdf(file =paste("correlation_matrix_plot",TF_chosen,".pdf"),width = 10,height = 8)
print(p1)
dev.off()

end_time=Sys.time()
capture.output(end_time-start_time,file=paste("runtime",date_analysis,".txt"))
capture.output(sessionInfo(),file=paste("sessionInfo",date_analysis,".txt"))
