#!>=R4.0
rm(list = ls())
setwd("/Users/fcb/Desktop/Human TF dynamics from multiomics/resource/mouse/GRN/dorothea/mouse_GRN_classAB_activation_only")
library(dorothea)
data(dorothea_mm, package = "dorothea")
dorothea_mm_selected=dorothea_mm[which((dorothea_mm$confidence %in% c("A","B"))*
                                         (dorothea_mm$mor==1)==1),]

TF_list_dorothea=names(table(dorothea_mm_selected[,"tf"]))
TF_list_num=length(TF_list_dorothea)

target_list_dorothea=names(table(dorothea_mm_selected[,"target"]))
target_list_num=length(target_list_dorothea)

dorothea_logic_matrix=matrix(0,TF_list_num,target_list_num)
rownames(dorothea_logic_matrix)=TF_list_dorothea
colnames(dorothea_logic_matrix)=target_list_dorothea

for(i_dorothea_logic_matrix in TF_list_dorothea)
{
  for(j_dorothea_logic_matrix in target_list_dorothea)
  {
    if(sum((dorothea_mm_selected$tf==i_dorothea_logic_matrix)*(dorothea_mm_selected$target==j_dorothea_logic_matrix))>0)
    {
        dorothea_logic_matrix[i_dorothea_logic_matrix,j_dorothea_logic_matrix]=1
    }
  }
}

write.csv(dorothea_logic_matrix,paste("mouse_dorothea_logic_matrix classAB_activation_only_analyzed20200521_2.csv",sep=""))
