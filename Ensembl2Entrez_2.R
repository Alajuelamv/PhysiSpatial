library("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library(dplyr)
library(tidyjson)
# res_transposed <- read.csv("/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/Pipelines/uwu.csv")
# res_transposed <- read.csv("/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/pruebas/Stutilityvsscanpy/scanpy/uwu.csv")
res_transposed <- read.csv("uwu5.txt",sep = "\t")
res_transposed$Entrez_Gene_Id <- mapIds(org.Hs.eg.db, keys = res_transposed$Entrez_Gene_Id, keytype="ENSEMBL", column = "ENTREZID")
colnames(res_transposed)[1] <-"Hugo_Symbol"

# write.table(res_transposed,file = "PROFILE-master/Data/Visium/uwu2.txt",sep="\t",row.names = FALSE,quote = FALSE)
write.table(res_transposed,file = "input_PROFILE.txt",sep="\t",row.names = FALSE,quote = FALSE)

