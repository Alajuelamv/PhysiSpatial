library("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library(dplyr)
library(tidyjson)
# res_transposed <- read.csv("/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/Pipelines/uwu.csv")
# res_transposed <- read.csv("/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/pruebas/Stutilityvsscanpy/scanpy/uwu.csv")
res_transposed <- read.csv("input_Ensembl2Entrez.csv")
res_transposed$Entrez_Gene_Id <- mapIds(org.Hs.eg.db, keys = res_transposed$Entrez_Gene_Id, keytype="ENSEMBL", column = "ENTREZID")
colnames(res_transposed)[1] <-"Hugo_Symbol"

# este apply es para la ultima row que marca el cluster al cual nos referimos pero ya veremos si podemos tal uwu
# res_transposed <- res_transposed[-nrow(res_transposed),] no lo ponemos porque para el setup lo necesitamos
final <- dim(res_transposed)[1]

res_transposed1 <- res_transposed[1:10000,]
res_transposed2 <- res_transposed[10001:final,]

res_transposed1 <- apply(res_transposed1,2,as.character) # esto importante para poder escribir la tabla sino hay algunas columnas que son lista que no deja
res_transposed2 <- apply(res_transposed2,2,as.character)

res_transposed <- rbind(res_transposed1,res_transposed2)

# write.table(res_transposed,file = "PROFILE-master/Data/Visium/uwu2.txt",sep="\t",row.names = FALSE,quote = FALSE)
write.table(res_transposed,file = "input_ESTIMATE.txt",sep="\t",row.names = FALSE,quote = FALSE)

