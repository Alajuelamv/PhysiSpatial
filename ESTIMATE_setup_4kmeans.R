library(estimate)

filterCommonGenes(input.f="input_ESTIMATE.txt", output.f="visium_breast.gct", id="GeneSymbol")

estimateScore("visium_breast.gct", "visium_estimate_score.gct", platform="affymetrix")

plotPurity(scores="visium_estimate_score.gct", samples="all_samples", platform="affymetrix",output.dir ="/home/alejandro/Escritorio/MasterBioinformatica/TFM/TFM2022AlejandroMadrid/Spatial2Physi/output_estimate" )

library(cluster)
library(CePa)
library("dplyr")
library(ggfortify)

clusters_raw <- read.csv("input_ESTIMATE.txt", sep="\t")

estimate <- read.gct("visium_estimate_score.gct")
estimate <- as.data.frame(estimate)
clusters <- as.data.frame(tail(clusters_raw, n=1))
clusters <- clusters[,-c(1,2)]

estimateT <- t(estimate)
clustersT <- t(clusters)
data_frame_merge <- merge(estimateT, clustersT, by = 'row.names', all = TRUE)

colnames(data_frame_merge)[6] <- "Cluster"

pca_res <- prcomp(data_frame_merge[2:5], scale. = TRUE)
autoplot(pca_res, data = data_frame_merge, colour = "Cluster")

summary(pca_res)

library(cluster)
library(factoextra)

iris_transform = as.data.frame(-pca_res$x[,1:2])
kmeans_iris = kmeans(iris_transform, centers = 4, nstart = 50)
fviz_cluster(kmeans_iris, data = iris_transform)

data_frame_merge$kmeans <- kmeans_iris$cluster

meta_clusters <- data_frame_merge %>% select(Cluster,kmeans)

table <- table(meta_clusters$Cluster,meta_clusters$kmeans)

set_up <- data_frame_merge %>% select(Row.names,kmeans)

write.csv(set_up, file = "uwu4_1.csv",row.names = FALSE,quote = FALSE,col.names = FALSE)
