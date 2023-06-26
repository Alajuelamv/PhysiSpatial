original_path <- normalizePath(getwd())
base_folder <- dirname(dirname(original_path))

knitr::opts_knit$set(root.dir = base_folder)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

if (!require("pacman")) install.packages("pacman")
list.of.packages <- c("tidyverse","VennDiagram","networkD3","knitr","gridExtra","mclust","diptest","moments","magrittr")
pacman::p_load(list.of.packages, character.only = TRUE)

tbl_to_df <- function(x){
  x <- x %>% as.data.frame %>%  remove_rownames %>% column_to_rownames(var="PATIENT_ID")
} 
rename <- dplyr::rename
select <- dplyr::select

#First import META clinical data for further use to correlate with logical modelling results
META_clin <- read_delim("PROFILE-master/Data/METABRIC/data_clinical_supp_patient.txt", delim = "\t") %>% select(PATIENT_ID, OS_STATUS, OS_MONTHS, CLAUDIN_SUBTYPE) %>% rename(PAM50=CLAUDIN_SUBTYPE) %>% mutate(PAM50=factor(PAM50))

#Due to inconsistencies in genes names (DEC1 replaced by 1-Dec for instance), Entrez ID are used to recover proper HUGO names. Import HUGO/Entrez table for later qualiity check
#Some tables have to be transposed
HUGO_Entrez <- read_delim("PROFILE-master/Data/Common/HUGO_Entrez.txt", delim = "\t") %>% na.omit
entrez_to_hugo <- function(df_input){
  df_input <- df_input %>% mutate(Hugo_Symbol=unlist(map2(.$Hugo_Symbol,.$Entrez_Gene_Id, function(x,y) if(y %in% HUGO_Entrez$`Entrez Gene ID` & !(y %in% .$Entrez_Gene_Id[duplicated(.$Entrez_Gene_Id) | duplicated(.$Entrez_Gene_Id, fromLast = T)])) HUGO_Entrez$`Approved Symbol`[which(HUGO_Entrez$`Entrez Gene ID`==y)] else x ))) 
  dupli <- duplicated(df_input$Hugo_Symbol) | duplicated(df_input$Hugo_Symbol, fromLast = T)
  df_input$Hugo_Symbol[dupli] <- paste(df_input$Hugo_Symbol[dupli],df_input$Entrez_Gene_Id[dupli],sep = "-")
  df_output <- df_input %>% select(-Entrez_Gene_Id) %>% select(-matches("^X[0-9]+$"))
}

tibble_transpose <- function(df_input){
  df_output <- df_input %>% gather(var, value, -Hugo_Symbol) %>%
    spread(Hugo_Symbol, value) %>%
    rename(PATIENT_ID=var) %>% 
    type_convert
}

# #Import META mutations data
# META_mut <- read_delim("PROFILE-master/Data/METABRIC/data_mutations_extended.txt", delim = "\t", skip = 2) %>% select(Tumor_Sample_Barcode, Hugo_Symbol,Variant_Classification, `MA:protein.change`, SIFT, PolyPhen)  %>% rename(PATIENT_ID=Tumor_Sample_Barcode)
# 
# META_mut_patients <- readLines("PROFILE-master/Data/METABRIC/data_mutations_extended.txt",n=1) %>% strsplit(split=" ") %>% unlist %>% tail(length(.)-1)
# 
# #Import and process other META omics data
# META_CNA <- read_delim("PROFILE-master/Data/METABRIC/data_CNA.txt", delim = "\t") %>% entrez_to_hugo %>% tibble_transpose


META_RNA <- read_delim("input_PROFILE.txt", delim = "\t") %>% entrez_to_hugo %>% tibble_transpose

#Additional imports: PAM50 gene list
PAM50 <- read_delim("PROFILE-master/Data/Common/pam50_centroids.txt", delim = "\t") %>% rename(Gene=X1) %>% select(Gene)

#Genes involved in Fumia model
genenames <- read.table("PROFILE-master/Models/Fumia2013/Fumia_namesToHugo_curated.txt",header=T,sep="\t")
geneindex <- strsplit(as.character(genenames[,2]), split = ",") %>% sapply(function(l){gsub(" ","",l)})
geneindex <- data.frame(V1 = rep(genenames[,1], sapply(geneindex, length)), V2 = unlist(geneindex))
model_nodes_HUGO <- unique(geneindex[,2]) %>% sub("^\\s+", "", .)

#Create new variables with only genes related to the model
# METAmodel_mut <- META_mut %>% filter(Hugo_Symbol %in% model_nodes_HUGO)
# METAmodel_CNA <- META_CNA %>% select(PATIENT_ID, one_of(model_nodes_HUGO))
METAmodel_RNA <- META_RNA %>% select(PATIENT_ID, one_of(model_nodes_HUGO))

# flog.threshold(ERROR)
# 
# grid.draw(venn.diagram(list(CNA=META_CNA$PATIENT_ID, RNA=META_RNA$PATIENT_ID, Mut=META_mut_patients, Clinical=META_clin$PATIENT_ID), filename=NULL, col = "transparent", fill = c("cornflowerblue", "green", "yellow", "darkorchid1"), alpha = 0.5, label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"), cex = 1.5, fontface = "bold", cat.default.pos = "text", cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5, cat.pos = 0, main = "METABRIC - Patients and available data types and clinical information"))
# 
# grid.newpage()
# grid.draw(venn.diagram(list(Model=model_nodes_HUGO, Mutants = unique(META_mut$Hugo_Symbol), CNA=META_CNA %>% select(-PATIENT_ID) %>% colnames, RNA=META_RNA %>% select(-PATIENT_ID) %>% colnames), filename=NULL, col = "transparent", fill = c("cornflowerblue", "green", "yellow", "darkorchid1"), alpha = 0.5, label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"), cex = 1.5, fontface = "bold", cat.default.pos = "text", cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5, cat.pos = 0, main = "Model genes and available data types"))

ggplot(META_clin)+geom_bar(aes(x=PAM50, fill=PAM50)) + ggtitle("Distribution of PAM50 subtypes in META cohort")

pca_data <- META_RNA %>% select(intersect(PAM50$Gene,colnames(.))) %>% prcomp %>% .$x %>% as.tibble %>%
  mutate(PATIENT_ID=substr(META_RNA$PATIENT_ID,1,12)) %>% full_join(META_clin, by="PATIENT_ID")

ggplot(pca_data, aes(x=PC1, y=PC2, colour=PAM50))+geom_point()+ggtitle("PAM50-Subtype and PAM50-genes PCA projection (based on RNAseq)")

#function to compute the Bimodality Index (BI) described in Wang et al. (2009)
BI <- function(dataset) {
  x <- dataset
  mc <- Mclust(na.omit(x), G = 2, modelNames = "E", verbose = FALSE)
  if (is.null(mc)) {
    bi <- NA
  } else {
    sigma <- sqrt(mc$parameters$variance$sigmasq)
    delta <- abs(diff(mc$parameters$mean))/sigma
    pi <- mc$parameters$pro[1]
    bi <- delta * sqrt(pi*(1-pi))
  }
  bi
}


#function to binarize the tails of the distribution, based on inter-quartile range (IQR), similar to methods described in teh outlier-sum statistic (Tibshirani and Hastie, 2007). Can be called with a reference dataset
OSclass <- function(exp_dataset, ref_dataset=exp_dataset) {
  classif <-rep(NA,length(exp_dataset))
  q25 <- quantile(ref_dataset,0.25, na.rm = T)
  q75 <- quantile(ref_dataset,0.75, na.rm = T)
  IQR <- q75 - q25 #InterQuartile Range
  classif[exp_dataset>IQR+q75] <- 1
  classif[exp_dataset<q25-IQR] <- 0
  return(classif)
}

#function to to binarize bimodal distributions based on a 2-modes gaussian mixture model (with equal variances). Can be called with a reference dataset
BIMclass <- function(exp_dataset, ref_dataset=exp_dataset) {
  mc <- Mclust(na.omit(ref_dataset), modelNames = "E", G=2, verbose = FALSE)
  classif <- rep(NA,length(exp_dataset))
  if (diff(mc$parameters$mean)>0){
    thresh_down <- max(mc$data[mc$classification==1 & mc$uncertainty <= 0.05])
    thresh_up <- min(mc$data[mc$classification==2 & mc$uncertainty <= 0.05])
    classif[exp_dataset<=thresh_down] <- 0
    classif[exp_dataset>=thresh_up] <- 1
  } else if (diff(mc$parameters$mean)<0){
    thresh_down <- max(mc$data[mc$classification==2 & mc$uncertainty <= 0.05])
    thresh_up <- min(mc$data[mc$classification==1 & mc$uncertainty <= 0.05])
    classif[exp_dataset<=thresh_down] <- 0
    classif[exp_dataset>=thresh_up] <- 1
  }
  return(classif)
}

#function for normalization of zero-inflated data
norm_fun_lin <- function(xdat, reference = xdat){
  x_proc <- (xdat-quantile(reference, 0.01, na.rm = T))/quantile(xdat-quantile(reference, 0.01, na.rm = T), 0.99, na.rm = T)
  x_proc[x_proc<0] <- 0
  x_proc[x_proc>1] <- 1
  x_proc
}

#function for normalization of unimodal data
norm_fun_sig <- function(xdat, reference = xdat){
  xdat <- xdat - median(reference, na.rm = T)
  lambda <- log(3)/mad(reference, na.rm = T)
  transformation <- function(x){
    y <- 1/(1+exp(-lambda*x))
    y
  }
  transformation(xdat) 
}

#function for normalization of unimodal data
norm_fun_bim <- function(xdat, reference = xdat) {
  not_na_xdat <- !is.na(xdat)
  not_na_ref <- !is.na(reference)
  mc <- Mclust(reference[not_na_ref], modelNames = "E", G=2, verbose = FALSE)
  pred <- predict.Mclust(mc,xdat[not_na_xdat])
  normalization <- rep(NA,length(xdat))
  if (diff(mc$parameters$mean)>0){
    normalization[not_na_xdat] <- pred$z[,2]
  } else if (diff(mc$parameters$mean)<0){
    normalization[not_na_xdat] <- pred$z[,1]
  }
  normalization
}

#Here we compute all statistical tools and criteria needed to perform the classification of distributions in the following categories: discarded, zero-inflated, unimodal and bimodal
compute_criteria <- function(exp_dataset){
  exp_dataset <- exp_dataset %>% select(-PATIENT_ID)
  criteria <- tibble(Gene=colnames(exp_dataset), Dip=NA, BI=NA, Kurtosis=NA, DropOutRate=NA, MeanNZ=NA, DenPeak=NA, Amplitude=NA)
  
  #Compute
  pb = txtProgressBar(min = 1, max = ncol(exp_dataset), initial = 1) 
  for (i in 1:ncol(exp_dataset)){
    x <- na.omit(unlist(exp_dataset[,i]))
    criteria$Amplitude[i] <- max(x)-min(x)
    
    if (criteria$Amplitude[i] !=0){
      criteria$Dip[i] <- dip.test(x)$p.value
      criteria$BI[i] <- BI(x)
      criteria$Kurtosis[i] <- kurtosis(x)-3
      criteria$DropOutRate[i] <- sum(x==0)/length(x)
      criteria$MeanNZ[i] <- sum(x)/sum(x!=0)
      den <- density(x, na.rm = T)
      criteria$DenPeak[i] <- den$x[which.max(den$y)]
    }
    
    setTxtProgressBar(pb,i)
  }
  
  threshold <- median(criteria$Amplitude)/10
  criteria <- criteria %>% 
    mutate(Category=ifelse(Amplitude<threshold | DropOutRate>0.95, "Discarded", NA)) %>%
    mutate(Category=ifelse(is.na(Category) & (BI>1.5 & Dip<0.05 & Kurtosis < 1),"Bimodal",Category)) %>%
    mutate(Category=ifelse(is.na(Category) & DenPeak<threshold, "ZeroInf", Category)) %>%
    mutate(Category=ifelse(is.na(Category), "Unimodal", Category))
  
  return(criteria)
}

criteria_META <- compute_criteria(META_RNA)

# print("META assignments:")
# kable(t(table(criteria_META$Category)))
# 
# print("META assignments for model-related nodes:")
# kable(t(table(criteria_META %>% filter(Gene %in% model_nodes_HUGO) %>% select(Category))))

criteria_META %>% filter(Category=="Bimodal") %>% select(Gene) %>% slice(sample(nrow(.),5)) %>% unlist %>% select(META_RNA,.) %>%  gather %>%
  ggplot(mapping = aes(x = value)) + geom_histogram(bins = 30) + facet_wrap(~key, scales = 'free') + ggtitle("20 random Bimodal genes")

criteria_META %>% filter(Category=="Unimodal") %>% select(Gene) %>% slice(sample(nrow(.),20)) %>% unlist %>% select(META_RNA,.) %>%  gather %>%
  ggplot(mapping = aes(x = value)) + geom_histogram(bins = 30) + facet_wrap(~key, scales = 'free') + ggtitle("20 random Unimodal genes")


binarize_exp <-  function(exp_dataset, ref_dataset, ref_criteria, gene, show.plot=F){
  if(!missing(gene)){
    
    gene_cat <- ref_criteria %>% filter(Gene==gene) %>% select(Category) %>% unlist
    x <- unlist(select(exp_dataset,gene))
    x_ref <- unlist(select(ref_dataset,gene))
    
    if (gene_cat=="Discarded"){
      stop("Discarded gene")
      
    } else if (gene_cat=="Bimodal"){
      gene_bin <- BIMclass(x,x_ref)
      
    } else {
      gene_bin <- OSclass(x,x_ref)
    }
    names(gene_bin) <- exp_dataset$PATIENT_ID
    if(show.plot==T){
      if(all(is.na(gene_bin))){
        tibble(Continuous=x) %>% ggplot(aes(x=Continuous))+geom_histogram(bins=30)+ggtitle(gene)
      } else {
        tibble(Continuous=x, Discrete=factor(gene_bin)) %>% ggplot(aes(x=Continuous, fill=Discrete))+geom_histogram(bins=30)+ggtitle(gene)
      }
    } else {
      return(gene_bin)
    }
    
  } else {
    exp_dataset <- tbl_to_df(exp_dataset) 
    ref_dataset <- tbl_to_df(ref_dataset)
    if(dim(exp_dataset)[2] != dim(ref_criteria)[1]){stop("Different number of genes")}
    logi_dis <- ref_criteria$Category=="Discarded"
    logi_OS <- ref_criteria$Category=="Unimodal" | ref_criteria$Category=="ZeroInf"
    logi_bim <- ref_criteria$Category=="Bimodal"
    exp_dataset[,logi_dis] <- lapply(exp_dataset[,logi_dis], function(x) rep(NA, length(x)))
    exp_dataset[,logi_OS] <- mapply(function(x,y) OSclass(x,y), as.data.frame(exp_dataset[,logi_OS]), as.data.frame(ref_dataset[,logi_OS]))
    exp_dataset[,logi_bim] <- mapply(function(x,y) BIMclass(x,y), as.data.frame(exp_dataset[,logi_bim]), as.data.frame(ref_dataset[,logi_bim]))
    
    return(exp_dataset)
  }
  
}

print("Bimodal example:")
binarize_exp(META_RNA,META_RNA, criteria_META, "RPL14", T)
print("Unimodal example:")
binarize_exp(META_RNA,META_RNA, criteria_META, "ERBB2", T)


META_RNA_prof <- binarize_exp(META_RNA,META_RNA, criteria_META)

META_RNA_prof %>% select(one_of(model_nodes_HUGO)) %>% is.na %>% `!` %>% rowSums %>% sort(decreasing=TRUE) %>% plot(ylab="# of bin. RNA", xlab="Patients", main="Distribution of binarized RNA per patients in META cohort")

META_RNA_prof %>% select(one_of(model_nodes_HUGO)) %>% is.na %>% `!` %>% colSums %>% sort(decreasing=TRUE) %>%  head(n=50) %>% barplot(las=2, ylab="# of mutations across cohort", main="Distribution of mutations for top 50 genes in META cohort")

#Profiles
METAmodel_RNA_prof <- META_RNA_prof %>% select(one_of(model_nodes_HUGO))


normalize_exp <-  function(exp_dataset, ref_dataset, ref_criteria, gene, show.plot=F){
  if(!missing(gene)){
    
    gene_cat <- ref_criteria %>% filter(Gene==gene) %>% select(Category) %>% unlist
    x <- unlist(select(exp_dataset,gene))
    x_ref <- unlist(select(ref_dataset,gene))
    
    if (gene_cat=="Discarded"){
      stop("Discarded gene")
      
    } else if (gene_cat=="Bimodal"){
      gene_bin <- norm_fun_bim(x,x_ref)
      
    } else if (gene_cat=="Unimodal"){
      gene_bin <- norm_fun_sig(x,x_ref)
      
    } else {
      gene_bin <- norm_fun_lin(x,x_ref)
    }
    names(gene_bin) <- exp_dataset$PATIENT_ID
    
    if(show.plot==T){
      gene_bin %>% unlist %>% as.data.frame %>% ggplot(aes(x=.)) + geom_histogram(bins=30)+xlab(gene)
    } else {
      return(gene_bin)
    }
    
  } else {
    exp_dataset <- tbl_to_df(exp_dataset) 
    ref_dataset <- tbl_to_df(ref_dataset)
    if(dim(exp_dataset)[2] != dim(ref_criteria)[1]){stop("Different number of genes")}
    logi_dis <- ref_criteria$Category=="Discarded"
    logi_uni <- ref_criteria$Category=="Unimodal"
    logi_zero <- ref_criteria$Category=="ZeroInf"
    logi_bim <- ref_criteria$Category=="Bimodal"
    exp_dataset[,logi_dis] <- lapply(exp_dataset[,logi_dis], function(x) rep(NA, length(x)))
    exp_dataset[,logi_uni] <- mapply(function(x,y) norm_fun_sig(x,y), as.data.frame(exp_dataset[,logi_uni]), as.data.frame(ref_dataset[,logi_uni]))
    exp_dataset[,logi_zero] <- mapply(function(x,y) norm_fun_lin(x,y), as.data.frame(exp_dataset[,logi_zero]), as.data.frame(ref_dataset[,logi_zero]))
    exp_dataset[,logi_bim] <- mapply(function(x,y) norm_fun_bim(x,y), as.data.frame(exp_dataset[,logi_bim]), as.data.frame(ref_dataset[,logi_bim]))
    
    return(exp_dataset)
  }
  
}

META_RNA_prof_norm <- normalize_exp(META_RNA, META_RNA, criteria_META)

METAmodel_RNA_prof_norm <- META_RNA_prof_norm %>% select(one_of(model_nodes_HUGO))

inverse_mapping <- function(dataset, datatype = 'bin', geneindex_var){
  
  geneindex_var %<>% mutate(V1=as.character(V1), V2=as.character(V2))
  inv_dataset <- data.frame(row.names = rownames(dataset))
  patients <- rownames(dataset)
  nodes_list <- geneindex_var$V1 %>% unique %>% as.character
  
  for (node in nodes_list){
    mapped_genes <- geneindex_var$V2[geneindex_var$V1 == node]
    if (any(mapped_genes != "") & any(mapped_genes %in% colnames(dataset))){
      if (length(mapped_genes) == 1){
        new_column <- select(dataset, one_of(mapped_genes))
      } else if (datatype == 'norm'){
        new_column <-  data.frame(rowMeans(select(dataset, one_of(mapped_genes)), na.rm = T))
      } else if (datatype == 'bin'){
        interm <- rowMeans(select(dataset, one_of(mapped_genes)), na.rm = T)
        interm[interm!=0 & interm !=1] <- NA
        interm[is.nan(interm)] <- NA
        new_column <- data.frame(interm)
      }
      inv_dataset[patients,node] <- new_column[patients,]
    }
  }
  return(inv_dataset)
}

if (!dir.exists("PROFILE-master/Results/Profiles")){
  dir.create("PROFILE-master/Results/Profiles")
}


write.csv(inverse_mapping(METAmodel_RNA_prof, 
                          geneindex_var = geneindex),
          "PROFILE-master/Results/Profiles/Visium_META_RNA.csv")

write.csv(inverse_mapping(METAmodel_RNA_prof_norm,
                          datatype='norm', 
                          geneindex_var = geneindex),
          "PROFILE-master/Results/Profiles/Visium_META_RNA_norm.csv")
























