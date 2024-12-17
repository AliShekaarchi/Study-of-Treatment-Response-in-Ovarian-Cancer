
############################################################################################
############################################################################################
############################################################################################

library(GEOquery)

setwd("C:/My Files/BIO project/New folder")

#GSE28739

#RAW Data

get_my_GSE <- getGEOSuppFiles("GSE28739") #Download
list.files("GSE28739") #Get File Names from Downloaded GSE Folder

untar("GSE28739/GSE28739_RAW.tar" , exdir = "GSE28739/CEL") #Extract zip


############################################################################################

#MetaData


sample_names <- c("GSM7985880_con1.txt" , "GSM7985881_con2.txt",
                  "GSM7985882_con3.txt","GSM7985883_uVB1.txt",
                  "GSM7985884_uVB2.txt","GSM7985885_uVB3.txt")


dir <- dir("C:/My Files/BIO project/New folder")
sample_names <- dir

sample_type <- c("sensitive", "sensitive", "sensitive", "sensitive",
                 "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
                 "sensitive", "sensitive", "sensitive", "sensitive",
                 "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
                 "sensitive", "sensitive", "sensitive", "sensitive",
                 "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
                 "sensitive", "sensitive", "sensitive", "sensitive",
                 "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
                 "sensitive", "sensitive", "sensitive", "sensitive",
                 "resistant", "resistant", "resistant", "resistant", "resistant", "resistant")
MetaData <- data.frame(sample_names,sample_type)
View(MetaData)

colnames(MetaData) <- c("FileName" , "Type")
############################################################################################

library(limma)


data <- read.maimages(files = MetaData , source = "agilent", green.only = T)

?read.maimages

hist(data$E , breaks = 10000 , xlim = c(0,250), main = "Histogram of foreground signals light before normalization" )
hist(data$Eb , breaks = 10000 , xlim = c(0,100), main = "Histogram of background signals light before normalization" )

png("foreground Histogram.png", width = 1920, height = 1080)
hist(data$E , breaks = 10000 , xlim = c(0,250), main = "Histogram of foreground signals light before normalization" )
dev.off()
png("background Histogram.png", width = 1920, height = 1080)
hist(data$Eb , breaks = 10000 , xlim = c(0,100), main = "Histogram of background signals light before normalization" )
dev.off()

sum(data$E == 0)
sum(data$E < 10)
sum(data$E < 20)


length(data$E[,1])
sqrt(length(data$E[,1]))

length(data$E[,1]) /246

mat <- matrix(data$E[,1], nrow = 246, byrow = TRUE)
image(mat , col = c("black","white"))

library(pheatmap)

mycols <- c("black" ,"grey15","grey20",
            "grey25","grey30","grey35","grey40","grey45","grey50",
            "grey55","grey60","grey65","grey70","grey75","grey80",
            "grey85","grey90","grey95","grey100","white")
pheatmap(mat , cluster_rows = F , cluster_cols = F , color = mycols)
        

BiocManager::install("AgiMicroRna")

library(AgiMicroRna)

boxplotMicroRna(log2(data$E), maintitle = "Boxplot of foreground signals light before normalization") #foreground light
boxplotMicroRna(log2(data$Eb), maintitle = "Boxplot of background signals light before normalization") #background light

png("foreground boxplot before normalization.png", width = 1920, height = 1080)
boxplotMicroRna(log2(data$E), maintitle = "Boxplot of foreground signals light before normalization") #foreground light
dev.off()
png("background boxplot before normalization.png", width = 1920, height = 1080)
boxplotMicroRna(log2(data$Eb), maintitle = "Boxplot of background signals light before normalization") #background light
dev.off()


plotDensityMicroRna(log2(data$E), maintitle="Density Histogram before normalization")

png("Density Histogram before normalization.png", width = 1920, height = 1080)
plotDensityMicroRna(log2(data$E), maintitle="Density Histogram before normalization")
dev.off()


hierclusMicroRna(object = log2(data$E), methdis = "euclidean", methclu = "complete", sel = F)
title(main = "\n\nBefore Normalization")
      
png("Hierarchical Cluster MicroRna before normalization.png", width = 1920, height = 1080)
hierclusMicroRna(object = log2(data$E), methdis = "euclidean", methclu = "complete", sel = F)
title(main = "\n\nBefore Normalization")
dev.off()


data_corrected <- backgroundCorrect(RG = data , method = "normexp" , offset = 15)
View(data_corrected)

?backgroundCorrect

data_norm <- normalizeBetweenArrays(object = data_corrected , method = "quantile")
View(data_norm)
View(data_norm$genes)



boxplotMicroRna(data_norm$E, , maintitle = "Boxplot of signals light after background correction & normalization")

png("Corrected Normalized Data Boxplot.png", width = 1920, height = 1080)
boxplotMicroRna(data_norm$E, , maintitle = "Boxplot of signals light after background correction & normalization")
dev.off()

plotDensityMicroRna(data_norm$E, maintitle="Density Histogram after backgroung correction & normalization")

png("Normalized Data Density Histogram.png", width = 1920, height = 1080)
plotDensityMicroRna(data_norm$E, maintitle="Density Histogram after backgroung correction & normalization")
dev.off()

hierclusMicroRna(object = data_norm$E, methdis = "euclidean", methclu = "complete", sel = F)
title(main = "\n\n After background correction & Normalization")

png("Normalized Data Hierarchical Clusters MicroRna.png", width = 1920, height = 1080)
hierclusMicroRna(object = data_norm$E, methdis = "euclidean", methclu = "complete", sel = F)
title(main = "\n\n After background correction & Normalization")
dev.off()

data_norm_E <- data_norm$E
View(data_norm_E)

data_norm_E <- na.omit(data_norm_E)

hist(data_norm_E , breaks = 100 , xlim = c(0,10), , main = "Histogram of signals light after background correction & normalization")

png("Corrected & Normalized Data Histogram of signals light.png", width = 1920, height = 1080)
hist(data_norm_E , breaks = 100 , xlim = c(0,10), , main = "Histogram of signals light after background correction & normalization")
dev.off()

data_norm_average <- avereps(x = data_norm , ID = data_norm$genes$GeneName)   

?avereps

data.frame(data_norm_average)
View(data_norm_average)
View(data_norm_average$E)

write.table(data_norm_average , "Back ground corrected & normalized data.csv" , quote = F , sep = "\t")


#####

Groups <-  c("sensitive", "sensitive", "sensitive", "sensitive",
             "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
             "sensitive", "sensitive", "sensitive", "sensitive",
             "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
             "sensitive", "sensitive", "sensitive", "sensitive",
             "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
             "sensitive", "sensitive", "sensitive", "sensitive",
             "resistant", "resistant", "resistant", "resistant", "resistant", "resistant",
             "sensitive", "sensitive", "sensitive", "sensitive",
             "resistant", "resistant", "resistant", "resistant", "resistant", "resistant")

Design <- model.matrix(~0+factor(Groups))  

colnames(Design) <- c("sensitive" , "resistant")

fit <- lmFit(data_norm_average ,design =  Design)

cont <- makeContrasts(contrasts = "sensitive - resistant" , levels = Design)
fit2 <- contrasts.fit(fit = fit , contrasts = cont)

fit3 <- eBayes(fit2 , 0.01)

Result <- topTable(fit3, adjust="fdr", sort.by="B", number=Inf)
View(Result)

Result <- Result[order(Result$adj.P.Val , decreasing = F) , ]
View(Result)


write.table(Result , "Contrast sensitive- resistant Result.xlsx" , quote = F , sep = "\t")
write.xlsx(Result, file = "Contrast_sensitive_resistant_Result.xlsx", quote = FALSE)


