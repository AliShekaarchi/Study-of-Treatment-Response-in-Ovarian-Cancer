
############################################################################################
############################################################################################
############################################################################################

# Method 1 by Affy


install.packages("BiocManager")
library(BiocManager)

BiocManager::install("affy")
library(affy)


setwd("C:/My Files/BIO project/GSE30161")


my_Cel_files <- ReadAffy()
data <- exprs(my_Cel_files)
View(data)
head(data)
max(data)

write.table(data , "RAW data.txt" , quote = F , sep = "\t")


### Chip Image

image(my_Cel_files[2])

png("Selected Chip.png", width = 1920, height = 1080)
image(my_Cel_files[1])
dev.off()

### RNA Degradation


my_RNA_DEG <- AffyRNAdeg(my_Cel_files)
plotAffyRNAdeg(my_RNA_DEG)

png("RNA Degradation Plot.png", width = 1920, height = 1080)
plotAffyRNAdeg(my_RNA_DEG)
dev.off()



### Box Plot

boxplot(my_Cel_files)
boxplot(data)

View(data)


png("box plot before normaliziation.png", width = 1920, height = 1080)
boxplot(data)
dev.off()

png("box plot after normaliziation.png", width = 1920, height = 1080)
boxplot(my_Cel_files)
dev.off()



### Histogram

hist(my_Cel_files, main="Histogram of Samples", xlab="genes", ylab="Freq", col="blue")

png("Histogram.png", width = 1920, height = 1080)
hist(my_Cel_files, main="Histogram of Samples", xlab="genes", ylab="Freq", col="blue")
dev.off()


### PDF File


pdf("QC.pdf", width = 6, height = 4)
image(my_Cel_files[1])
plotAffyRNAdeg(my_RNA_DEG)
hist(my_Cel_files, main="Histogram of Samples", xlab="genes", ylab="Freq", col="blue")
boxplot(my_Cel_files)
dev.off()

### Loop for getting all chips photo


dir.create("chip images")
png_name <- c()

for (k in 1:ncol(data)){
  
  png_name <- paste0("chip images/" , "photo-" , k , ".png")
  png(png_name, width = 1920, height = 1080)
  image(my_Cel_files[k])
  dev.off()
  
  
}

############################################################################################
############################################################################################
############################################################################################

