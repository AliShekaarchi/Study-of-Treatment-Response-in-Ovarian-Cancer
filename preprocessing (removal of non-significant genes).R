
############################################################################################
############################################################################################
############################################################################################

library(affy)

setwd("C:/My Files/BIO project/GSE30161_RAW")


### Data Normalizing
### Method 1 


# Data Loading


my_Cel_files <- ReadAffy()
data <- exprs(my_Cel_files)

head(data)
View(data)


hist(data , breaks = 20000 , xlim = c(0 , 200))

hist(data , breaks = 10000 , xlim = c(0 , 100))
hist(data[,2])
hist(data[,1])

boxplot(my_Cel_files)



hist(data[,1], breaks = 10000,xlim = c(0,200))

data_filtered <- data[which(rowSums(data)<1000),]
data_filtered <- data[which(rowMeans(data) < 100),]

nrow(data_filtered)
nrow(data)

hist(data_filtered, breaks = 10,xlim = c(0,2000) )

sum(data == 0)
sum(data < 40)

nrow(data)


# RMA Normalization


norm_data_with_rma <- affy::rma(my_Cel_files) # In order not to be confused with the Oligo package

rma_dataframe <- exprs(norm_data_with_rma)
View(rma_dataframe)
head(rma_dataframe)
max(rma_dataframe)

sum(rma_dataframe == 0)

hist(rma_dataframe)

hist(rma_dataframe[,1] , breaks = 100 , xlim = c(0,14))

boxplot(rma_dataframe)

png("RMA Normalized.png", width = 1920, height = 1080)
boxplot(rma_dataframe)
dev.off()

nrow(data)
nrow(rma_dataframe)

write.table(rma_dataframe , "RMA Normalized data.txt" , quote = F , sep = "\t")

# mas5 Normalization

norm_data_with_mas5 <- mas5(my_Cel_files)
mas5_dataframe <- exprs(norm_data_with_mas5)

View(mas5_dataframe)
head(mas5_dataframe)

mas5_dataframe_log <- log2(mas5_dataframe)

View(mas5_dataframe_log)

boxplot(mas5_dataframe_log)

png("mas5 Normalized.png", width = 1920, height = 1080)
boxplot(mas5_dataframe_log)
dev.off()

# RMA with mas5 comparison

pdf("RMA vs mas5 .pdf", width = 12, height = 8)

boxplot(rma_dataframe)
boxplot(mas5_dataframe_log)

dev.off()

# Other Types

library(tkWidgets)

normal_data <- expresso(my_Cel_files , widget = T)
n_data <- exprs(normal_data)

### Annotation :

View(rma_dataframe)

annotation_file <- read.csv("GPL570-55999.csv")
View(annotation_file)

nrow(annotation_file)
nrow(rma_dataframe)


final_data <- rma_dataframe
View(final_data)

n <- 54675

rownames(final_data)[n] == annotation_file[n,1]

symbol <- annotation_file[1:54675,11]

#null_symb <- rep("Null" , 195)

#symbol <- append(symbol,null_symb)

length(symbol)
nrow(final_data)

rownames(final_data) <- symbol
View(final_data)


write.table(final_data , "final data.txt" , quote = F , sep = "\t")

### Annotation : 


View(annotation_file)

annot_df <- c()
annot_df <- cbind(annot_df ,annotation_file$Gene.Symbol)
annot_df <- cbind(annot_df ,annotation_file$ENTREZ_GENE_ID)
rownames(annot_df) <- annotation_file$ID

colnames(annot_df) <- c("Symbol" , "Entrez")

View(annot_df)
annot_df <- as.data.frame(annot_df)


annot_df["1053_at",1]


Prime.View <- function(x){
  
 return(annot_df[x,1])
  
}

Prime.View("1053_at")

View(rma_dataframe)
rownames(rma_dataframe)[1:10]

df <- rma_dataframe[200:300 , ]
View(df)

Prime.View(rownames(df))
rownames(df) <- Prime.View(rownames(df))

df_test <- final_data[200:300,]
View(df_test)
