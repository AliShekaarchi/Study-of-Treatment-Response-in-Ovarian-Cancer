
############################################################################################
############################################################################################
############################################################################################


setwd("C:/My Files/BIO project/GSE30161_RAW")



data <- read.table("final data.txt" , sep = "\t" , header = T, row.names = NULL)
View(data)

data_2 <- data[data$row.names!="Null",]
View(data_2)

data_3 <- data_2[!(data_2$row.names=="" | is.na(data_2$row.names)),]
View(data_3)

# Example #

data_exmp <- data.frame(
  character = c("A", "B", "A", "C", "B"),
  num1 = c(10, 20, 30, 40, 50),
  num2 = c(15, 25, 35, 45, 55)
)
View(data_exmp)


data_exmp_result<- aggregate(. ~ character, data = data_exmp, FUN = mean)
View(data_exmp_result)


######

data_4 <- aggregate(. ~ row.names, data = data, FUN = mean)
View(data_4)

data_5 <- subset(data_4, select=-c(GSM746878.CEL, GSM746882.CEL, GSM746917.CEL, GSM746890.CEL, GSM746896.CEL, GSM746898.CEL, GSM746885.CEL))



rownames(data_5) <- data_5[,1]
final_data <- data_5[,-1]
View(final_data)


write.csv(final_data , "final table.csv")
### Limma :

library(limma)

Groups <- Groups <- c(
  "PR (partial response)",
  "CR (completely response)",
  "PR (partial response)",
  "PR (partial response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "PR (partial response)",
  "PR (partial response)",
  "CR (completely response)",
  "PR (partial response)",
  "PR (partial response)",
  "CR (completely response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "PR (partial response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "PR (partial response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "CR (completely response)",
  "PR (partial response)"
)

Design <- model.matrix(~0+factor(Groups))  

View(Design)

dim(final_data)
dim(Design)
 
colnames(Design) <- c("CR" , "PR")

fit <- lmFit(final_data ,design =  Design)
View(fit)


cont <- makeContrasts(contrasts = "CR - PR" , levels = Design)
fit2 <- contrasts.fit(fit = fit , contrasts = cont)

fit3 <- eBayes(fit2 , 0.01)

Result <- topTable(fit3, adjust="fdr", sort.by="B", number=Inf)
View(Result)

Result <- Result[order(Result$adj.P.Val , decreasing = F) , ]
View(Result)

###

write.table(Result , "Final Result.txt" , quote = F , sep  = "\t")

###

Groups <- c("Normal" ,"Normal","Normal",
            "Treated" , "Treated" ,"Tumor" )

Design <- model.matrix(~0+factor(Groups))  
colnames(Design) <- c("Normal" , "Treated" ,"Tumor" )

fit <- lmFit(final_data ,design =  Design)
View(fit)

groups_cont <- make.names(c("Normal" , "Treated" ,"Tumor" ))
cts <- paste(groups_cont, c(tail(groups_cont, -1), head(groups_cont, 1)), sep="-")

cont <- makeContrasts(contrasts = cts  , levels = Design)

fit2 <- contrasts.fit(fit = fit , contrasts = cont)

fit3 <- eBayes(fit2 , 0.01)
?eBayes()

Result <- topTable(fit3, adjust="fdr", sort.by="B", number=Inf)
View(Result)

Result <- Result[order(Result$adj.P.Val , decreasing = F) ,]
View(Result)

Result["ETV7",]

###

write.table(Result , "Final Result.txt" , quote = F , sep  = "\t")
write.csv(Result , "Final Result.csv")

View(final_data)
View(Result)

rownames(Result)[1:100]

top_genes <- final_data[rownames(Result)[1:100],]
View(top_genes)

write.csv(top_genes , "Top genes.csv")
