#working with results of transdecoder to find best hit for 1000 DEgenes based on uniref90 .tab.gz file is the result of uniprof90)

library(tidyverse)

GO <- read_csv("GO_Term.txt", col_names = F)

Trans <- read_table2("transcript_names.txt", col_names = F)

Trans$X1 <- sub("::.*","",Trans$X1)

Trans$X2 <- sub("UniRef90_","",Trans$X2)



New <- inner_join(Trans,GO, by = c("X2" = "X1"))

New <- New[,-2]



write_csv(New, "Genes_and_Terms.tab" , col_names =F)

system("sed 's/NA//g' Genes_and_GO_Terms.tab |            
       
       sed 's/,/\t/' > Gene_and_GO_Terms.tab")



#062017 Working for creating S9 file like in BIS180L

library(tidyverse)

temp <- read_delim("uniref_GO_terms.tab", "\t", col_names = F)

temp$X2 <- gsub("\\[GO:[0-9]{7}\\]","",temp$X2)

temp$X2 <- gsub(";","/",temp$X2)

trans <- read_delim("transcript_names.txt", "\t", col_names = F)

trans$X1 <- sub("::.*","",trans$X1)

trans$X2 <- sub("UniRef90_","",trans$X2)

combined <- merge(trans,temp, by.x = "X2", by.y = "X1")

combined$X2 <- sub(" $","",combined$X2)

write.table(combined, "Genes_and_Terms.tab", sep = "\t", row.names = F, col.names = F, quote = F)

#distribution of contig length of trinity result

test1=read.table("Trinity.fasta.bed",header=T,row.names=1)
View(test1)
colnames(test1)[1] <- "target_id"
colnames(test1)[2] <- "length"
hist(test1)
hist(test1$length)
hist(test1$length,breaks=50)
test_short <- test1[test1$length <=1000, ]
dim(test_short)
test_short <- test1[test1$length <=400, ]
test_short <- test1[test1$length <=300, ]
hist(test_short$length,xlim=c(0,1000))
###add discription of BUSCO for the result of running BUSCO
info <- read_delim("embryophyta_3193_OrthoDB9_orthogroup_info.txt", delim = "\t", col_names = T)
complete <- read_delim("Complete_BUSCO_name.txt", col_names = F, delim = "\t")
info <- info[,1:3]
test <- left_join(complete, info, by = c("X1" = "OrthoGroupID"))