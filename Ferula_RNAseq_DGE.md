# Ferula_RNAseq_DEG
Hajar  
June 8, 2017  


#import data


```r
counts <- read.csv("~/Documents/Ferula_RNAseq_Rbased/input/combined_counts_V2.csv", header=T, row.names="target_id")
#head(counts) for pooling two library of BR3 (they are not replication)
counts$BR3.2 <- counts$BR3 + counts$BR3.1
counts$BR3 <- NULL
counts$BR3.1 <- NULL
colnames(counts)[19] <- "BR3"
counts$BF3.2 <- counts$BF3 + counts$BF3.1
counts$BF3 <- NULL
counts$BF3.1 <- NULL
colnames(counts)[18] <- "BF3"
#colnames(counts)[1] 
#rownames(counts)
dim(counts) # 158617 18 158617 is num of transcripts in Reference
```

```
## [1] 158617     18
```

```r
write.csv(counts, file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/input/Ferula_RNAseq.coutns_v2.csv")
```
#filter based on read count, assign group, normalize, design matrix

```r
hist(colSums(counts,na.rm=TRUE))
```

![](Ferula_RNAseq_DGE_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
colSums(counts,na.rm=TRUE)
```

```
##      DS6      DF6      DL6      DR3      DS3      DF3      DL3      DR2 
## 22185775 22738489 20446520 21387648 19793593  9549581 22938854 17169702 
##      DS2      DF2      DL2      DR6      NF3      NF6      NR3      NR6 
## 26336386 21407169 14960410 21625910 20917468 17632293 20384896 16792343 
##      BR3      BF3 
## 14958214 14782562
```

```r
#general threshold
colSums(counts,na.rm=TRUE) > 1000000 # all samples are true
```

```
##  DS6  DF6  DL6  DR3  DS3  DF3  DL3  DR2  DS2  DF2  DL2  DR6  NF3  NF6  NR3 
## TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE 
##  NR6  BR3  BF3 
## TRUE TRUE TRUE
```

```r
counts.nolow <- counts[,colSums(counts,na.rm=TRUE) > 1000000]
dim(counts.nolow) #all of samples has more than 1000000 counts
```

```
## [1] 158617     18
```

```r
#sample description
samples <- data.frame(file=colnames(counts),
                      facility=factor(sub("(B|D|N)(S|F|L|R)(2|3|6)","\\1",colnames(counts))),
                      trt=factor(sub("(B|D|N)(S|F|L|R)(2|3|6)","\\2",colnames(counts))),
                    
                      genotype=factor(sub("(B|D|N)(S|F|L|R)(2|3|6)","\\3",colnames(counts)))) 
head(samples) 
```

```
##   file facility trt genotype
## 1  DS6        D   S        6
## 2  DF6        D   F        6
## 3  DL6        D   L        6
## 4  DR3        D   R        3
## 5  DS3        D   S        3
## 6  DF3        D   F        3
```

```r
#convert NA to zero
counts[is.na(counts)]<-0
# eliminating genes with low expression levels by retaining genes with > 10 reads in >= 3 samples
counts.small <-counts[rowSums(counts > 10) >= 3,] 
dim(counts.small) #61234 18 
```

```
## [1] 62147    18
```

```r
dim(counts) #158617 18
```

```
## [1] 158617     18
```

```r
write.csv(counts.small, file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/input/Ferula_RNAseq.coutns.small_v2.csv") # use this for other analysis
samples.small <- data.frame(file=colnames(counts.small),                     facility=factor(sub("(B|D|N)(S|F|L|R)(2|3|6)","\\1",colnames(counts.small))),
                      trt=factor(sub("(B|D|N)(S|F|L|R)(2|3|6)","\\2",colnames(counts.small))),
                    
                      genotype=factor(sub("(B|D|N)(S|F|L|R)(2|3|6)","\\3",colnames(counts.small)))) 
head(samples.small)
```

```
##   file facility trt genotype
## 1  DS6        D   S        6
## 2  DF6        D   F        6
## 3  DL6        D   L        6
## 4  DR3        D   R        3
## 5  DS3        D   S        3
## 6  DF3        D   F        3
```

```r
save(samples.small,file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/input/Ferula_RNAseq.samples.small_v2.Rdata")
#assign group by combining all the experimental factors into one combined factor
Genotype<-levels(samples.small$genotype)
samples.small$group <- paste(samples.small$genotype,samples.small$trt,samples.small$facility,sep=".")
samples.small$genotype<-as.character(samples.small$genotype)
```

#install.packages("edgeR")

```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
dge <- DGEList(counts=counts.small, group=samples.small$group) 
length(colnames(counts.small)) # 20
```

```
## [1] 18
```

```r
dge<-calcNormFactors(dge, method = "TMM")
# look at the normalization factors
nrow(dge$samples) # 18 
```

```
## [1] 18
```

```r
hist(dge$samples[,3]) # norm.factor > 3 seem outlier? I don't want to remove sample 
```

![](Ferula_RNAseq_DGE_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
plot(log10(dge$sample[,"lib.size"]),dge$sample[,"norm.factors"]) 
```

![](Ferula_RNAseq_DGE_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
#when we want to change the ref from first icon to others
samples.small$genotype <- as.factor(samples.small$genotype)
samples.small$genotype <- relevel(samples.small$genotype,ref="3")
samples.small$trt <- as.factor(samples.small$trt)
samples.small$trt <- relevel(samples.small$trt,ref="R")

#design model for each factor (genotype & trt factors)
design1 <- model.matrix(~genotype+trt + facility, data=samples.small)
colnames(design1)
```

```
## [1] "(Intercept)" "genotype2"   "genotype6"   "trtF"        "trtL"       
## [6] "trtS"        "facilityD"   "facilityN"
```

```r
#First the overall dispersion
dge <- estimateGLMCommonDisp(dge,design1, verbose = T) #Disp = 1.15 , BCV =  1.076 
```

```
## Disp = 0.89357 , BCV = 0.9453
```

```r
dge <- estimateGLMTrendedDisp(dge,design1)
dge <- estimateGLMTagwiseDisp(dge,design1)
save(dge,file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/dge_v2.Rdata")
plotBCV(dge)
```

![](Ferula_RNAseq_DGE_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
dev.copy(png,'BCV.png')
```

```
## quartz_off_screen 
##                 3
```

```r
mds.dge <- plotMDS(dge, method = "bcv",labels = dge$samples$group)
dge.fit <- glmFit(dge, design1)
#colnames(dge.fit)
```

#To find genes that are differentially expressed in gt 2 & 6 vs 3


```r
dge.lrt <- glmLRT(dge.fit,coef = c("genotype2","genotype6"))
#the top 10 most differentially expressed genes
#topTags(dge.lrt)
#summary(decideTestsDGE(dge.lrt,p=0.05)) why error?
#Extract genes with a FDR < 0.01 (could also use 0.05)
DEgenes1 <- topTags(dge.lrt,n = Inf)$table[topTags(dge.lrt,n = Inf)$table$FDR<0.05,]
dim(DEgenes1) #2096   6
```

```
## [1] 2700    6
```

```r
colnames(DEgenes1)
```

```
## [1] "logFC.genotype2" "logFC.genotype6" "logCPM"          "LR"             
## [5] "PValue"          "FDR"
```

```r
#head(DEgenes1)
#save to a file
write.csv(DEgenes1,file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/Ferula_RNAseq.DEgenes1_v2.csv")
#To find genes that are differentially expressed in gt 2 & 6 vs 3 seperately
dge.lrt.gt2 <- glmLRT(dge.fit,coef = c("genotype2"))
#topTags(dge.lrt.gt2)
summary(decideTestsDGE(dge.lrt.gt2, p=0.05))
```

```
##    [,1] 
## -1   455
## 0  61406
## 1    286
```

```r
dge.lrt.gt6 <- glmLRT(dge.fit,coef = c("genotype6"))
#topTags(dge.lrt.gt6)
summary(decideTestsDGE(dge.lrt.gt6, p=0.05))
```

```
##    [,1] 
## -1   903
## 0  60640
## 1    604
```

```r
DEGs1 <- data.frame ("gt2" = as.data.frame(summary(decideTestsDGE(dge.lrt.gt2, p=0.05)))$Freq,
                                     "gt6" = as.data.frame(summary(decideTestsDGE(dge.lrt.gt6, p=0.05)))$Freq)
rownames(DEGs1) <- c("down", "no", "up")
DEGs1 <- DEGs1[c("down", "up"),]
DEGs1
```

```
##      gt2 gt6
## down 455 903
## up   286 604
```

```r
save(DEGs1, file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/Ferula_RNAseq.DEGs1_v2.Rdata")
library(reshape2)
```

```
## Warning: package 'reshape2' was built under R version 3.2.5
```

```r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 3.2.5
```

```r
DEGs1.melt <- melt(DEGs1)
```

```
## No id variables; using all as measure variables
```

```r
DEGs1.melt$DE <- rep(c("down", "up"), 2)
colnames(DEGs1.melt) <- c("genotype", "number", "DE")
DEGs1.melt
```

```
##   genotype number   DE
## 1      gt2    455 down
## 2      gt2    286   up
## 3      gt6    903 down
## 4      gt6    604   up
```

```r
# reorder: up 1st down 2nd 
DEGs1.melt$DE <- factor(DEGs1.melt$DE, levels = c("up", "down"))
DEGs1.melt <- DEGs1.melt[order(DEGs1.melt$DE),]
DEGs1.melt
```

```
##   genotype number   DE
## 2      gt2    286   up
## 4      gt6    604   up
## 1      gt2    455 down
## 3      gt6    903 down
```

```r
DEGs1.melt$gt <- gsub("(X)(\\.)(S|L|F)(\\.)", "\\1",DEGs1.melt$genotype)

DEGs1.melt$trt <- gsub("(X)(\\.)(R|S|L|F)(\\.)", "\\2",DEGs1.melt$genotype)


### making ggplot for DEGs
library(ggplot2)
p.DEGs1 <- ggplot(data = DEGs1.melt)
p.DEGs1 <- p.DEGs1 + geom_bar(mapping = aes(fill=DE, x = factor(DE), y = number) , stat= "identity")
#p.DEGs1 <- p.DEGs1 + facet_grid(~genotype) 
p.DEGs1 <- p.DEGs1 + labs(y = "number of differentially expressed genes", x = "")
p.DEGs1
```

![](Ferula_RNAseq_DGE_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggsave(p.DEGs1, file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/Ferula_RNAseq.p.DEG1_v2.png")
```

```
## Saving 7 x 5 in image
```

#To find genes that are differentially expressed in trt S, F & L vs R 


```r
dge.lrt.trt <- glmLRT(dge.fit,coef = c("trtF","trtL", "trtS"))
#the top 10 most differentially expressed genes
#topTags(dge.lrt.trt)
#summary(decideTestsDGE(dge.lrt,p=0.05)) why error?
#Extract genes with a FDR < 0.01 (could also use 0.05)
DEgenes2 <- topTags(dge.lrt.trt,n = Inf)$table[topTags(dge.lrt.trt,n = Inf)$table$FDR<0.05,]
write.csv(DEgenes2,file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/Ferula_RNAseq.DEgenes1_v2.csv")
dim(DEgenes2) #20140 7
```

```
## [1] 12687     7
```

```r
#colnames(DEgenes2)
head(DEgenes2)
```

```
##                          logFC.trtF logFC.trtL logFC.trtS   logCPM
## TRINITY_DN59790_c0_g1_i1   8.990086  1.5099392  1.8899683 2.964240
## TRINITY_DN15353_c0_g2_i1  10.000930  1.0079972  0.4142302 3.419545
## TRINITY_DN67795_c0_g1_i1  10.819517  5.5040013  0.7663922 3.058659
## TRINITY_DN21262_c0_g1_i2   8.891950 -2.5158477 -3.2335133 4.173338
## TRINITY_DN13576_c0_g1_i1   9.763360  3.8651982  1.7732271 4.454639
## TRINITY_DN32850_c0_g3_i1   9.708344 -0.7907408 -1.3949423 1.793940
##                                LR       PValue          FDR
## TRINITY_DN59790_c0_g1_i1 147.6280 8.558943e-32 5.319126e-27
## TRINITY_DN15353_c0_g2_i1 145.8443 2.075646e-31 6.449757e-27
## TRINITY_DN67795_c0_g1_i1 140.4021 3.095758e-30 5.463040e-26
## TRINITY_DN21262_c0_g1_i2 140.1456 3.516205e-30 5.463040e-26
## TRINITY_DN13576_c0_g1_i1 132.9266 1.265682e-28 1.573167e-24
## TRINITY_DN32850_c0_g3_i1 124.2949 9.168993e-27 9.497090e-23
```

```r
#create .txt file for creating gff3
names <- rownames(DEgenes2[1:1000,])
names <- paste0("<",names)
names <- data.frame(names)
write.table(names, "top_transcripts.txt", col.names = F, row.names = F, quote = F)

#for all DE genes
names2 <-rownames(DEgenes2)
names2 <- paste0("<",names2)
names2 <- data.frame(names2)
write.table(names, "All_transcripts.txt", col.names = F, row.names = F, quote = F)
#To find genes that are differentially expressed in trt S, F & L vs R seperately
dge.lrt.trtF <- glmLRT(dge.fit,coef = c("trtF"))
#topTags(dge.lrt.trtF)
summary(decideTestsDGE(dge.lrt.trtF, p=0.05))
```

```
##    [,1] 
## -1  2271
## 0  54876
## 1   5000
```

```r
dge.lrt.trtL <- glmLRT(dge.fit,coef = c("trtL"))
#topTags(dge.lrt.trtL)
summary(decideTestsDGE(dge.lrt.trtL, p=0.05))
```

```
##    [,1] 
## -1  3872
## 0  55421
## 1   2854
```

```r
dge.lrt.trtS <- glmLRT(dge.fit,coef = c("trtS"))
#topTags(dge.lrt.trtS)
summary(decideTestsDGE(dge.lrt.trtS, p=0.05))
```

```
##    [,1] 
## -1   180
## 0  61558
## 1    409
```

```r
DEGs2 <- data.frame ("trtF" = as.data.frame(summary(decideTestsDGE(dge.lrt.trtF, p=0.05)))$Freq,
                                     "trtL" = as.data.frame(summary(decideTestsDGE(dge.lrt.trtL, p=0.05)))$Freq,
                                     "trtS"= as.data.frame( summary(decideTestsDGE(dge.lrt.trtS, p=0.05)))$Freq)
rownames(DEGs2) <- c("down", "no", "up")
DEGs2 <- DEGs2[c("down", "up"),]
DEGs2
```

```
##      trtF trtL trtS
## down 2271 3872  180
## up   5000 2854  409
```

```r
save(DEGs2, file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/Ferula_RNAseq.DEGs2_v2.Rdata")
#library(reshape2)
#library(ggplot2)
DEGs2.melt <- melt(DEGs2)
```

```
## No id variables; using all as measure variables
```

```r
DEGs2.melt$DE <- rep(c("down", "up"), 3)
colnames(DEGs2.melt) <- c("genotype", "number", "DE")
DEGs2.melt
```

```
##   genotype number   DE
## 1     trtF   2271 down
## 2     trtF   5000   up
## 3     trtL   3872 down
## 4     trtL   2854   up
## 5     trtS    180 down
## 6     trtS    409   up
```

```r
# reorder: up 1st down 2nd 
DEGs2.melt$DE <- factor(DEGs2.melt$DE, levels = c("up", "down"))
DEGs2.melt <- DEGs2.melt[order(DEGs2.melt$DE),]
DEGs2.melt
```

```
##   genotype number   DE
## 2     trtF   5000   up
## 4     trtL   2854   up
## 6     trtS    409   up
## 1     trtF   2271 down
## 3     trtL   3872 down
## 5     trtS    180 down
```

```r
DEGs2.melt$gt <- gsub("(X)(\\.)(S|L|F)(\\.)", "\\1",DEGs2.melt$genotype)

DEGs2.melt$trt <- gsub("(X)(\\.)(R|S|L|F)(\\.)", "\\2",DEGs2.melt$genotype)
### making ggplot for DEGs
library(ggplot2)
p.DEGs2 <- ggplot(data = DEGs2.melt)
p.DEGs2 <- p.DEGs2 + geom_bar(mapping = aes(fill=DE, x = factor(DE), y = number) , stat= "identity")
p.DEGs2 <- p.DEGs2 + facet_grid(~genotype) 
p.DEGs2 <- p.DEGs2 + labs(y = "number of differentially expressed genes", x = "")
p.DEGs2
```

![](Ferula_RNAseq_DGE_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
ggsave(p.DEGs2, file="/Users/hajaramini/Documents/Ferula_RNAseq_Rbased/output/Ferula_RNAseq.p.DEG2_v2.png")
```

```
## Saving 7 x 5 in image
```

