##fullfull is for testing

#Creates proper matrix before using topic modeling and performing great analysis

library(Amelia)
library(tidyverse)


Amp <- read.table(snakemake@input[["CpG"]], header=TRUE, sep="\t")
Coords <- read.table(snakemake@input[["coords"]])

mergedCol <- unite(Coords, comb, V2, V3, sep = "-", remove = TRUE, na.rm = FALSE)
mergedChrDone <- unite(mergedCol, Coords, V1, comb, sep = ":", remove = TRUE, na.rm = FALSE)

names(Amp)[names(Amp) == "X"] <- "probes"
Amp$probes <- as.character(Amp$probes)


names(mergedChrDone)[names(mergedChrDone) == "V4"] <- "probes"
mergedChrDone$probes <-as.character(mergedChrDone$probes)

merged <- left_join(mergedChrDone, Amp, by = c("probes"))

mergedSansprobes <- within(merged, rm("probes"))

mergedSansprobes <- mergedSansprobes[which(rowMeans(!is.na(mergedSansprobes)) > 0.5), ]

idvars = c('Coords')

### Follows is a loop for bounds to prevent AmeliaII to make some negative values in the dataframe
nc = ncol(mergedSansprobes)

loopyBound = rbind(c(2, 0, 1))

n=3
while (n <= nc)
{
loopyBound = rbind(loopyBound, c(n, 0, 1))
n=n+1
}

Amp2 <- amelia(mergedSansprobes, m = 1, idvars = idvars, bound = loopyBound)

imputedMatrix <- Amp2$imputations[[1]]

imputedMatrix2 <- imputedMatrix[,-1]
rownames(imputedMatrix2) <- imputedMatrix[,1]

write.csv(imputedMatrix2, file=snakemake@output[[1]])
#############PrepMeta###############

samplemeta <- read.table(snakemake@input[["premeta"]], header=TRUE, sep="\t")


###OTHERS THAN BRCA::: _Subtype_mRNA', 'Subtype_Selected' & 'Subtype_Immune_Model_Based
#metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_Selected, Subtype_Immune_Model_Based)


###BRCA_FULL
metaframe <- samplemeta %>% select(donor_id, PAM50, ER.Status)


#replace NA with string
metaframe.c <- metaframe %>% mutate_all(as.character)
metaframe.c <- metaframe.c %>% replace(is.na(.), "NA")

#replace numbers as such: 1=Negative, 2=Positive, and LuminalA and LuminalB to "Luminal A" and "Luminal B"
metaframe.c$ER.Status <- gsub("1","Negative", metaframe.c$ER.Status)
metaframe.c$ER.Status <- gsub("2","Positive", metaframe.c$ER.Status)
metaframe.c$PAM50 <- gsub("LuminalA","Luminal A", metaframe.c$PAM50)
metaframe.c$PAM50 <- gsub("LuminalB","Luminal B", metaframe.c$PAM50)

#If redundancy in data (E.g same donors multiple times)
#metaSD <- metaframe.c %>% distinct(icgc_donor_id, .keep_all = TRUE)
metaSD <- metaframe.c %>% distinct(donor_id, .keep_all = TRUE)

#set rownames
meta2 <- metaSD
rownames(meta2) <- metaSD[,1]


meta <- meta2 %>% select(-1)


###########RemoveNonMetaCpGRows############

drop <- c(rownames(meta))
CpG.sample.tab = imputedMatrix2[,(names(imputedMatrix2) %in% drop)]


############CISTOPIC########
suppressWarnings(library("cisTopic"))

#Name the PDF after what cancer type and TF:
PDF.name <- snakemake@output[[2]]

#Removing NAs by deletion: (CHECK IF SHOULD USE AMPUTATION!!!
#cpg <- CpG.sample.tab[complete.cases(CpG.sample.tab), ]
#Finding what to use from meta:
metacolnames <- colnames(meta)

cisTopicObject <- createcisTopicObject(is.acc=0.5, min.cells=0, min.regions=0, count.matrix=data.frame(CpG.sample.tab)) #Set is.acc=0 or is.acc=0.01 | min.cells=0, min.regions=0

cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(4:18), seed=123, nCores=4, addModels=FALSE)
#cisTopicObject <- runModels(cisTopicObject, topic=c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60), seed=123, nCores=20, addModels=FALSE)
#cisTopicObject <- runModels(cisTopicObject, topic=c(2,4,6,8), seed=123, nCores=20, addModels=FALSE)

pdf(PDF.name, height=10, width=15)
#cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject)
cisTopicObject <- runUmap(cisTopicObject, target ='cell')
cellTopicHeatmap(cisTopicObject, col.low="blue", col.mid="white",col.high="red", colorBy=c("PAM50","ER.Status"))
cellTopicHeatmap(cisTopicObject, method = "Probability",col.low="blue", col.mid="white",col.high="red", colorBy=c("PAM50","ER.Status"))

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Colored by metadata", pos=4)

#par(mfrow=c(1,2)) |changed intervals to 20
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status"), cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, cex.dot = 1.5)
###

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Topic by probability", pos=4)

#par(mfrow=c(3,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE,col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)
###
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Topic by Z-score", pos=4)

#par(mfrow=c(3,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE,col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Regionscores by Z-score", pos=4)


#############################################################
cisTopicObject <- getRegionsScores(cisTopicObject, method="Z-score", scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by probability", pos=4)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by Z-score", pos=4)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

#############################################################################
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Regionscore by NormTop", pos=4)

cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by probability", pos=4)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by Z-score", pos=4)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

###################################################
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"regionScore by probability", pos=4)

cisTopicObject <- getRegionsScores(cisTopicObject, method='Probability', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by probability", pos=4)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by Z-score", pos=4)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)





dev.off()

#### Write out topic assignments to the patients
write.csv(cisTopicObject@selected.model$document_expects, file=snakemake@output[[3]])

#### Region scores per topic (normalized umap)
write.csv(cisTopicObject@region.data, file=snakemake@output[[4]])

#### Unnormalized region assignments
write.csv(cisTopicObject@selected.model$topics, file=snakemake@output[[5]])

#cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=2)
#### UmapRegScore
#write.csv(cisTopicObject@dr$region, file=snakemake@output[[6]])

#### UmapPatScore
#write.csv(cisTopicObject@dr$cell, file=snakemake@output[[7]])


saveRDS(cisTopicObject, file=snakemake@output[[6]])
#getBedFiles(cisTopicObject, path='CisTopicBed_BRCA_SNA')


     snakemake@output[[7]])

write.csv(meta, file=snakemake@output[[8]])
