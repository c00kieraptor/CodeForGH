#Creates proper matrix before using topic modeling and performing great analysis
library(Amelia)
library(tidyverse)


methTable <- read.table(snakemake@input[["CpG"]], header=TRUE, sep="\t")
Coords <- read.table(snakemake@input[["coords"]])
#Coords <- read.table("/storage/mathelierarea/processed/petear/SnakemakeInputFiles/BRCA-US_FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed")

mergedCol <- unite(Coords, comb, V2, V3, sep = "-", remove = TRUE, na.rm = FALSE)
mergedChrDone <- unite(mergedCol, Coords, V1, comb, sep = ":", remove = TRUE, na.rm = FALSE)
names(methTable)[names(methTable) == "X"] <- "probes"
methTable$probes <- as.character(methTable$probes)


names(mergedChrDone)[names(mergedChrDone) == "V4"] <- "probes"
mergedChrDone$probes <-as.character(mergedChrDone$probes)

methTable <- left_join(mergedChrDone, methTable, by = c("probes"))

methTable <- within(methTable, rm("probes"))

methTable <- methTable[which(rowMeans(!is.na(methTable)) > 0.5), ]


idvars = c('Coords')

### Follows is a loop for bounds to prevent AmeliaII to make some negative values in the dataframe
nc = ncol(methTable)

loopyBound = rbind(c(2, 0, 1))

n=3
while (n <= nc)
{
loopyBound = rbind(loopyBound, c(n, 0, 1))
n=n+1
}

methTable <- amelia(methTable, m = 1, idvars = idvars)

methTable <- methTable$imputations[[1]]

methTable2 <- methTable[,-1]
rownames(methTable2) <- methTable[,1]

write.csv(methTable2, file=snakemake@output[[1]])
#############PrepMeta###############

samplemeta <- read.table(snakemake@input[["premeta"]], header=TRUE, sep="\t")
#samplemeta <- read.table("/storage/mathelierarea/processed/petear/SnakemakeInputFiles/Meta/sampleinfo_TCGA_RNA_seq_cluster.txt", header=TRUE, sep="\t")

###OTHERS THAN BRCA::: _Subtype_mRNA', 'Subtype_Selected' & 'Subtype_Immune_Model_Based
#metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_Selected, Subtype_Immune_Model_Based)


###BRCA_FULL
metaframe <- samplemeta %>% select(donor_id, PAM50, ER.Status, PAM50_genefu, PR.Status, HER2.Final.Status)



#replace numbers as such: 1=Negative, 2=Positive, and LuminalA and LuminalB to "Luminal A" and "Luminal B"
metaframe$ER.Status <- gsub("1","Negative", metaframe$ER.Status)
metaframe$ER.Status <- gsub("2","Positive", metaframe$ER.Status)
metaframe$PAM50 <- gsub("LuminalA","Luminal A", metaframe$PAM50)
metaframe$PAM50 <- gsub("LuminalB","Luminal B", metaframe$PAM50)

#If redundancy in data (E.g same donors multiple times)
#metaSD <- metaframe %>% distinct(icgc_donor_id, .keep_all = TRUE)
metaSD <- metaframe %>% distinct(donor_id, .keep_all = TRUE)

#set rownames
meta2 <- metaSD
rownames(meta2) <- metaSD[,1]


meta <- meta2 %>% select(-1)


###########RemoveNonMetaCpGRows############
#Cange to add metadata rows with NA for CpGs without metadata
#Find CpGs without row in meta
#
#som e i col men ikke i row.

MTList <- colnames(methTable2)
metaList <- rownames(meta)

for (item in MTList)
{
    if (item %in% metaList==TRUE)
    {
        MTList <- MTList[! MTList %in% c(item)]
    }
}

#MTList = list of rows needing to be added to metadata
rn <- row.names(meta)
meta[nrow(meta) + seq_along(MTList), ] <- NA 
row.names(meta) <- c(rn, MTList)


#replace NA with string
#meta <- meta %>% mutate_all(as.character)
meta <- meta %>% replace(is.na(.), "NA")

#Sort meta
meta <- meta[order(row.names(meta)),]

MTList <- colnames(methTable2)

#remove rows not in methTable2 (to be used in topicBina.py)
meta <- meta[row.names(meta) %in% MTList, ]


#free memory:
rm(methTable)
rm(Coords)
rm(samplemeta)
rm(metaframe)
rm(metaSD)


############CISTOPIC########
suppressWarnings(library("cisTopic"))

#Name the PDF after what cancer type and TF:
PDF.name <- snakemake@output[[2]]

#Removing NAs by deletion: (CHECK IF SHOULD USE AMPUTATION!!!
#cpg <- CpG.sample.tab[complete.cases(CpG.sample.tab), ]

cisTopicObject <- createcisTopicObject(is.acc=0.5, min.cells=0, min.regions=0, count.matrix=data.frame(z)) #Set is.acc=0 or is.acc=0.01 | min.cells=0, min.regions=0

cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(4:18), seed=123, nCores=4, addModels=FALSE)
#cisTopicObject <- runModels(cisTopicObject, topic=c(8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60), seed=123, nCores=20, addModels=FALSE)
#cisTopicObject <- runModels(cisTopicObject, topic=c(2,4,6,8), seed=123, nCores=20, addModels=FALSE)

pdf(PDF.name, height=10, width=15)
#cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject)

########## fICA ########
library(fastICA)

t_cisTo <- t(cisTopicObject@selected.model$document_expects)
donorMat  <- matrix(rownames(t_cisTo))
fica <- fastICA(t_cisTo, n.comp= 2, row.norm=FALSE, method = "C")


ficaS <- fica$S

colnames(ficaS) <- c("V1", "V2")
rownames(ficaS) <- c(donorMat)

metaficaS <- meta[rownames(meta) %in% rownames(ficaS), ]

metaficaS <- as.data.frame(metaficaS)
ficaS <- as.data.frame(ficaS)
ficaSMerge <- merge(metaficaS, ficaS, by = 0)

library(RColorBrewer)
################# GGPLOT2 for PAM50 ################
ficaSMergePAMSansNA <- ficaSMerge[!(ficaSMerge$PAM50=="NA"),]

figficaSPAM <- ggplot(ficaSMerge, aes(V1, V2))
figficaSPAM + geom_point(aes(colour = factor(PAM50)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSPAMSansNA <- ggplot(ficaSMergePAMSansNA, aes(V1, V2))
figficaSPAMSansNA + geom_point(aes(colour = factor(PAM50)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

################# GGPLOT2 for ER.Status ################

ficaSMergeERSansNA <- ficaSMerge[!(ficaSMerge$ER.Status=="NA"),]

figficaSER <- ggplot(ficaSMerge, aes(V1, V2))
figficaSER + geom_point(aes(colour = factor(ER.Status)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSERSansNA <- ggplot(ficaSMergeERSansNA, aes(V1, V2))
figficaSERSansNA + geom_point(aes(colour = factor(ER.Status)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################# GGPLOT2 for PAM50_genefu ################
ficaSMergePGSansNA <- ficaSMerge[!(ficaSMerge$PAM50_genefu=="NA"),]

figficaSPG <- ggplot(ficaSMerge, aes(V1, V2))
figficaSPG + geom_point(aes(colour = factor(PAM50_genefu)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSPGSansNA <- ggplot(ficaSMergePGSansNA, aes(V1, V2))
figficaSPGSansNA + geom_point(aes(colour = factor(PAM50_genefu)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################# GGPLOT2 for PR.Status ################
ficaSMergePRSansNA <- ficaSMerge[!(ficaSMerge$PR.Status=="NA"),]

figficaSPR <- ggplot(ficaSMerge, aes(V1, V2))
figficaSPR + geom_point(aes(colour = factor(PR.Status)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSPRSansNA <- ggplot(ficaSMergePRSansNA, aes(V1, V2))
figficaSPRSansNA + geom_point(aes(colour = factor(PR.Status)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


################# GGPLOT2 for HER2.Final.Status ################
ficaSMergeHER2SansNA <- ficaSMerge[!(ficaSMerge$HER2.Final.Status=="NA"),]

figficaSHER2 <- ggplot(ficaSMerge, aes(V1, V2))
figficaSHER2 + geom_point(aes(colour = factor(HER2.Final.Status)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

figficaSHER2SansNA <- ggplot(ficaSMergeHER2SansNA, aes(V1, V2))
figficaSHER2SansNA + geom_point(aes(colour = factor(HER2.Final.Status)), size = 2) + scale_color_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())


######################################################
#####Umap and further cisTopic
######################################################
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"regionScore by NormTop", pos=4)


cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
#par(mfrow=c(3,5))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status","PAM50_genefu","PR.Status","HER2.Final.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by probability (by NormTop)", pos=4)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by Z-score (by NormTop)", pos=4)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by NormTop (by NormTop)", pos=4)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


######################################################
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"regionScore by probability", pos=4)

cisTopicObject <- getRegionsScores(cisTopicObject, method='Probability', scale=TRUE)
#par(mfrow=c(3,5))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
cisTopicObject <- runUmap(cisTopicObject, target="cell", n_components=3)
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c("PAM50","ER.Status","PAM50_genefu","PR.Status","HER2.Final.Status"), cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by probability (by probability)", pos=4)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by Z-score (by probability)", pos=4)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Region by NormTop (by probability)", pos=4)

cisTopicObject <- runUmap(cisTopicObject, target ='region')

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)

cisTopicObject <- runUmap(cisTopicObject, target ='region', n_components=3)

plotFeatures(cisTopicObject, method='Umap', target='region', topic_contr='NormTop', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=3, legend=FALSE, col.low='blue', col.mid='white', col.high='red', cex.dot = 1.5)


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

#Not needed
#write.csv(fica, file=snakemake@output[[7]])

#getBedFiles(cisTopicObject, path='CisTopicBed_BRCA_SNA')
 
write.csv(meta, file=snakemake@output[[7]])

#Create empty dataframe:
df <- data.frame(matrix(ncol=2, nrow=1))
colnames(df) <- c("chrPos","TopicX")

#merge dataframes
for (attr in attributes(cisTopicObject@binarized.cisTopics)$names)
{
    makeDF <- data.frame(cisTopicObject@binarized.cisTopics[attr])
    rowColDF <- tibble::rownames_to_column(makeDF, "chrPos")
    df <- merge(df, rowColDF, by="chrPos", all=TRUE)
}
df$TopicX <- NULL

write.csv(df@binarized.cisTopics, file=snakemake@output[[7]])

write.csv(meta, file=snakemake@output[[8]])