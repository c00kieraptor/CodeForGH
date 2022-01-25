DfForHM <- read.csv("ProbeTopicScore.csv" )

library(ggplot2)
library(reshape)

DfForHM <- subset(DfForHM, select = -X )
pdf("testsfgege.pdf", height=10, width=15)

dev.off()