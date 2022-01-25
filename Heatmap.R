DfForHM <- read.csv("ProbeTopicScore.csv" )

library(ggplot2)
library(reshape)

DfForHM <- subset(DfForHM, select = -X )

DfMelted <- melt(DfForHM)

names(DfMelted)[names(DfMelted) == "Unnamed..0"] <- "Coordinate"
names(DfMelted)[names(DfMelted) == "variable"] <- "Topics"

pdf("HeatmapTop20percentProbes.pdf", height=10, width=15)

heatmap <- ggplot(DfMelted, aes(Topics, Coordinate, fill= value)) + 
 geom_tile()

heatmap

dev.off()