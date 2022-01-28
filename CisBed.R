suppressWarnings(library("cisTopic"))
args = commandArgs(trailingOnly=TRUE)

argument1 <- args[1]
cisTopicObject <- readRDS(argument1)

pathcmd <- args[2]
getBedFiles(cisTopicObject, path=pathcmd)
