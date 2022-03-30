
#Load libraries
library(ComplexHeatmap)


#Load data
df <- read.table(snakemake@input[["py_DF"]], header=TRUE, sep="\t")
meta <- read.table(snakemake@input[["py_meta"]], header=TRUE, sep="\t")

# df <- read.table("py_df.csv", header=TRUE, sep=",")
# meta <- read.table("py_df.csv", header=TRUE, sep="\t")