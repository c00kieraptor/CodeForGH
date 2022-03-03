
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


    head(df2)
    
    df <- merge(df, cisTopicObject@binarized.cisTopics[attr], by=rownames, all=TRUE)
     
    df <- df[, !duplicated(colnames(df), fromLast = TRUE)]
    
    filename = sprintf("This is where a %s goes.", a)
    write.table(cisTopicObject@binarized.cisTopics[attr], file=attr)    

print(attr)
print(cisTopicObject@binarized.cisTopics[attr]) 

df <- merge(df, cisTopicObject@binarized.cisTopics$Topic1, by=0, all=TRUE)
