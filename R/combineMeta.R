combineMeta <- function(df, meta, 
                    cloneTypes=c(None = 0, Single = 1, Small = 5, Medium = 10, Large = 20, Hyperexpanded = 150)) {
    
        data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", "CTaa", "ID")]))
        data2 <- data2[data2[,"barcode"] %in% rownames(meta), ]
        data2 <- as.data.frame(data2 %>% group_by(data2[,"CTaa"], 
                                                  data2[,"ID"]) %>% summarise(Frequency = n()))
        colnames(data2)[c(1,2)] <- c("CTaa", "condition")
        x <- unique(data[,"ID"])
        Con.df  <- NULL
        for (i in seq_along(x)) {
            sub1 <- subset(data, data[,"ID"] == x[i])
            sub2 <- subset(data2, data2[,"condition"] == x[i])
            merge <- merge(sub1, sub2, by="CTaa")
            Con.df <- rbind.data.frame(Con.df, merge) } 
        
        
    Con.df$cloneType <- NA
    for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <- 
        paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], 
               ' < X <= ', cloneTypes[x], ')') }
    for (i in 2:length(cloneTypes)) { Con.df$cloneType <- 
        ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency 
               <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
    PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt", 
                                "CTaa", "CTstrict", "Frequency", "cloneType")])
    rownames(PreMeta) <- PreMeta$barcode
    
    merged <- merge(PreMeta, meta, by = "row.names")
    rownames(merged) <- merged$barcode
    merged <- merged[,-c(1:2)]
    return(merged) 
}
