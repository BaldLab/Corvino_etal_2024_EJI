getCircles <- function(data, proportion = FALSE, clonotypesOnly = FALSE) {

    test <- data[, c("CTaa", "seurat_clusters")]
    dTest <- reshape2::dcast(test, CTaa ~ seurat_clusters)
    dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
    dTest <- dTest[-1]
    total <- nrow(dTest)
    ##This will prevent counting clonotypes by cell number
    if (clonotypesOnly == TRUE) {
        dTest[dTest > 1] <- 1
    }
    matrix_out <- matrix(ncol = ncol(dTest), nrow = ncol(dTest))
    for (x in seq_len(ncol(dTest))) {
        for (y in seq_len(ncol(dTest)) ){
            matrix_out[y,x] <- length(which(dTest[,x] >= 1 & dTest[,y] >= 1))
        }
    }
    colnames(matrix_out) <- colnames(dTest)
    rownames(matrix_out) <- colnames(dTest)
    #Need to subtract extra cells - will take the difference of the sum of the column minus and the respective cell and subtract that from the respective cell
    for (y in seq_len(ncol(matrix_out))) {
        matrix_out[y,y] <- matrix_out[y,y] - (sum(matrix_out[,y])-matrix_out[y,y])
        if (matrix_out[y,y] < 0) {
            matrix_out[y,y] <- 0
        }
    }
    # Reduces the clonotypes in half - this will allow for accurate depiction by total number of cells in cluster
    for (i in seq_len(ncol(matrix_out))) {
        for (j in seq_len(ncol(matrix_out))) {
            matrix_out[i,j] <- as.integer(matrix_out[i,j]/2)
        }
    }
    output <- data.frame(from = rep(rownames(matrix_out), times = ncol(matrix_out)),
                         to = rep(colnames(matrix_out), each = nrow(matrix_out)),
                         value = as.vector(matrix_out),
                         stringsAsFactors = FALSE)
    # Reorder columns to eliminate redundant comparisons
    for (k in 1:nrow(output)) {
        max <- order(output[k,1:2])[1] #which is first alphabetically
        max <- output[k,max]
        min <- order(output[k,1:2])[2] #which is second alphabetically
        min <- output[k,min]
        output[k,1] <- max
        output[k,2] <- min
    }
    unique <- rownames(unique(output[,1:2])) #removing redundant comparisons
    output <- output[rownames(output) %in% unique, ]
    if (proportion == TRUE) {
        output$value <- output$value/total
    } 
    
    return(output)
}

getIntegratedCircle <- function(data, proportion = FALSE, clonotypesOnly = FALSE) {
    output <- NULL
    test <- data[, c("CTaa", "seurat_clusters", "condition")]
    totalUS <- table(subset(test, !is.na(CTaa))$condition)[1]
    totalStim <- table(subset(test, !is.na(CTaa))$condition)[2]
    
    dTest <- reshape2::dcast(test, CTaa ~ condition + seurat_clusters)
    dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
    dTest <- dTest[,-1]
    ##This will prevent counting clonotypes by cell number
    if (clonotypesOnly == TRUE) {
        dTest[dTest > 1] <- 1
    }
    matrix_out <- matrix(ncol = ncol(dTest), nrow = ncol(dTest))
    for (x in seq_len(ncol(dTest))) {
        for (y in seq_len(ncol(dTest)) ){
            matrix_out[y,x] <- length(which(dTest[,x] >= 1 & dTest[,y] >= 1))
        }
    }
    colnames(matrix_out) <- colnames(dTest)
    rownames(matrix_out) <- colnames(dTest)
    #Need to subtract extra cells - will take the difference of the sum of the column minus and the respective cell and subtract that from the respective cell
    for (y in seq_len(ncol(matrix_out))) {
        matrix_out[y,y] <- matrix_out[y,y] - (sum(matrix_out[,y])-matrix_out[y,y])
        if (matrix_out[y,y] < 0) {
            matrix_out[y,y] <- 0
        }
    }
    # Reduces the clonotypes in half - this will allow for accurate depiction by total number of cells in cluster
    for (i in seq_len(ncol(matrix_out))) {
        for (j in seq_len(ncol(matrix_out))) {
            matrix_out[i,j] <- as.integer(matrix_out[i,j]/2)
        }
    }
    output <- data.frame(from = rep(rownames(matrix_out), times = ncol(matrix_out)),
                         to = rep(colnames(matrix_out), each = nrow(matrix_out)),
                         value = as.vector(matrix_out),
                         stringsAsFactors = FALSE)
    # Reorder columns to eliminate redundant comparisons
    for (k in 1:nrow(output)) {
        max <- order(output[k,1:2])[1] #which is first alphabetically
        max <- output[k,max]
        min <- order(output[k,1:2])[2] #which is second alphabetically
        min <- output[k,min]
        output[k,1] <- max
        output[k,2] <- min
    }
    unique <- rownames(unique(output[,1:2])) #removing redundant comparisons
    output <- output[rownames(output) %in% unique, ]
    #output <- output[which(output$from != output$to),]
    
    
    output$from_group <- stringr::str_split(output$from, "_", simplify = T, n=2)[,1]
    output$to_group <- stringr::str_split(output$to, "_", simplify = T, n=2)[,1]
    output$from<- stringr::str_split(output$from, "_", simplify = T, n=2)[,2]
    output$to <- stringr::str_split(output$to, "_", simplify = T, n=2)[,2]
    if (proportion == TRUE) {
        output$value <- ifelse(output$from == "US", output$value/totalUS, output$value/totalStim)
    }   
        
    return(output)
}
