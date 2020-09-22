<<<<<<< HEAD
getCircles <- function(data, proportion = FALSE) {

=======
getCircles <- function(data, proportion = FALSE, clonotypesOnly = FALSE) {
# browser()
>>>>>>> fd1656f521318b87ab3107f261047301404a6970
    test <- data[, c("CTaa", "seurat_clusters")]
    dTest <- reshape2::dcast(test, CTaa ~ seurat_clusters)
    dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
    dTest <- dTest[-1]
    total <- nrow(dTest)
<<<<<<< HEAD

    matrix_out <- matrix(ncol = ncol(dTest), nrow = ncol(dTest))
=======
    
    # This will prevent counting clonotypes by cell number
    if (clonotypesOnly == TRUE) {
        dTest[dTest > 1] <- 1
    }
    
    # Create matrix of clonotype overlap between clusters - this reduces the dataset to a # of clonotypes overlap and ignores cell # info!!!
    matrix_out <- matrix(0, ncol = ncol(dTest), nrow = ncol(dTest))
>>>>>>> fd1656f521318b87ab3107f261047301404a6970
    for (x in seq_len(ncol(dTest))) {
        for (y in seq_len(ncol(dTest)) ){
            matrix_out[y,x] <- length(which(dTest[,x] >= 1 & dTest[,y] >= 1)) # need to consider adding up values obtained from which() to keep cell # info
        }
    }
    
    colnames(matrix_out) <- colnames(dTest)
    rownames(matrix_out) <- colnames(dTest)
<<<<<<< HEAD
    
    #Need to subtract extra cells - will take the difference of the sum of the column minus and the respective cell and subtract that from the respective cell
=======
    # Need to subtract extra cells - will take the difference of the sum of the column minus and the respective cell and subtract that from the respective cell
>>>>>>> fd1656f521318b87ab3107f261047301404a6970
    for (y in seq_len(ncol(matrix_out))) {
        matrix_out[y,y] <- matrix_out[y,y] - (sum(matrix_out[,y])-matrix_out[y,y])
        if (matrix_out[y,y] < 0) {
            matrix_out[y,y] <- 0
        }
    }
<<<<<<< HEAD

=======
    # Reduces the clonotypes in half - this will allow for accurate depiction by total number of cells in cluster
    # Is a loop necessary - can't just dataset <- dataset/2
    for (i in seq_len(ncol(matrix_out))) {
        for (j in seq_len(ncol(matrix_out))) {
            matrix_out[i,j] <- as.integer(matrix_out[i,j]/2)
        }
    }
    
>>>>>>> fd1656f521318b87ab3107f261047301404a6970
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
    unique.rows <- rownames(unique(output[,1:2])) #removing redundant comparisons
    output <- output[rownames(output) %in% unique.rows, ]
    if (proportion == TRUE) {
        output$value <- output$value/total
    } 
    
    return(output)
}

<<<<<<< HEAD
getIntegratedCircle <- function(data, proportion = FALSE) {
=======
getIntegratedCircle <- function(data, proportion = FALSE, clonotypesOnly = FALSE) {
  #browser()
>>>>>>> fd1656f521318b87ab3107f261047301404a6970
    output <- NULL
    test <- data[, c("CTaa", "seurat_clusters", "condition")]
    totalUS <- table(subset(test, !is.na(CTaa))$condition)[1]
    totalStim <- table(subset(test, !is.na(CTaa))$condition)[2]
    
    dTest <- reshape2::dcast(test, CTaa ~ condition + seurat_clusters)
    dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
    dTest <- dTest[,-1]
<<<<<<< HEAD

    matrix_out <- matrix(0,ncol = ncol(dTest), nrow = ncol(dTest))
    for (x in seq_len(ncol(dTest))) {
        for (y in seq_len(ncol(dTest))){
=======
    ##This will prevent counting clonotypes by cell number
    if (clonotypesOnly == TRUE) {
        dTest[dTest > 1] <- 1
    }
    matrix_out <- matrix(0, ncol = ncol(dTest), nrow = ncol(dTest))
    for (x in seq_len(ncol(dTest))) {
        for (y in seq_len(ncol(dTest))){
            matrix_out[y,x] <- matrix_out[y,x] + length(which(dTest[,x] >= 1 & dTest[,y] >= 1 & rowSums(dTest[,-c(x,y)]) == 0))
            if (dTest[,x] >= 1 & dTest[,y] >= 1 & rowSums(dTest[,-c(x,y)]) != 0) {
                others <- which(dTest[,x] >= 1 & dTest[,y] >= 1 & rowSums(dTest[,-c(x,y)]) != 0)
                for (z in seq_along(others)) {
                    float_integer <- which(dTest[z,] > 1)
                    float_integer <- na.omit(float_integer[float_integer != c(x,y)])
                    if (x == y) {
                        for (a in seq_along(float_integer)) {
                        matrix_out[a,x] = matrix_out[a,x] + 1
                        matrix_out[x,a] = matrix_out[x,a] + 1
                        }
                    } else {
                    for (a in seq_along(float_integer)) {
                        matrix_out[a,x] = matrix_out[a,x] + 1
                        matrix_out[a,y] = matrix_out[a,y] + 1
                        matrix_out[x,a] = matrix_out[x,a] + 1
                        matrix_out[y,a] = matrix_out[y,a] + 1
                    }
                    } 
                }
            }
        
        if (x == y) {
>>>>>>> fd1656f521318b87ab3107f261047301404a6970
            matrix_out[y,x] <- length(which(dTest[,x] >= 1 & dTest[,y] >= 1))
        }
    }
    colnames(matrix_out) <- colnames(dTest)
    rownames(matrix_out) <- colnames(dTest)
    
    # Need to subtract extra cells - will take the difference of the sum of the column minus and the respective cell and subtract that from the respective cell
    for (y in seq_len(ncol(matrix_out))) {
        matrix_out[y,y] <- matrix_out[y,y] - (sum(matrix_out[,y])-matrix_out[y,y])
        if (matrix_out[y,y] < 0) {
            matrix_out[y,y] <- 0
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
}}
