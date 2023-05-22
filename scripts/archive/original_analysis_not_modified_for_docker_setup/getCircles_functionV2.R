

getCircles <- function(data, proportion = FALSE, clonotypesOnly = FALSE) {

  # Reformat data to generate a clonotype by cluster matrix
  temp <- data[, c("CTaa", "Clusters_l1")]
  dTest <- reshape2::dcast(temp, CTaa ~ Clusters_l1)
  dTest <- dTest[apply(dTest[, -1], 1, function(x) !all(x == 0)), ]
  dTest <- dTest[-1]
  total <- nrow(dTest)

  # Initialize output matrix
  matrix_out <- matrix(0, ncol = ncol(dTest), nrow = ncol(dTest))

  # Calculate overlap of clonotypes or # of cells with overlapping clonotypes for each clusterX_clusterY pair
  if (clonotypesOnly == TRUE) {

    # Reduce data to a binary matrix representing presence or absence of each clonotype in a cluster
    dTest[dTest > 1] <- 1

    for (x in seq_len(ncol(dTest))) {
      for (y in seq_len(ncol(dTest))) {
        if (x == y) {
          # when x == y calculate the number of clonotypes which are unique to the cluster and have no overlap
          matrix_out[x, x] <- sum(dTest[, x] >= 1 & rowSums(dTest[, -x]) == 0)
        } else {
          # in all other instances, count the number of clonotypes which are present in both clusters under consideration
          matrix_out[y, x] <- sum(dTest[, x] >= 1 & dTest[, y] >= 1)
        }
      }
    }
  } else {
    for (x in seq_len(ncol(dTest))) {
      for (y in seq_len(ncol(dTest))) {
        temp.vec <- NULL # ensure no carry over from previous loop // unsure if needed
        if (x == y) {
          # when x == y calculate the number of clonotypes which are unique to the cluster and have no overlap
          temp.vec <- which(dTest[, x] >= 1 & rowSums(dTest[, -x]) == 0)
          if (length(temp.vec) > 0) { # require a catch if() to detect if temp.vec is empty
            matrix_out[x, x] <- sum(dTest[temp.vec, x])
          } else {
            matrix_out[y, x] <- 0
          }
        } else {
          # get the every row where clonotypes overlap between x & y
          temp.vec <- which(dTest[, x] >= 1 & dTest[, y] >= 1)
          if (length(temp.vec) > 0) { # require a catch if() to detect if temp.vec is empty
            matrix_out[y, x] <- sum(dTest[temp.vec, c(x, y)]) # Sum the number of cells across both X and Y columns
          } else {
            matrix_out[y, x] <- 0 # there are no overlaps between the clusters
          }
        }
      }
    }
  }
  colnames(matrix_out) <- colnames(dTest)
  rownames(matrix_out) <- colnames(dTest)

  output <- data.frame(
    from = rep(rownames(matrix_out), times = ncol(matrix_out)),
    to = rep(colnames(matrix_out), each = nrow(matrix_out)),
    value = as.vector(matrix_out),
    stringsAsFactors = FALSE
  )

  # Reorder columns to eliminate redundant comparisons
  for (k in 1:nrow(output)) {
    max <- order(output[k, 1:2])[1] # which is first alphabetically
    max <- output[k, max]
    min <- order(output[k, 1:2])[2] # which is second alphabetically
    min <- output[k, min]
    output[k, 1] <- max
    output[k, 2] <- min
  }
  unique.rows <- rownames(unique(output[, 1:2])) # removing redundant comparisons
  output <- output[rownames(output) %in% unique.rows, ]
  
  if (proportion == TRUE) {
    output$value <- output$value / total
  }

  return(output)
}


getIntegratedCircle <- function(data, proportion = FALSE, clonotypesOnly = FALSE) {
  #browser()
  output <- NULL
  test <- data[, c("CTaa", "Clusters_l1", "condition")]
  totalUS <- table(subset(test, !is.na(CTaa))$condition)[1]
  totalStim <- table(subset(test, !is.na(CTaa))$condition)[2]

  dTest <- reshape2::dcast(test, CTaa ~ condition + Clusters_l1)
  dTest <- dTest[apply(dTest[, -1], 1, function(x) !all(x == 0)), ]
  dTest <- dTest[, -1]



  # Initialize output matrix
  matrix_out <- matrix(0, ncol = ncol(dTest), nrow = ncol(dTest))

  # Calculate overlap of clonotypes or # of cells with overlapping clonotypes for each clusterX_clusterY pair
  if (clonotypesOnly == TRUE) {

    # Reduce data to a binary matrix representing presence or absence of each clonotype in a cluster
    dTest[dTest > 1] <- 1

    for (x in seq_len(ncol(dTest))) {
      for (y in seq_len(ncol(dTest))) {
        if (x == y) {
          # when x == y calculate the number of clonotypes which are unique to the cluster and have no overlap
          matrix_out[x, x] <- sum(dTest[, x] >= 1 & rowSums(dTest[, -x]) == 0)
        } else {
          # in all other instances, count the number of clonotypes which are present in both clusters under consideration
          matrix_out[y, x] <- sum(dTest[, x] >= 1 & dTest[, y] >= 1)
        }
      }
    }
  } else {
    for (x in seq_len(ncol(dTest))) {
      for (y in seq_len(ncol(dTest))) {
        temp.vec <- NULL # ensure no carry over from previous loop // unsure if needed
        if (x == y) {
          # when x == y calculate the number of clonotypes which are unique to the cluster and have no overlap
          temp.vec <- which(dTest[, x] >= 1 & rowSums(dTest[, -x]) == 0)
          if (length(temp.vec) > 0) { # require a catch if() to detect if temp.vec is empty
            matrix_out[x, x] <- sum(dTest[temp.vec, x])
          } else {
            matrix_out[y, x] <- 0
          }
        } else {
          # get the every row where clonotypes overlap between x & y
          temp.vec <- which(dTest[, x] >= 1 & dTest[, y] >= 1)
          if (length(temp.vec) > 0) { # require a catch if() to detect if temp.vec is empty
            matrix_out[y, x] <- sum(dTest[temp.vec, c(x, y)]) # Sum the number of cells across both X and Y columns
          } else {
            matrix_out[y, x] <- 0 # there are no overlaps between the clusters
          }
        }
      }
    }
  }


  colnames(matrix_out) <- colnames(dTest)
  rownames(matrix_out) <- colnames(dTest)



  output <- data.frame(
    from = rep(rownames(matrix_out), times = ncol(matrix_out)),
    to = rep(colnames(matrix_out), each = nrow(matrix_out)),
    value = as.vector(matrix_out),
    stringsAsFactors = FALSE
  )


  # Reorder columns to eliminate redundant comparisons
  for (k in 1:nrow(output)) {
    max <- order(output[k, 1:2])[1] # which is first alphabetically
    max <- output[k, max]
  min <- order(output[k, 1:2])[2] # which is second alphabetically
    min <- output[k, min]
    output[k, 1] <- max
    output[k, 2] <- min
  }
  unique.rows <- rownames(unique(output[, 1:2])) # removing redundant comparisons
  output <- output[rownames(output) %in% unique.rows, ]
  

  output$from_group <- stringr::str_split(output$from, "_", simplify = T, n = 2)[, 1]
  output$to_group <- stringr::str_split(output$to, "_", simplify = T, n = 2)[, 1]
  output$from <- stringr::str_split(output$from, "_", simplify = T, n = 2)[, 2]
  output$to <- stringr::str_split(output$to, "_", simplify = T, n = 2)[, 2]
  
  if (proportion == TRUE) {
    output$value <- ifelse(output$from_group == "US", output$value / totalUS, output$value / totalStim)
  }

  return(output)
}
