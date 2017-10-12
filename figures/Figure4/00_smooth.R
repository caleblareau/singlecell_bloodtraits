library(RANN)
library(irlba)

# Quick call to irlba for SVD/PCA
run_pca <- function(data, n_components){
 pr <- irlba::prcomp_irlba(data, n = n_components) 
 return(pr$rotation)
}

compute_markov <- function(data, knn=10, epsilon=1, distance_metric='euclidean', knn_autotune = 10){
  
  # KNN Graph; Should add more distance metrics here?
  if(distance_metric == "euclidean"){
    di <- RANN::nn2(data, k = knn+1)
  }
  
  # Auto tune
  if(knn_autotune > 0){
    col <- min(knn, knn_autotune) + 1
    di$nn.dists <- sweep(di$nn.dists, 1, di$nn.dists[,col], FUN="/")
  }
  
  # Set up matrix; include zeros added later in python version
  if(epsilon > 0){
    W <- Matrix::sparseMatrix(i = rep(di$nn.idx[,1],knn+1), j = as.vector(di$nn.idx), x = as.vector(di$nn.dists))
  } else {
    W <- Matrix::sparseMatrix(i = rep(di$nn.idx[,1],knn+1), j = as.vector(di$nn.idx), x = rep(1,(knn + 1)*dim(di$nn.idx)[1]))
  }
  
  # Symmetry
  W <- W + t(W)
  
  # Apply Gaussian Kernel
  if(epsilon > 0){
    W <-  Matrix::sparseMatrix(i = summary(W)$i, j = summary(W)$j, x = exp(-1*(summary(W)$x/(epsilon^2))))
  }
  
  # Return normalized transition matrix
  return(sweep(W, 1, Matrix::rowSums(W), FUN="/"))
  
}

"%^%" <- function(base,power){
    out <- diag(nrow(base))
    while(power > 1){
        if(power %% 2 == 1){
            out <- out %*% base
        }
        base <- base %*% base
        power <- power %/% 2
    }
    out %*% base
}

impute_fast <- function(data, L, t, rescale_percent=99, L_t=NULL, tprev=NULL){
  cat('Starting Imputation\n')
  if(is.null(L_t)){
    L_t <- L%^%t
  } else {
    L_t <- L_t %*% L%^%(t-tprev)
  }

  
  newdata <- data %*% L_t
  return(newdata)
  if(rescale_percent != 0){
    # Per-cell rescale_%tile; new maxs
    rp_cells <- apply(data, 2, quantile, probs = rescale_percent, na.rm = TRUE)
    nm_cells <- apply(newdata, 2, max)
    newdata <- sweep(newdata, 2, (rp_cells/nm_cells), FUN="/")
  }
  
  
  
}

# Main Method
smoothPCAmat <- function(pcmat, kernel='gaussian', n_pca_components=20, 
                  t=6, knn=30, knn_autotune=10, epsilon=1, rescale=0.995, k_knn=50, perplexity=30){
 
  # No current use for: random_pca, perplexity in current implementation
  stopifnot(kernel == 'gaussian')

  # Compute Markov Matrix
  if(kernel == "gaussian"){
    L <- compute_markov(pcmat, knn=knn, epsilon=epsilon, distance_metric='euclidean', knn_autotune=knn_autotune) 
  }
   return(L)
}

