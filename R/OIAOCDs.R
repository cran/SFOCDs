#' Orthogonal Integer Array Method for Optimal Covariate Designs (OCDs)
#'
#' @param design Input a design in matrix format and block size k such that it is an odd number.
#'@description Consider OIA of order same as the block size of the required design. Superimpose each array separately into the incidence matrix (N) to get the W matrices. The maximum number of W matrices depends on the column order of OIA.
#'
#' @returns Generates W matrices and Inter product sums of W matrices.
#'@author
#'Neethu RS
#'
#'Cini Varghese
#'
#'Mohd Harun
#'
#'Anindita Datta
#'
#'Ashutosh Dalal
#' @references
#' Das, K., N. K. Mandal, and B. K. Sinha. (2003) <https://doi.org/10.1016/S0378-3758(02)00151-9>. Optimal experimental designs for models with covariates. Journal of Statistical Planning and Inference 115(1): 273-285.
#'
#' Bansal, N., and D. K. Garg. (2022)<https://doi.org/10.1007/s42519-022-00244-0>. Optimum covariate designs for three associate PBIB designs. Journal of Statistical Theory and Practice 16(3): 1-15.
#'
#' @examples
#' library(SFOCDs)
#' mat<-matrix(c(1,2,3,1,2,4,1,2,5,1,3,4,1,3,5,1,4,5,2,3,4,2,3,5,2,4,5,3,4,5),nrow=10,byrow=TRUE)
#' OIAOCDs(mat)
#' @export
OIAOCDs<-function(design){
design<-as.matrix(design)
block <- nrow(design)
trt <- max(design)
incidence_matrix <- matrix(0, nrow = block, ncol = trt)
for (i in 1:block) {
  for (j in 1:ncol(design)) {
    incidence_matrix[i, design[i, j]] <- 1
  }
}
incidence_matrix_transpose <- t(incidence_matrix)
inc<-incidence_matrix_transpose

#####oia
orthogonal_integer_array <- function(n) {
  # Construct orthogonal basis with row sum zero
  # Use Helmert matrix
  H <- matrix(0, n, n-1)
  for (k in 1:(n-1)) {
    H[1:k, k]   <- 1
    H[k+1, k]   <- -k
  }
  return(H)
}
oia<-orthogonal_integer_array(ncol(design))

####superimpose
num_w <- ncol(oia)

# List to store W matrices
W_mats <- vector("list", ncol(oia))

for (w_idx in 1:ncol(oia)) {
  W <- inc  # start with incidence matrix

  # Loop over blocks and positions
  for (block in 1:nrow(design)) {
    treatments <- design[block, ]
    for (pos in seq_along(treatments)) {
      treatment <- treatments[pos]
      W[treatment, block] <- oia[pos, w_idx]
    }
  }

  W_mats[[w_idx]] <- W
}
interproductsum<-list()
#indices<-list()
for(i in 1:(length(W_mats)-1)){
  for(j in (i+1):length(W_mats)){
    AxB<-W_mats[[i]]*W_mats[[j]]
    intprdsum<-AxB%*%rep(1,ncol(AxB))
    interproductsum<-append(interproductsum,list(intprdsum,c(i,j)))
    #indices<-append(indices,)
  }
}
new_list1 <- lapply(W_mats, function(W) {
  list(
    w = W,
    row_sum = rowSums(W),
    col_sum = colSums(W)
  )
})

new_list<-split(interproductsum,rep(1:(length(interproductsum)/2),each=2))
lm<-list("W matrices",new_list1,"Interproduct Sums",new_list)
return(lm)
}
