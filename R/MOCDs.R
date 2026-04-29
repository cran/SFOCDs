#' M Array Method for Optimal Covariate Designs (OCDs)
#'
#' @param design Input a design in matrix format and block size k such that (k+1) is a prime number.
#'@description Consider Mutually Orthogonal Latin squares of order s x s, extract first column from it and make a new matrix called Initial block sequence matrix L of order s x s-1. Remove the last row from the L matrix and obtain the incidence matrix of it keeping zeros to the positions corresponding to the elements that were present in the deleted row of L, and then remove the row that contains only non-zero elements, the square matrix thus formed is the M matrix. From the columns of M matrix, choose((s-1),2) pairs are possible. Each of these column pairs is then superimposed to N. Through this method choose((s-1),2) W matrices can be developed.
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
#' mat<-matrix(c(1,2,3,4,1,2,4,5,1,3,4,5,2,3,4,5),nrow=4,byrow=TRUE)
#' MOCDs(mat)
#' @export
MOCDs<-function(design){

  ####incidence matrix
  design<-as.matrix(design)
  #design<-rearrange_backtrack(design)
  # the input complentary design
  #####incidence matrix
  k=ncol(design)
  block <- nrow(design)
  trt <- max(design)
  #######incidence matrix
  incidence_matrix <- matrix(0, nrow = block, ncol = trt)
  for (i in 1:block) {
    for (j in 1:ncol(design)) {
      incidence_matrix[i, design[i, j]] <- 1
    }
  }
  incidence_matrix_transpose <- t(incidence_matrix)
  inc<-incidence_matrix_transpose

  ####initial block
  MOLS<-function(v){
    MOLS<-list()
    for(i in 1:(v-1)){
      mols<-NULL
      seq<-c(1:v)
      j=0
      repeat{
        mols<-rbind(mols,c(seq+j))
        j=j+i
        if(nrow(mols)==v){
          mols=mols%%v
          mols[mols==0]<-v
          MOLS<-append(MOLS,list(mols))
          break
        }
      }
    }
    return(MOLS)
  }
  mols<-MOLS(k+1)
  del_col_mols<-lapply(mols,function(mat)mat[,1])
  L<-do.call(cbind,del_col_mols)

  L_del <- L[-(nrow(L)-1), ]


  # Number of treatments (v) and blocks (b)

  b <- ncol(L)
  v <- max(L)
  # Initialize  matrix (v x b)####### m matrix
  m1 <- matrix(0, nrow = v, ncol = b)

  # Fill incidence matrix
  for (j in 1:b) {
    for (i in 1:nrow(L_del)) {
      treatment <- L_del[i, j]
      m1[treatment, j] <- 1
    }
  }
  m2 <- m1[apply(m1, 1, function(x) any(x == 0)), ]

  #print("Matrix after deleting rows with no zeros:")
  #print(m2)
  m <- ifelse(m2 == 0, -1, m2)
  col_pairs <- combn(ncol(m), b/2, simplify = FALSE)
  mpair <- do.call(cbind, lapply(col_pairs, function(cols) {
    apply(m[, cols], 1, prod)
  }))


  W_mats<-list()
  pos_ist<-list()
  for(i in 1:ncol(mpair)){
    minus_pos<-which(mpair[,i]==-1)
    pos_ist<-append(pos_ist,list(design[,c(minus_pos)]))
  }
  incid<-incidence_matrix_transpose

  for(i in 1:length(pos_ist)){
    for(j in 1:nrow(pos_ist[[i]])){
      incid[c(pos_ist[[i]][j,]),j]<--1
    }
    W_mats<-append(W_mats,list(incid))
    incid<-incidence_matrix_transpose
  }

  for (i in seq_along(m)) {

    # skip if m[[i]] is NULL, has no rows, or not a data.frame/matrix
    if (is.null(m[[i]]) || !is.matrix(m[[i]]) && !is.data.frame(m[[i]]) || nrow(as.data.frame(m[[i]])) == 0) {
      next
    }

    for (j in 1:nrow(as.data.frame(m[[i]]))) {
      all_pos_incidence <- abs(m[[i]][j, ])
      incd[c(all_pos_incidence), j] <- incd[c(all_pos_incidence), j] * -1
    }

    W_mats <- append(W_mats, list(incd))
  }

  interproductsum <- list()

  k1 <- ncol(m)              # order of m matrix
  r1 <- k1 / 2                # number of W matrices to combine

  # all combinations of W matrices of size r
  comb_idx <- combn(length(W_mats), r1)

  for(t in 1:ncol(comb_idx)){

    idx <- comb_idx[, t]

    # multiply all selected W matrices
    AxB <- Reduce("*", W_mats[idx])

    # row-wise sum
    intprdsum <- AxB %*% rep(1, ncol(AxB))

    interproductsum <- append(interproductsum,
                              list(intprdsum, idx))
  }

  interproductsum <- list()

  for(t in 1:ncol(comb_idx)){

    idx <- comb_idx[, t]

    AxB <- Reduce("*", W_mats[idx])

    intprdsum <- AxB %*% rep(1, ncol(AxB))

    # store as ONE object (not two separate entries)
    interproductsum[[t]] <- list(
      sum = intprdsum,
      idx = idx
    )
  }

  total_count <- length(interproductsum)
  # Each column of mpair corresponds to one W
  W_sets <- lapply(1:ncol(mpair), function(i) which(mpair[, i] == -1))
  fixed_index <- 1

  valid_indices <- which(sapply(W_sets, function(x) fixed_index %in% x))

  valid_W <- W_mats[valid_indices]

  length(valid_W)
  length(W_mats)
  lm<-list("W matrices",valid_W,"Interproduct Sums",interproductsum,length(valid_W),length(W_mats))
  return(lm)
}

