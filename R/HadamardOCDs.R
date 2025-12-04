#' Hadamard Method for Optimal Covariate Designs (OCDs)
#'
#' @param design Input a design in matrix format and block size k is multiple of 4.
#'@description Construct Hadamard matrix H_k = (1, h_1 , h_2,..., h_k-1) where k is the block size of the required design. Then superimpose each columns of H_k leaving the first column which is in natural order separately into the N matrix to get the W matrices. The maximum number of W matrices will be k-1.
#'
#' @returns Generates W matrices and Inter product sums of W matrices.
#' @export
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
#' HadamardOCDs(mat)
HadamardOCDs<-function(design){
  ####incidence matrix
  design<-as.matrix(design)
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
  ######

  ##Hadamard##
  # Function to generate Hadamard matrix of order k
  hadamard_matrix <- function(k) {
    if (k == 1) {
      return(matrix(1, nrow = 1, ncol = 1))  # Base case, H_1 = [1]
    } else {
      # Recursively build the Hadamard matrix for order 2^n
      H <- hadamard_matrix(k / 2)

      # Combine matrices to form the Hadamard matrix of order k
      top_left <- H
      top_right <- H
      bottom_left <- H
      bottom_right <- -H

      # Combine the blocks into a larger matrix
      return(rbind(
        cbind(top_left, top_right),
        cbind(bottom_left, bottom_right)
      ))
    }
  }
  k <- ncol(design)
  H1 <- hadamard_matrix(k)
  H<-H1[,-1]


  #####super imposing part
  W_mats<-list()
  pos_ist<-list()
  for(i in 1:ncol(H)){
    minus_pos<-which(H[,i]==-1)
    pos_ist<-append(pos_ist,list(design[,c(minus_pos)]))
  }
  #######
  incid<-incidence_matrix_transpose

  for(i in 1:length(pos_ist)){
    for(j in 1:nrow(pos_ist[[i]])){
      incid[c(pos_ist[[i]][j,]),j]<--1
    }
    W_mats<-append(W_mats,list(incid))
    incid<-incidence_matrix_transpose
  }
  ######################
  in_H_negative_pos<-NULL
  for(i in 1:ncol(H)){
    in_H_negative_pos<-rbind(in_H_negative_pos, which(H[,i]==-1))
  }
  ####################
  position_mat_list<-list()
  for(j in 1:nrow(in_H_negative_pos)){
    pos_mat<-design[,c(in_H_negative_pos[j,])]
    position_mat_list<-append(position_mat_list,list(pos_mat))
  }
  ##################3
  all_w_mats<-list()
  for(i in 1:length(position_mat_list)){
    incid1<-incid
    for(j in 1:nrow(position_mat_list[[i]])){
      #for(k in 1:nrow(position_mat_list[[i]])){
      incid1[c(position_mat_list[[i]][j,]),j]<--1
    }
    all_w_mats<-append(all_w_mats,list(incid1))
  }
  ##################################
  W_mats<-all_w_mats
  interproductsum<-list()
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

#hadmardOCDs(mat)




