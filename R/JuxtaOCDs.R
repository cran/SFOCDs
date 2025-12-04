#' Juxtaposed Method for Optimal Covariate Designs (OCDs)
#'
#' @param design Input a design in matrix format and block size k such that (k+1) is a prime number.
#'@description
#'Consider L matrix, construct resolvable sets by grouping columns into pairs that have the same ordered set of elements. For each pair, arrange the two column-sets horizontally (2(s-1)) and change the signs of any two sets. This new setup form the P_i matrix of order 2(s-1) x 2. Then superimpose the first column of P_i onto N and that produce W_i, where first set of order v x b will be W_i^11 and other set below is W_i^21. Likewise use second column of P_i  to  get W_i^12 and W_i^22.  Repeat for every P_i to get collection of W_i's. The grand total of Hadamard product of all W_i^ij will be zero provided a foldover of any one of the W_i^ij is taken.
#'
#' @returns Generates W matrices and Inter product sums of W matrices.
#' @references
#' Das, K., N. K. Mandal, and B. K. Sinha. (2003) <https://doi.org/10.1016/S0378-3758(02)00151-9>. Optimal experimental designs for models with covariates. Journal of Statistical Planning and Inference 115(1): 273-285.
#'
#' Bansal, N., and D. K. Garg. (2022)<https://doi.org/10.1007/s42519-022-00244-0>. Optimum covariate designs for three associate PBIB designs. Journal of Statistical Theory and Practice 16(3): 1-15.
#'
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
#'
#' @examples
#' library(SFOCDs)
#' mat1<-matrix(c(
#' 1,2,3,4,5,6,
#' 7,8,9,10,11,1,
#' 12,13,14,15,2,7,
#' 16,17,18,3,8,12,
#' 19,20,4,9,13,16,
#' 21,5,10,14,17,19,
#' 6,11,15,18,20,21),nrow=7,byrow=TRUE)
#' JuxtaOCDs(mat1)
JuxtaOCDs<-function(design){
  design<-as.matrix(design)
  block <- nrow(design)
  trt <- max(design)
  mat=design
  incidence_mat<-function(mat){
    mat<-as.matrix(mat)
    incidence_matrix <- matrix(0, nrow = nrow(mat), ncol = max(c(mat)))
    for (i in 1:nrow(mat)) {
      #for (j in 1:max(c(mat))) {
      incidence_matrix[i, c(mat[i, ])] <- 1
      #}
    }
    return(incidence_matrix)
  }
  incidence_matrix_transpose <- t(incidence_mat(design))
  inc<-incidence_matrix_transpose
  #write.table(inc,"clipboard",sep="\t")

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
  k=ncol(design)
  mols<-MOLS(k+1)
  del_col_mols<-lapply(mols,function(mat)mat[,1])
  L<-do.call(cbind,del_col_mols)
  L
  L<-L[-1,]
  L<-L-1
  #####pmatrix
  same_pairs <- function(A, B) {
    sort_rows <- function(mat) {
      t(apply(mat, 1, sort))
    }

    A<-sort_rows(A)
    B<-sort_rows(B)
    # Convert rows to character strings (unique row identifiers)
    A_rows <- apply(A, 1, paste, collapse = "_")
    B_rows <- apply(B, 1, paste, collapse = "_")

    # Compare sorted row sets
    identical(sort(A_rows), sort(B_rows))
  }
  ###difference vector
  diff<-function(mat){
    A<-mat[,2]-mat[,1]
    return(A)
  }
  ##first pair of column fix for P1
  # first_col<-as.matrix(L[,1])
  # second_col<-as.matrix(L[,2])
  first_fixed_pair<-L[,1:2]
  res_cols1<-c(3:ncol(L))
  comb1<-t(combn(res_cols1,2))
  first_diff<-diff(first_fixed_pair)
  ############
  for(i in 1:nrow(comb1)){
    ref_mat1<-L[,c(comb1[i,])]
    if(same_pairs(first_fixed_pair,ref_mat1)==TRUE){
      if(identical(first_diff,diff(ref_mat1))==TRUE){
        break
      }
    }
  }

  ##second pair of column fix for P2
  # first_col<-as.matrix(L[,1])
  # third_col<-as.matrix(L[,3])
  second_fixed_pair<-L[,c(1,3)]
  res_cols2<-c(2,4:ncol(L))
  second_diff<-diff(second_fixed_pair)
  comb2<-t(combn(res_cols2,2))
  ############
  for(i in 1:nrow(comb2)){
    ref_mat2<-L[,c(comb2[i,])]
    if(same_pairs(second_fixed_pair,ref_mat2)==TRUE){
      if(identical(second_diff,diff(ref_mat2))==TRUE){
        break
      }
    }
  }
  #########P1 and P2 mats
  P1<-rbind(cbind(first_fixed_pair[,1],ref_mat1[,2]*-1),
            cbind(first_fixed_pair[,2]*-1,ref_mat1[,1]))
  P2<-rbind(cbind(second_fixed_pair[,1],ref_mat2[,2]*-1),
            cbind(second_fixed_pair[,2]*-1,ref_mat2[,1]))
  ############################
  ###Superimposing
  design
  inc
  ##for first case
  a<-first_fixed_pair[,1]
  b<-first_fixed_pair[,2]*-1
  d<-ref_mat1[,1]
  c<-ref_mat1[,2]*-1
  mat1<-cbind(a,b,c,d)
  ###superimposing in first case
  w1_list<-list()
  for(j in 1:ncol(mat1)){
    inc2<-inc
    for(i in 1:nrow(design)){
      inc2[c(design[i,]),i]<-mat1[,j]
    }
    w1_list<-append(w1_list,list(inc2))
  }
  ##for second case
  A<-second_fixed_pair[,1]
  B<-second_fixed_pair[,2]*-1
  D<-ref_mat2[,1]
  C<-ref_mat2[,2]*-1
  mat2<-cbind(A,B,C,D)
  ############
  ###superimposing in second case
  w2_list<-list()
  for(j in 1:ncol(mat2)){
    inc2<-inc
    for(i in 1:nrow(design)){
      inc2[c(design[i,]),i]<-mat2[,j]
    }
    w2_list<-append(w2_list,list(inc2))
  }
  #############Now W's final form
  W1<-cbind(rbind(w1_list[[1]],w1_list[[2]]),
            rbind(w1_list[[3]],w1_list[[4]]))
  W2<-cbind(rbind(w2_list[[1]],w2_list[[2]]),
            rbind(w2_list[[3]],w2_list[[4]]))

  ###############
  W_mats<-list(W1,W2)
  #apply(W1,2,sum)
  #interproductsum<-list()
  #indices<-list()
  w11_1<-w1_list[[1]]
  w12_1<-w1_list[[3]]
  w21_1<-w1_list[[2]]
  w22_1<-w1_list[[4]]
  ############
  w11_2<-w2_list[[1]]
  w12_2<-w2_list[[3]]
  w21_2<-w2_list[[2]]
  w22_2<-w2_list[[4]]
  ###################
  interproductsum<-w11_1*w11_2+w12_1*-w12_2+w21_1*w21_2+w22_1*-w22_2
  interproductsum<-apply(interproductsum,1,sum)
  new_list1 <- lapply(W_mats, function(W) {
    list(
      w = W,
      row_sum = rowSums(W),
      col_sum = colSums(W)
    )
  })

  #new_list<-split(interproductsum,rep(1:(length(interproductsum)/2),each=2))
  lm<-list("W matrices",new_list1,"Interproduct Sums",interproductsum)
  return(lm)
}
