
#' Treatment Position Vs Treatment Scatter Plot
#'
#' @param design Input a design (base or reshuffled) in matrix format
#'@description
#'The function will give the scatter plot showing the spread offered by design points in the experimental space. The x axis of the plot represent treatments and y axis the position of treatments in each block. Different colors in the dot represent the blocks.
#'@importFrom graphics axis
#'@importFrom utils combn
#' @returns Generates scatter plot of treatment position Vs treatment
#' @export
#'
#' @examples
#' library(SFOCDs)
#' mat<-matrix(c(  1,  4,  2,  5,
#' 2,  5,  3,  6,
#' 3,  6,  1,  4,
#' 4,  1,  5,  2,
#' 5,  2,  6,  3,
#' 6,  3,  4,  1),nrow=6,byrow=TRUE)
#' MaxDot(mat)
MaxDot<-function(design){

  mat<-as.matrix(design)
  #mat<-rearrange_backtrack(mat)
  final_mat<-NULL
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      final_mat=rbind(final_mat,c(mat[i,j],j,i))
    }
  }
  final_mat
  colnames(final_mat)<-c("Treatment","Position","Block")
  Array<-final_mat

  #plot

  aray<-Array[order(Array[,3]),]
  # Define colors for each category
  b=nrow(mat)
  color_the_blocks<-function(matrix,b){
    matrix<-as.matrix(matrix)
    k<-nrow(matrix)/b
    colors <- c("red", "blue", "green", "purple", "orange","violet","wheat","yellow",
                "turquoise","skyblue","salmon","black","pink","peru","palegreen","navy","hotpink","gold",
                "firebrick","azure")
    seq<-colors[1:b]
    new_seq<-rep(seq,each=k)
    pch_seq<-rep(1:b,each=k)
    # Create a scatter plot with different colors for different categories
    plot(matrix[,2], matrix[,1], xlab="Position", ylab="Treatment", col = new_seq,xaxt="n", pch = 16,cex=2,bg=, main = "Scatter Plot ")
    axis(1,at=seq(1,b,by=1),labels=seq(1,b,by=1))
  }
  color_the_blocks(aray,b)

}



