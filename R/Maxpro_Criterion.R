#' Maxpro Criterion for Treatment vs Position
#'
#' @param design Input a design in matrix format
#'@description
#'User input should be the original design and this function
#'automatically convert the design in treatment vs position 2D array and then
#'print the Maxpro Criterion value.
#'
#' @returns Maxpro criterion value for a given design
#' @export
#'
#' @examples
#' library(SFOCDs)
#' mat<-matrix(c(
#' 1,2,3,
#' 2,1,4,
#' 3,4,1,
#' 4,3,2),nrow=4,byrow=TRUE)
#' Maxpro_Criterion(mat)

Maxpro_Criterion<-function(design){
  d<-as.matrix(design)
  #Convert design to proper format
  convert_des<-function(design){
    mat<-as.matrix(design)
    final_mat<-NULL
    for(i in 1:nrow(mat)){
      for(j in 1:ncol(mat)){
        final_mat=rbind(final_mat,c(mat[i,j],j,i))
      }
    }
    return(final_mat[,-ncol(final_mat)])
  }
  ##########

  #####Maxpro_criterion

  MP<-function(design){
    #design<-read.table("clipboard",sep='\t')
    d<-as.matrix(design)
    v=max(d)
    k=ncol(d)
    sumsq<-NULL
    for(i in 1:(nrow(d)-1)){
      j=i+1
      while(j<=nrow(d)){
        sumsq=c(sumsq,(1/prod((abs(d[i,1]-d[j,1])+(1/v))^2,((abs(d[i,2]-d[j,2])+(1/k)))^2)))
        j=j+1
      }
    }
    value<-(sum(sumsq)/choose(nrow(d),2))^(1/2)
    value
    return(value)
  }
  ######input the design for which maxpro want to find
  des<-convert_des(d)

  #########Find Maxpro value of converted design
  return(MP(des))
}
