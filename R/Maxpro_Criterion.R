#' Maxpro Criterion
#'
#' @param design Input a design in matrix format
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
return(value)
}
