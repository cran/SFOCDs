#' Reshuffle Design
#'
#' @param design Input any base design
#'@param type Input 1 or 2. By default set 1. If type =2 then reshuffled for JuxtaOCDs.
#' @returns Reshuffled the base design such that each column has minimum possible repeatations of symbols.
#' @export
#'
#' @examples
#' library(SFOCDs)
#' mat<-matrix(c(1,2,3,4,1,2,4,5,1,3,4,5,2,3,4,5),nrow=4,byrow=TRUE)
#' reshuffle_des(mat,type=1)

reshuffle_des<-function(design,type=1){
  design=as.matrix(design)
  if(type==1){
    reshuffled_mat <- function(design) {

      b <- nrow(design)
      k <- ncol(design)

      # replication of each treatment
      trt <- as.vector(design)
      r <- table(trt)[1]   # assuming equal replication

      beta <- r / k

      result <- matrix(NA, b, k)

      fill_row <- function(i) {
        if (i > b) return(TRUE)

        perms <- gtools::permutations(k, k, design[i, ])

        for (p in 1:nrow(perms)) {

          candidate <- perms[p, ]
          valid <- TRUE

          for (j in 1:k) {

            count <- sum(result[, j] == candidate[j], na.rm = TRUE)

            if (count >= beta) {
              valid <- FALSE
              break
            }
          }

          if (valid) {
            result[i, ] <<- candidate

            if (fill_row(i + 1)) return(TRUE)

            result[i, ] <<- NA
          }
        }

        return(FALSE)
      }

      success <- fill_row(1)

      if (!success) stop("No valid arrangement exists")

      return(result)
    }
  }else{
reshuffled_mat <- function(design) {
  k <- ncol(design)
  target <- k + 1

  new_design <- matrix(NA, nrow=nrow(design), ncol=k)
  pos_list <- list()

  # First row fixed
  new_design[1, ] <- design[1, ]

  for (j in 1:k) {
    pos_list[[as.character(design[1, j])]] <- j
  }

  # Process remaining rows
  for (i in 2:nrow(design)) {
    row <- design[i, ]
    new_row <- rep(NA, k)
    used_pos <- rep(FALSE, k)

    # Step 1: try complementary placement
    for (j in 1:k) {
      trt <- row[j]
      trt_char <- as.character(trt)

      if (!is.null(pos_list[[trt_char]])) {
        pos <- target - pos_list[[trt_char]]

        if (!used_pos[pos]) {
          new_row[pos] <- trt
          used_pos[pos] <- TRUE
        }
      }
    }

    # Step 2: fill remaining safely
    for (j in 1:k) {
      trt <- row[j]
      if (!(trt %in% new_row)) {
        free_pos <- which(!used_pos)[1]
        new_row[free_pos] <- trt
        used_pos[free_pos] <- TRUE

        # store position if first time
        trt_char <- as.character(trt)
        if (is.null(pos_list[[trt_char]])) {
          pos_list[[trt_char]] <- free_pos
        }
      }
    }

    new_design[i, ] <- new_row
  }

  return(new_design)
}
  }
  return(reshuffled_mat(design))
}
