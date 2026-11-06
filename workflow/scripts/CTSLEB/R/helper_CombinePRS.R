#' Helper function for PRSscore()
#' @description
#' Combines PRS
#' Not exported
#' @param scores  description
#' @param pthres description
#' @param prs_p_other_  description
#' @return 'prs_mat' global variable. Matrix of combined PRSs with clumping and
#' pthres as column name
#'
helper_CombinePRS <- function(scores,
                              pthres,
                              prs_p_other_)
  {

  print("combining PRS results ...")
  prs_list <- list()
  temp <- 1
  names <- colnames(scores[3:ncol(scores)])
  for(k1 in 1:length(pthres)){
    for(k2 in 1:length(pthres)){
      prs_file <- paste0(prs_p_other_,k1,".p_tar_",k2,".sscore")
      #print(paste0("reading :", prs_file))
      prs_temp <- as.data.frame(fread(prs_file))
      score_cols <- grep("^SCORE[0-9]+_SUM$", colnames(prs_temp))
      if (length(score_cols) == 0L) {
        score_cols <- 5:ncol(prs_temp)
        warning("No SCORE*_SUM columns found; falling back to columns 5:ncol")
      }
      if (length(score_cols) != length(names)) {
        stop(sprintf("Expected %d score columns, found %d",
                     length(names), length(score_cols)))
      }
      prs_list[[temp]] <- prs_temp[, score_cols, drop = FALSE]
      colnames(prs_list[[temp]]) <- paste0(names,"_",
                                           "p_other_",
                                           pthres[k1],
                                           "_p_tar_",
                                           pthres[k2])
      temp <- temp + 1
    }
  }
  prs_mat <- as.data.frame(cbind(prs_temp[,1:2],bind_cols(prs_list)))
  print("prs_mat complete... ")
  return(prs_mat)
}
