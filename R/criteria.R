
# PS: 
# i <- which.min(mat_o[, "objfun"]), which is the first index that 
# gets the the same minimum


#' Hausdorff distance between vectors
#' 
#' Calculate the hausdorff distance between vectors
#' 
#' @param B Vector, B_{m} the true change-points locations
#' @param eB Vector, eB_{n} the estimated change-points locations.
#'   For example \code{eB = c(20, 40, 60)}, without 0 and nX included.
#' @param nX Interger, the length of X matrix
#' @return A list of OE (over estimated) and UE(under estimated)
hausdorff_dist_vector <- function(B, eB, nX) {
  if (length(eB) == 0) {
    return(list(OE=0, UE=nX))
  }
  m <- length(B)
  n <- length(eB)
  result <- matrix(0, m, n)
  for (i in 1:m) {
    for (j in 1:n) { result[i, j] = abs(B[i] - eB[j])}
  }
  OE <- max(apply(result, 1, min))
  UE <- max(apply(result, 2, min))
  return(list(OE = OE, UE = UE))
}


#' The distribution of change-point num
#' 
#' The criterion of the distribution of change-point num
#' 
#' @param K_mat Matrix, the estimated change-point matrix, the results of OP, PELT, DP.
#'   A n_sim * n_algorithms matrix, whose element is the value of K_hat by cv. 
#' @param Ktrue Vector, the true change-point vector
#' @return A matrix, crit_mat
criterion_cp_num <- function(K_mat, Ktrue) {
  n <- nrow(K_mat)
  cnames <- c("-3", "-2", "-1", "0", "1", "2", "3", "mean", "SD", "MSE")
  crit_mat <- matrix(0, ncol(K_mat), 10, dimnames = list(colnames(K_mat), cnames))
  for (kname in colnames(K_mat)) {  # kname = colnames(K_mat)[1]
    K_hat <- K_mat[, kname]
    crit_mat[kname, "-3"] <- sum(K_hat - Ktrue <= -3) / n * 100
    crit_mat[kname, "3"] <- sum(K_hat - Ktrue >= 3) / n * 100
    v <- c(-2, -1, 0, 1, 2)
    for (i in v) {
      crit_mat[kname, paste0(i)] = sum(K_hat - Ktrue == i) / n * 100
    }
    x0 <- K_hat - Ktrue
    crit_mat[kname, "mean"] <- mean(K_hat - Ktrue)
    crit_mat[kname, "SD"] <- sd(K_hat - Ktrue)
    crit_mat[kname, "MSE"] <- mean((x0 - mean(x0))^2)
  }
  return(crit_mat)
}


#' Criterion of Hausdorff distance
#' 
#' Calculate the criterion of Hausdorff distance between vectors.
#' 
#' @param hausdorff_OE vector, overestimated
#' @param hausdorff_UE vector, underestimated
#' @return A matrix, crit_mat
criterion_hausdorff <- function(hausdorff_OE, hausdorff_UE) {
  # cnames <- c( "OE_mean", "OE_SD", "OE_MSE", "UE_mean", "UE_SD", "UE_MSE")
  cnames <- c()
  eg_nm <- expand.grid(c("mean", "SD", "MSE"), c("OE", "UE"))
  for (i in 1:nrow(eg_nm)) {
    cnames <- c(cnames, paste(eg_nm[i, 2], eg_nm[i, 1], sep="_"))
  }
  crit_mat <- matrix(0, 1, 6)
  colnames(crit_mat) <- cnames
  crit_mat[1, "OE_mean"] <- mean(hausdorff_OE)
  crit_mat[1, "OE_SD"] <- sd(hausdorff_OE)
  crit_mat[1, "OE_MSE"] <- mean((hausdorff_OE - mean(hausdorff_OE))^2)
  crit_mat[1, "UE_mean"] <- mean(hausdorff_UE)
  crit_mat[1, "UE_SD"] <- sd(hausdorff_UE)
  crit_mat[1, "UE_MSE"] <- mean((hausdorff_UE - mean(hausdorff_UE))^2)
  return(crit_mat)
}


errMatTrim <- function(errMat) {
  ### To check if "pen" is in the columns of errMat
  if ("pen" %in% colnames(errMat)) {
    ### which means it is the result of the OP and PELT.
    ### drop duplicated "K_hat" and keep the one with the biggest penalty "pen".
    errMat <- errMat[order(errMat[, "K_hat"], -errMat[, "pen"]), ]
    errMat <- errMat[!duplicated(errMat[, "K_hat"]), ,drop=FALSE]
    ### for the special case with only one item, set drop=FALSE.
    if ("K_diff" %in% colnames(errMat)) {
      h_names <- c("pen", "objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
    } else {
      h_names <- c("pen", "objfit", "objfun", "K_hat")
    }
  } else {
    errMat <- errMat[order(errMat[, "K_hat"]), ]
    if ("K_diff" %in% colnames(errMat)) {
      h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
    } else {
      h_names <- c("objfit", "objfun", "K_hat")
    }
  }
  
  if ("K_diff" %in% colnames(errMat)) {
    errMat <- errMat[abs(errMat[, "K_diff"]) < 20, ,drop=FALSE]
  }
  ### remove the columns with all zeroes after trimming!
  errMat <- errMat[, colSums(errMat != 0, na.rm = TRUE) > 0]
  ## errMat2 <- errMat[, !colSums(errMat != 0, na.rm = TRUE) < 1e-4]
  ### rename the others as from tau_1, to tau_L
  k <- length(h_names)
  m <- dim(errMat)[2]
  colnames(errMat)[(k+1):m] <- paste("tau", c(1:(m-k)), sep="_")
  return(errMat)
}


CVMatConcat <- function(mat_o, mat_e){
  ## get the union of results of cv_o and cv_e.
  ## If mat_o or mat_e is not Trimmed, then Trim it!
  ## mat_o <- errMatTrim(mat_o)
  ## mat_e <- errMatTrim(mat_e)
  # if ("pen" %in% colnames(mat_o)) {
  #   if ("OE" %in% colnames(mat_o)) {
  #     h_names <- c("pen", "objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  #     by_names <- c("K_hat", "K_diff")
  #   } else {
  #     h_names <- c("pen", "objfit", "objfun", "K_hat")
  #     by_names <- c("K_hat")
  #   }
  # } else {
  #   if ("OE" %in% colnames(mat_o)) {
  #     h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  #     by_names <- c("K_hat", "K_diff")
  #   } else {
  #     h_names <- c("objfit", "objfun", "K_hat")
  #     by_names <- c("K_hat")
  #   }
  # }
  # df_o <- data.frame(mat_o[, h_names, drop=FALSE])
  # df_e <- data.frame(mat_e[, h_names, drop=FALSE])
  
  if ("OE" %in% colnames(mat_o)) {
    by_names <- c("K_hat", "K_diff")
  } else {
    by_names <- c("K_hat")
  }
  df_o <- data.frame(mat_o)
  df_e <- data.frame(mat_e)
  
  ### we do not combine the estimated change-points by sum or mean function
  ### it is meaningless to get the mean of (20, 40, 60) and (18, 46, 70)
  ### So we only use the h_names to combine all the columns including tau's
  # df_union <- merge(df_o, df_e, by=by_names, all=TRUE)  # union
  # df_inner <- merge(df_o, df_e, by=by_names, all=FALSE)  # inter
  df <- merge(df_o, df_e, by=by_names, all=TRUE)  # union
  
  ### concatenate type is "intersection", "union" respectively
  ### inter--means both of mat_o and mat_e have values;
  ### union--means any of mat_o and mat_e has value.
  tau_names <- colnames(mat_o)[grep("^tau", colnames(mat_o))]
  if ("OE" %in% colnames(mat_o)) {
    tmp_names <- c(c("objfit", "objfun", "OE", "UE"), tau_names)
  } else {
    tmp_names <- c(c("objfit", "objfun"), tau_names)
  }
  col_names <- colnames(df)
  drop_cols <- c()
  for (name in tmp_names) {  # name = tmp_names[1]
    ## paste0(name, c(".x$", ".y$"), collapse="|")
    sel_cols <- grep(
      paste0(name, c(".x$", ".y$"), collapse="|"), col_names, value=TRUE, invert = FALSE)
    if (length(sel_cols) == 0) {next}
    df[, paste0(name, "_inter")] <- apply(df[, sel_cols], 1, mean, na.rm=FALSE)
    df[, paste0(name, "_union")] <- apply(df[, sel_cols], 1, mean, na.rm=TRUE)
    drop_cols <- c(drop_cols, sel_cols)
  }
  
  # sel_cols <- grep(
  #   paste(drop_cols, collapse="|"), colnames(df), value=TRUE, invert = TRUE)
  # df2 <- df[, sel_cols]
  df2 <- df[ , -which(colnames(df) %in% drop_cols)]
  
  ### It is the same as the following but saved some energy
  # df["UE_inter"] <- apply(df[, c("UE.x", "UE.y")], 1, mean, na.rm=FALSE)
  # df["UE_union"] <- apply(df[, c("UE.x", "UE.y")], 1, mean, na.rm=TRUE)
  # df2 <- df[,-c(grep(".x$", colnames(df)), grep(".y$", colnames(df)))]
  
  ### for the intersection data
  df_inter <- df2[, -grep("_union", colnames(df2))]  # remove
  ### df_inter <- df_inter[complete.cases(df_inter), ,drop=FALSE]
  df_inter <- df_inter[order(df_inter[, "K_hat"]), ]
  col_old <- colnames(df_inter)  # removing the suffixes
  col_new <- gsub(pattern = "_inter", replacement = "", x = col_old)
  colnames(df_inter) <- col_new
  df_inter <- as.matrix(df_inter)
  ### for the union data
  df_union <- df2[, -grep("_inter", colnames(df2))]  # remove
  ### df_union <- df_union[complete.cases(df_union), ,drop=FALSE]
  df_union <- df_union[order(df_union[, "K_hat"]), ]
  col_old <- colnames(df_union)  # removing the suffixes
  col_new <- gsub(pattern = "_union", replacement = "", x = col_old)
  colnames(df_union) <- col_new
  df_union <- as.matrix(df_union)
  return(list(df_inter = df_inter, df_union = df_union))
}


#' Function converts errMat to optimal Khat result
#' 
#' Function converts errMat to optimal Khat result
#' 
#' @param res A list of mat_o and mat_e, both are error matrices.
#' @return A list of opt_cploc_e, opt_cploc_o, opt_skou, opt_Khat.
#'   opt_cploc_o, opt_cploc_e, are both arrays, the optimal location;
#'   opt_skou, (sick_or_not, K_hat, K_diff, OE, UE)
#'   opt_Khat, a list of integers, the optimal change-point number, 
#'   from (O, E, inter, union) data.
CVMat_optimal_Khat <- function(res){
  ## which() returns the first index that it touches the minimum
  mat_o <- res$mat_o
  mat_e <- res$mat_e
  mat_o <- errMatTrim(mat_o)
  mat_e <- errMatTrim(mat_e)
  
  df_list <- CVMatConcat(mat_o, mat_e)
  df_inter <- df_list$df_inter
  df_union <- df_list$df_union
  
  if ("OE" %in% colnames(mat_e)) {
    s_names <- c("K_hat", "K_diff", "OE", "UE")
  } else if ("K_diff" %in% colnames(mat_e)) {
    s_names <- c("K_hat", "K_diff")
  } else {
    s_names <- c("K_hat")
  }
  
  ### vo is an array of (sick_or_not, K_hat, OE, UE)
  ## For mat_o
  flag <- any(diff(mat_o[, "objfit"]) > 0)
  idx <- which.min(mat_o[, "objfun"])
  KE_o <- mat_o[idx, s_names, drop=FALSE]
  vo <- c(sick=flag, K_hat=KE_o)
  opt_o <- mat_o[idx, ,drop=FALSE]
  
  flag <- any(diff(mat_e[, "objfit"]) > 0)
  idx <- which.min(mat_e[, "objfun"])
  KE_e <- mat_e[idx, s_names, drop=FALSE]
  ve <- c(sick=flag, K_hat=KE_e)
  opt_e <- mat_e[idx, ,drop=FALSE]
  
  flag <- any(diff(df_union[, "objfit"]) > 0)
  idx <- which.min(df_union[, "objfun"])
  KE_u <- df_union[idx, s_names, drop=FALSE]
  vu <- c(sick=flag, K_hat=KE_u)
  
  flag <- any(diff(df_inter[, "objfit"]) > 0)
  idx <- which.min(df_inter[, "objfun"])
  KE_i <- df_inter[idx, s_names, drop=FALSE]
  vi <- c(sick=flag, K_hat=KE_i)
  
  ### To check if "pen" is in the columns of errMat
  ## length(h_names) is 6
  ## h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  if ("pen" %in% colnames(mat_o)) {
    mo <- 3 + length(s_names) + opt_o[1, "K_hat"]
    me <- 3 + length(s_names) + opt_e[1, "K_hat"]
  } else {
    mo <- 2 + length(s_names) + opt_o[1, "K_hat"]
    me <- 2 + length(s_names) + opt_e[1, "K_hat"]
  }
  
  KHAT <- unname(c(vo['K_hat'], ve['K_hat'], vi['K_hat'], vu['K_hat']))
  SKOU <- unname(c(vo, ve, vi, vu))  # (sick_or_not, K_hat, K_diff, OE, UE)
  opt_list <- list(
    opt_cploc_o = opt_o[, 1:mo], 
    opt_cploc_e = opt_e[, 1:me], 
    opt_skou = SKOU, 
    opt_Khat = KHAT)
  return(opt_list)
}



#' Optimal change-point number using high-dimension BIC 
#' 
#' Function gives the optimal change-point number using High-dimension BIC 
#' 
#' @param errMat The error matrix
#' @return A list of optimal location and optimal number.
#'   opt_cploc, array, the optimal location;
#'   opt_Khat, an integer, the optimal change-point number.
errMat_optimal_Khat <- function(errMat){
  ### Usually we feed in err_mat obtained from hBIC, or CVMatConcat
  ### We only need columns of "objfit", "objfun", "K_hat", "K_diff", "OE", "UE".
  ### if the value of flag is TRUE, then the fitting process is wrong.
  errMat <- errMat[order(errMat[, "K_hat"]), ]  # increasing "K_hat"
  flag <- any(diff(errMat[, "objfit"]) > 0)
  i <- which.min(errMat[, "objfun"])
  if ("OE" %in% colnames(errMat)) {
    h_names <- c("K_hat", "K_diff", "OE", "UE")
  } else if ("K_diff" %in% colnames(errMat)) {
    h_names <- c("K_hat", "K_diff")
  } else {
    h_names <- c("K_hat")
  }
  KE <- errMat[i, h_names]
  SKOU <- unname(c(sick=flag, KE))
  opt <- errMat[i, ,drop=FALSE]
  ### To check if "pen" is in the columns of errMat
  ## h_names <- c("pen", "objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  ## Basically, length(h_names) is 6.
  if ("pen" %in% colnames(errMat)) {
    m <- 3 + length(h_names)
    ### which means the errMat includes the estimated tau_hat.
    if (ncol(errMat) > m) {
      m <- 3 + length(h_names) + opt[1, "K_hat"]
    }
  } else {
    m <- 2 + length(h_names)
    if (ncol(errMat) > m) {
      m <- 2 + length(h_names) + opt[1, "K_hat"]
    }
  }
  
  opt_list <- list(
    opt_cploc = opt[, 1:m], opt_skou = SKOU, opt_Khat = errMat[i, "K_hat"])
  return(opt_list)
}


#' Optimal change-point number using CV
#' 
#' Function gives the optimal change-point number using CV
#' 
#' It is the same as \code{\link{CVMat_optimal_Khat}}.
#' 
#' @param res A list of mat_o and mat_e, both are error matrices.
#' @return A list of opt_cploc_e, opt_cploc_o, opt_skou, opt_Khat.
#'   opt_cploc_o, opt_cploc_e, are both arrays, the optimal location;
#'   opt_skou, (sick_or_not, K_hat, K_diff, OE, UE)
#'   opt_Khat, a list of integers, the optimal change-point number, 
#'    from (O, E, inter, union) data.
hCV_optimal_Khat <- function(res) {
  # res is the result of the copss_oppelt_cv, or copss_dp_cv
  # res = op_res; # res = pelt_res; # res = dp_res;
  mat_o <- res$mat_o
  mat_e <- res$mat_e
  if ("pen" %in% colnames(mat_o)) {
    ### which means it is the result of the OP and PELT.
    ### drop duplicated "K_hat" and keep the one with the biggest "pen"
    mat_o <- mat_o[order(mat_o[, "K_hat"], -mat_o[, "pen"]), ]
    mat_e <- mat_e[order(mat_e[, "K_hat"], -mat_e[, "pen"]), ]
    mat_o <- mat_o[!duplicated(mat_o[, "K_hat"]), ,drop=FALSE]
    mat_e <- mat_e[!duplicated(mat_e[, "K_hat"]), ,drop=FALSE]
    ### for the special case of only one item, set drop=FALSE.
  } else {
    mat_o <- mat_o[order(mat_o[, "K_hat"]), ]
    mat_e <- mat_e[order(mat_e[, "K_hat"]), ]
  }
  mat_e <- mat_e[abs(mat_e[, "K_diff"]) < 20, ,drop=FALSE]
  mat_o <- mat_o[abs(mat_o[, "K_diff"]) < 20, ,drop=FALSE]
  
  # plot(mat_o[, "K_hat"], mat_o[, "objfun"])
  # plot(mat_e[, "K_hat"], mat_e[, "objfun"])
  # plot(mat_e[, "K_hat"], mat_e[, "objfit"])
  
  ### The flags are used to check if the fitting function is decreasing.
  flag_o <- any(diff(mat_o[, "objfit"]) > 0)
  ### which is the first index that gets the the same minimum
  io <- which.min(mat_o[, "objfun"])
  KE_o <- mat_o[io, c("K_hat", "K_diff", "OE", "UE")]
  vo <- c(sick=flag_o, KE_o)
  opt_o <- mat_o[io, ,drop=FALSE]  # optimal of mat_o
  # mat_o <- mat_o[order(mat_o[, "K_hat"], decreasing = FALSE), ,drop=FALSE]
  # mat_o <- mat_o[order(mat_o[, "objfun"], decreasing = FALSE), ,drop=FALSE]
  
  flag_e <- any(diff(mat_e[, "objfit"]) > 0)
  ie <- which.min(mat_e[, "objfun"])
  KE_e <- mat_e[ie, c("K_hat", "K_diff", "OE", "UE")]
  ve <- c(sick=flag_e, KE_e)
  opt_e <- mat_e[ie, ,drop=FALSE]
  
  ### get the union of results of cv_o and cv_e.
  if ("pen" %in% colnames(mat_o)) {
    h_names <- c("pen", "objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
    by_names <- c("pen", "K_hat")
  } else{
    h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
    by_names <- c("K_hat")
  }
  df_o <- data.frame(mat_o[, h_names, drop=FALSE])
  df_e <- data.frame(mat_e[, h_names, drop=FALSE])
  
  # df_union <- merge(df_o, df_e, by=by_names, all=TRUE)  # union
  # df_inner <- merge(df_o, df_e, by=by_names, all=FALSE)  # inter
  df <- merge(df_o, df_e, by=by_names, all=TRUE)  # union
  
  ### concatenate type is "intersection", "union" respectively
  ### inter--means both of mat_o and mat_e have values;
  ### union--means any of mat_o and mat_e has value.
  tmp_names <- c("objfit", "objfun", "K_diff", "OE", "UE")
  col_names <- colnames(df)
  drop_cols <- c()
  for (name in tmp_names) {  # name = tmp_names[1]
    sel_cols <- grep(name, col_names, value=TRUE, invert = FALSE)
    if (length(sel_cols) == 0) {next}
    df[, paste0(name, "_inter")] <- apply(df[, sel_cols], 1, mean, na.rm=FALSE)
    df[, paste0(name, "_union")] <- apply(df[, sel_cols], 1, mean, na.rm=TRUE)
    drop_cols <- c(drop_cols, sel_cols)
  }
  # sel_cols <- grep(
  #   paste(drop_cols, collapse="|"), colnames(df), value=TRUE, invert = TRUE)
  # df2 <- df[, sel_cols]
  df2 <- df[ , -which(colnames(df) %in% drop_cols)]
  
  # df["objfit_inter"] <- apply(df[, c("objfit.x", "objfit.y")], 1, mean, na.rm=FALSE)
  # df["objfit_union"] <- apply(df[, c("objfit.x", "objfit.y")], 1, mean, na.rm=TRUE)
  # df2 <- df[,-c(grep(".x$", colnames(df)), grep(".y$", colnames(df)))]
  
  df_inter <- df2[, -grep("_union", colnames(df2))]
  df_inter <- df_inter[complete.cases(df_inter), ,drop=FALSE]
  df_inter <- df_inter[order(df_inter[, "K_hat"]), ]
  ## rename the columns 
  col_old <- colnames(df_inter)  # removing the suffixes
  col_new <- gsub(pattern = "_inter", replacement = "", x = col_old)
  colnames(df_inter) <- col_new
  df_inter <- data.matrix(df_inter)
  
  flag_i <- any(diff(df_inter[, "objfit"]) > 0)
  ii <- which.min(df_inter[, "objfun"])
  KE_i <- df_inter[ii, c("K_hat", "K_diff", "OE", "UE")]
  vi <- c(sick=flag_i, KE_i)
  ### for the union data
  df_union <- df2[, -grep("_inter", colnames(df2))]
  df_union <- df_union[complete.cases(df_union), ,drop=FALSE]
  df_union <- df_union[order(df_union[, "K_hat"]), ]
  ## rename the columns
  
  col_old <- colnames(df_union)  # removing the suffixes
  col_new <- gsub(pattern = "_union", replacement = "", x = col_old)
  colnames(df_union) <- col_new
  df_union <- data.matrix(df_union)
  
  flag_u <- any(diff(df_union[, "objfit"]) > 0)
  uu <- which.min(df_union[, "objfun"])
  KE_u <- df_union[uu, c("K_hat", "K_diff", "OE", "UE")]
  vu <- c(sick=flag_u, KE_u)
  
  ### To check if "pen" is in the columns of errMat
  ## length(h_names) is 6
  ## h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  if ("pen" %in% colnames(mat_o)) {
    mo <- 7 + opt_o[1, "K_hat"]
    me <- 7 + opt_e[1, "K_hat"]
  } else {
    mo <- 6 + opt_o[1, "K_hat"]
    me <- 6 + opt_e[1, "K_hat"]
  }
  ### vo is an array of (sick_or_not, K_hat, OE, UE)
  KHAT <- unname(c(vo['K_hat'], ve['K_hat'], vi['K_hat'], vu['K_hat']))
  SKOU <- unname(c(vo, ve, vi, vu))  # vo is (sick, K_hat, K_diff, OE, UE)
  opt_list <- list(
    opt_cploc_e = opt_o[, 1:mo], 
    opt_cploc_o = opt_e[, 1:me], 
    opt_skou = SKOU, 
    opt_Khat = KHAT)
  return(opt_list)
}













