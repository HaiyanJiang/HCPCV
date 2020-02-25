
# tau2n The true locations of change-points, which can be the ratio of 
#  \eqn{\tau / n} or the \eqn{\tau}'s. e.g, \code{tau2n=c(0.2, 0.4, 0.8)}, 
# or \code{tau2n = c(20, 40, 80)}, without 0, 1 or 0, n included.

copss_oppelt_hausdorff <- function(
    X, W, tau2n, cp_alg = cp_hdmean_op, penaltyform = "CV", normalized = TRUE){
  ## The copss method with OP or PELT algorithm to find the optimal number of 
  ## change-points and to give the hausdorff distance. 
  ## The change-point detection algorithm can be OP or PELT.
  ## cp_alg = cp_hdmean_pelt; normalized = TRUE
  Xo <- X[seq(1, nrow(X), 2), ]  # max(seq(1, nrow(X), 2))
  Xe <- X[seq(2, nrow(X), 2), ]  # max(seq(2, nrow(X), 2))
  if (!is.matrix(Xo)) {Xo = matrix(Xo, length(Xo), 1)}
  if (!is.matrix(Xe)) {Xe = matrix(Xe, length(Xe), 1)}
  
  n_o <- nrow(Xo)
  n_e <- nrow(Xe)
  B_o <- floor(tau2n * n_o)  # convert it to tau_true
  B_e <- floor(tau2n * n_e)
  
  #  GIVEN a GOOD sequence of penalty values
  alph_list <- c(0, 0.1, 0.3, 0.4, 0.5, 1, 1.1)
  c0_list <- seq(0.1, 10, by=0.1)
  eg <- expand.grid(alph_list, c0_list)
  n_pen <- nrow(eg)
  n <- nrow(X)
  p <- ncol(X)
  if (normalized) {traceR2 <- trace_R_square(X)} else {traceR2 <- p}
  pen_vec <- rep(0, n_pen)
  for (i in 1:n_pen) {
    alph <- eg[i, 1]
    c0 <- eg[i, 2]
    pen_vec[i] <- p + c0 * traceR2^(0.5) * (log(n/2)^(1.0 + alph))
    # pen_vec[i] <- calculate_penalty(X, normalized, penaltyform="CV", alph, c0)
  }
  pen_vec <- sort(pen_vec)  # increasing order, IMPORTANT.
  
  ### find some good ones larger than the first 30 pens.
  tmp <- pen_vec[30]
  Ko_list <- c()
  pen_poten <- c()  # potential penalty
  k <- 1
  while (TRUE) {
    pen_poten <- c(pen_poten, tmp)
    op_list <- cp_alg(X, W, pen = tmp)
    tau_hat <- op_list$cptau
    Lo <- length(tau_hat)
    # length(tau_hat); (Lo %in% Ko_list)
    if (!(Lo %in% Ko_list)) {Ko_list <- c(Ko_list, Lo)}
    tmp <- tmp * 1.1
    if ((tmp > pen_vec[n_pen]) | (Lo == 0)) {break}
    k <- k + 1
  }
  
  ### special penalties
  pen_special <- c()
  for (alph in c(0, 0.1, 0.3, 0.5)) {
    pen <- p + 2.5 * traceR2^(0.5) * (log(n/2)^(1.0 + alph))
    pen_special <- c(pen_special, pen)
  }
  pen_list <- c(pen_vec[1:30], pen_poten, pen_special)
  pen_list <- sort(pen_list)
  # plot(pen_list)
  epsilon <- 1e-4
  n_param <- length(pen_list)
  h_names <- c("pen", "objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  nc <- length(h_names)
  err_o <- matrix(0, n_param, nc)
  err_e <- matrix(0, n_param, nc)
  # err_e matrix, where Xo is used to estimate, Xe is to test;
  # err_o matrix, where Xe is used to estimate, Xo is to test. 
  colnames(err_o) <- h_names
  colnames(err_e) <- h_names
  
  # obj_oe: Xo is used to estimate, Xe is to test
  # obj_eo: Xe is used to estimate, Xo is to test.
  Ko_max <- 0
  Ke_max <- 0
  for (i in 1:n_param) {  # i = 1
    pen <- pen_list[i]
    ### for the data with odd index.
    ### op_list_o <- cp_alg(Xo, W, pen = pen)
    op_list_o <- cp_alg(Xo, W, pen = pen)
    ### For X data itself, when calculating cost_matrix, normalized nor not.
    ### For n=p=200, if normalized = FALSE, the algorithm does not work.
    ### For n=200, p=4, whatever normalized is, the algorithm works well.
    tau_hat_o <- op_list_o$cptau  # if there's no change-point, it also works
    cp_o <- c(0, tau_hat_o, n_o)
    obj_o <- cv_objfun(Xo, Xo, W, cp_o)
    obj_oe <- cv_objfun(Xo, Xe, W, cp_o)  ### fitting on Xo, testing on Xe
    if (obj_o < 1e-4) {obj_o <- 0;}  # IMPORTANT!!!
    if (obj_oe < 1e-4) {obj_oe <- 0;}  # IMPORTANT!!!
    Lo <- length(tau_hat_o)
    K_diff <- Lo - length(B_o)  # K_hat - K
    hv <- hausdorff_dist_vector(B_o, tau_hat_o, n_o)
    err_o[i, h_names] <- c(pen, obj_o, obj_oe, Lo, K_diff, hv$OE, hv$UE)
    n_diff <- max(Lo - Ko_max, 0)
    if (n_diff > 0) {
      err_o <- cbind(err_o, matrix(0, n_param, n_diff))
      Ko_max <- Ko_max + n_diff
    }
    if (Lo > 0) {err_o[i, (nc+1):(nc + Lo)] <- tau_hat_o}
    
    ### for the data with even index
    op_list_e <- cp_alg(Xe, W, pen = pen)
    tau_hat_e <- op_list_e$cptau
    cp_e <- c(0, tau_hat_e, n_e)
    obj_e <- cv_objfun(Xe, Xe, W, cp_e)
    obj_eo <- cv_objfun(Xe, Xo, W, cp_e)  # fitting on Xe testing on Xo
    if (obj_e < 1e-4) {obj_e <- 0;}  # IMPORTANT!!!
    if (obj_eo < 1e-4) {obj_eo <- 0;}  # IMPORTANT!!!
    Le <- length(tau_hat_e)
    K_diff <- Le - length(B_e)  # K_hat - K
    hv <- hausdorff_dist_vector(B_e, tau_hat_e, n_e)
    err_e[i, h_names] <- c(pen, obj_e, obj_eo, Le, K_diff, hv$OE, hv$UE)
    n_diff <- max(Le - Ke_max, 0)
    if (n_diff > 0) {
      err_e <- cbind(err_e, matrix(0, n_param, n_diff))
      Ke_max <- Ke_max + n_diff
    }
    if (Le > 0) {err_e[i, (nc+1):(nc + Le)] <- tau_hat_e}
  }
  # plot(err_e[, "K_hat"], err_e[, "objfun"])
  # plot(err_e[, "K_hat"], err_e[, "objfit"])
  # plot(err_o[, "K_hat"], err_o[, "objfun"])
  # plot(err_o[, "K_hat"], err_o[, "objfit"])
  res <- list(mat_e = err_e, mat_o = err_o)
  return(res)
}



copss_dp_hausdorff <- function(X, W, tau2n, Kmax = 15) {
  ## copss_dp_ErrMat_hausdorff is the same as copss_dp_hausdorff.
  # Xo, Xe, Xodd, Xeven
  Xo <- X[seq(1, nrow(X), 2), ]  # max(seq(1, nrow(X), 2))
  Xe <- X[seq(2, nrow(X), 2), ]  # max(seq(2, nrow(X), 2))
  if (!is.matrix(Xo)) {Xo = matrix(Xo, length(Xo), 1)}
  if (!is.matrix(Xe)) {Xe = matrix(Xe, length(Xe), 1)}
  
  n_o <- nrow(Xo)
  n_e <- nrow(Xe)
  B_o <- floor(tau2n * n_o)
  B_e <- floor(tau2n * n_e)
  
  Kmax <- min(Kmax, nrow(X) / 3)   ### tau_mat is Kmax*Kmax matrix
  tau_mat_o <- cp_hdmean_dp(Xo, W, Kmax = Kmax)
  tau_mat_e <- cp_hdmean_dp(Xe, W, Kmax = Kmax)
  ## for X data itself, when calculating cost_matrix, normalized nor not.
  
  t_names <- c()
  for (k in 1:Kmax) {t_names <- c(t_names, paste0("tau_", k))}
  h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  err_o <- matrix(NA, Kmax, length(h_names) + Kmax)
  err_e <- matrix(NA, Kmax, length(h_names) + Kmax)
  # err_o matrix, where Xo is used to estimate, Xe is to test;
  # err_e matrix, where Xe is used to estimate, Xo is to test.
  colnames(err_o) <- c(h_names, t_names)
  colnames(err_e) <- c(h_names, t_names)
  
  for (i in 1:Kmax) {  # i =1
    ### results fitting on Xo, i is the Khat
    cp_o <- c(0, tau_mat_o[i, 1:i], n_o)
    eB_o <- tau_mat_o[i, 1:i]
    hv <- hausdorff_dist_vector(B_o, eB_o, n_o)
    kdif <- i - length(tau2n)  # Khat - Ktrue
    obj_o <- cv_objfun(Xo, Xo, W, cp_o)
    if (obj_o < 1e-4) {obj_o <- 0;}  # IMPORTANT!!!
    obj_oe <- cv_objfun(Xo, Xe, W, cp_o)  # fit on Xo, test on Xe
    if (obj_oe < 1e-4) {obj_oe <- 0;}  # IMPORTANT!!!
    err_o[i, h_names] <- c(obj_o, obj_oe, i, kdif, hv$OE, hv$UE)
    err_o[i, t_names[1:i]] <- tau_mat_o[i, 1:i]
    ### results fitting on Xe
    cp_e <- c(0, tau_mat_e[i, 1:i], n_e)
    obj_e <- cv_objfun(Xe, Xe, W, cp_e)
    if (obj_e < 1e-4) {obj_e <- 0;}  # IMPORTANT!!!
    obj_eo <- cv_objfun(Xe, Xo, W, cp_e)
    if (obj_eo < 1e-4) {obj_eo <- 0;}  # IMPORTANT!!!
    eB_e <- tau_mat_e[i, 1:i]
    hv <- hausdorff_dist_vector(B_e, eB_e, n_e)
    err_e[i, h_names] <- c(obj_e, obj_eo,i, kdif, hv$OE, hv$UE)
    err_e[i, t_names[1:i]] <- tau_mat_e[i, 1:i]
  }
  # plot(err_e[, "K_hat"], err_e[, "objfun"])
  # plot(err_e[, "K_hat"], err_e[, "objfit"])
  # plot(err_o[, "K_hat"], err_o[, "objfun"])
  # plot(err_o[, "K_hat"], err_o[, "objfit"])
  res <- list(mat_e = err_e, mat_o = err_o)
  return(res)
}



hCV_oppelt_hausdorff <- function(
    X, W, tau2n, cp_alg = cp_hdmean_op, penaltyform = "CV", normalized = TRUE) {
  res <- copss_oppelt_hausdorff(
      X, W, tau2n, cp_alg = cp_alg, penaltyform = "CV", normalized = TRUE)
  mat_e <- res$mat_e
  mat_o <- res$mat_o
  p <- ncol(X)
  ## 1. if cv_objfun is sum, then mat_e[, "objfun"] is sum, we need
  ## mat_e[, "objfun"] = mat_e[, "objfun"] - mat_e[, "K_hat"] * p
  ## 2. if cv_objfun is mean, then mat_e[, "objfun"] is mean, we need
  ## mat_e[, "objfun"] = mat_e[, "objfun"] - mat_e[, "K_hat"] * p / n
  mat_e[, "objfun"] = mat_e[, "objfun"] - mat_e[, "K_hat"] * p
  mat_o[, "objfun"] = mat_o[, "objfun"] - mat_e[, "K_hat"] * p
  res <- list(mat_e = mat_e, mat_o = mat_o)
  return(res)
}



hCV_dp_hausdorff <- function(X, W, tau2n, Kmax = 15) {
  res <- copss_dp_hausdorff(X, W, tau2n, Kmax = 15)
  mat_e <- res$mat_e
  mat_o <- res$mat_o
  p <- ncol(X)
  ## 1. if cv_objfun is sum, then mat_e[, "objfun"] is sum, we need
  ## mat_e[, "objfun"] = mat_e[, "objfun"] - mat_e[, "K_hat"] * p
  ## 2. if cv_objfun is mean, then mat_e[, "objfun"] is mean, we need
  ## mat_e[, "objfun"] = mat_e[, "objfun"] - mat_e[, "K_hat"] * p / n
  mat_e[, "objfun"] = mat_e[, "objfun"] - mat_e[, "K_hat"] * p
  mat_o[, "objfun"] = mat_o[, "objfun"] - mat_e[, "K_hat"] * p
  res <- list(mat_e = mat_e, mat_o = mat_o)
  return(res)
}




hBIC_dp_hausdorff <- function(X, W, tau2n, alph, Kmax = 15, normalized=TRUE){
  ## pen <- calculate_penalty(X, normalized, "hBIC", alph = alph, c0 = 2.5)
  ## for X data itself, when calculating cost_matrix, normalized nor not.
  n <- nrow(X)
  B <- floor(tau2n * n)
  
  Kmax <- min(Kmax, nrow(X) / 3)   ### tauMat is Kmax*Kmax matrix
  tauMat <- cp_hdmean_dp(X, W, Kmax = Kmax)
  
  # "objfun" is the BIC obj_fun. "objfit" is the fitting error on X.
  # t_names <- c()
  # for (k in 1:Kmax) {t_names <- c(t_names, paste0("tau_", k))}
  # mapply(paste, "tau", sep='_', 1:Kmax, SIMPLIFY=TRUE, USE.NAMES = FALSE)
  # mapply(paste0, "tau_", 1:Kmax, SIMPLIFY=TRUE, USE.NAMES = FALSE)
  
  t_names = mapply(paste, "tau", sep='_', 1:Kmax, USE.NAMES = FALSE)
  h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  errMat <- matrix(NA, Kmax, length(h_names) + Kmax)  # fitting error on X.
  colnames(errMat) <- c(h_names, t_names)
  
  #### alph = 0.1; c0 = 2.5
  # calculate_penalty(X, normalized, penaltyform = "hBIC", alph, c0)
  ## pen <- calculate_penalty(X, normalized, "hBIC", alph = alph, c0 = 2.5)
  # p <- ncol(X)
  # n <- nrow(X)
  # p + c0 * trace_R_square(X)^(0.5) * log(n)^(1+alph)
  pen <- calculate_penalty(X, normalized, "hBIC", alph = alph, c0 = 2.5)
  for (k in 1:Kmax) {  # k = 1
    cp_tau <- c(0, tauMat[k, 1:k], n)
    eB <- tauMat[k, 1:k]
    hv <- hausdorff_dist_vector(B, eB, n)
    kdif <- k - length(B)
    obj_fit <- cv_objfun(X, X, W, cp_tau)
    obj_bic <- obj_fit + k * pen;  # L = k
    if (obj_fit < 1e-4) {obj_fit <- 0;}  # IMPORTANT!!!
    if (obj_bic < 1e-4) {obj_bic <- 0;}  # IMPORTANT!!!
    errMat[k, h_names] <- c(obj_fit, obj_bic, k, kdif, hv$OE, hv$UE)
    errMat[k, t_names[1:k]] <- tauMat[k, 1:k]
  }
  
  # errMat <- errMat[order(errMat[, "K_hat"]), ]
  # diff(errMat[,"objfit"]) > 0
  # plot(errMat[, "K_hat"], errMat[, "objfit"])
  return(errMat)
}


#' Calculate the fit error matrix using tau matrix obtained by DP algorithm 
#' 
#' Calculate the fit error matrix using tau matrix obtained by DP algorithm 
#' 
#' Calculate the fit error matrix using tau matrix obtained by DP algorithm
#' and it includes the hausdorff distance. 
#' The function is often used in BIC methods.
#' 
#' @param X N*p matrix
#' @param W vector with p-length as the global normalizer.
#' @param tau2n the true tau/n ratio.
#' @param tauMat candidate tau's matrix estimated by DP based on data Xtrain
#' @return An error matrix, with c("objfit", "K_hat", "K_diff", "OE", "UE") 
#'   and "tau_1", "tau_Kmax" as columns.
dpTauMat2FitErrMat <- function(X, W, tau2n, tauMat){
  Kmax <- nrow(tauMat)
  # mapply(paste, "tau", sep='_', 1:Kmax, SIMPLIFY=TRUE, USE.NAMES = FALSE)
  # mapply(paste0, "tau_", 1:Kmax, SIMPLIFY=TRUE, USE.NAMES = FALSE)
  t_names <- mapply(paste, "tau", sep='_', 1:Kmax, USE.NAMES = FALSE)
  h_names <- c("objfit", "K_hat", "K_diff", "OE", "UE")
  errMat <- matrix(NA, Kmax, length(h_names) + Kmax)  # fitting error on X.
  colnames(errMat) <- c(h_names, t_names)
  
  n <- nrow(X)
  B <- floor(tau2n * n)
  
  for (k in 1:Kmax) {  # k =1
    cp_tau <- c(0, tauMat[k, 1:k], n)
    eB <- tauMat[k, 1:k]
    haus <- hausdorff_dist_vector(B, eB, n)
    kdif <- k - length(B)
    obj_fit <- cv_objfun(X, X, W, cp_tau)
    if (obj_fit < 1e-4) {obj_fit <- 0;}  # IMPORTANT!!!
    errMat[k, h_names] <- c(obj_fit, k, kdif, haus$OE, haus$UE)
    errMat[k, t_names[1:k]] <- tauMat[k, 1:k]
  }
  # errMat <- errMat[order(errMat[, "K_hat"]), ]
  # diff(errMat[,"objfit"]) > 0
  # plot(errMat[, "K_hat"], errMat[, "objfit"])
  # return(list(mat_bic=errMat))
  return(errMat)
} 



hBIC_dp_hausdorff_bundle <- function(
    X, W, tau2n, c0, alph_list, Kmax = 15, normalized=TRUE) {
  ## pen <- calculate_penalty(X, normalized, "hBIC", alph = alph, c0 = 2.5)
  ## for X data itself, when calculating cost_matrix, normalized nor not.
  ## The best is c0 = 2.5, whatever alph is, the penalty is not too small
  ## An error matrix, with c("hBIC_0", "hBIC_0.1",..., "objfit", "K_hat", 
  ## "K_diff", "OE", "UE") and c("tau_1",  ..., "tau_Kmax") as columns.
  Kmax <- min(Kmax, nrow(X) / 3)   ### tauMat is Kmax*Kmax matrix
  tauMat <- cp_hdmean_dp(X, W, Kmax = Kmax)
  
  errMat <- dpTauMat2FitErrMat(X, W, tau2n, tauMat)
  
  # t_names = paste("tau", 1:6, sep='_')
  bic_names <- mapply(paste, "hBIC", sep='_', alph_list, USE.NAMES = FALSE)
  bicMat <- matrix(NA, Kmax, length(bic_names))  # fitting error on X.
  colnames(bicMat) <- bic_names
  #### alph = 0.1; c0 = 2.5
  # calculate_penalty(X, normalized, penaltyform = "hBIC", alph, c0)
  # pen <- calculate_penalty(X, normalized, "hBIC", alph = alph, c0 = 2.5)
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (normalized) {tr2 <- trace_R_square(X)} else {tr2 <- p}
  for (alph in alph_list) {
    pen <- p + c0 * tr2^(0.5) * log(n)^(1+alph)
    temp <- errMat[, "objfit"] + errMat[, "K_hat"] * pen
    bicMat[, paste("hBIC", alph, sep='_')] <- temp
  }
  bicMat <- cbind(bicMat, errMat)
  return(bicMat)
}



dpTauMat2cvErrMat <- function(Xtrain, Xtest, W, B, tauMat) {
  # tauMat = tauMat_e; B=Be; Xtrain=Xo; Xtest=Xe
  ## W: vector with p-length as the global normalizer.
  ## tauMat, B, are based on the Xtrain data
  Kmax <- nrow(tauMat)
  # t_names <- c()
  # for (k in 1:Kmax) {t_names <- c(t_names, paste0("tau_", k))}
  t_names <- paste("tau", c(1:Kmax), sep = "_")
  h_names <- c("objfit", "objfun", "K_hat", "K_diff", "OE", "UE")
  errMat <- matrix(NA, Kmax, length(h_names) + Kmax)
  colnames(errMat) <- c(h_names, t_names)
  n <- nrow(Xtrain)
  for (k in 1:Kmax) {  # k = 1
    ### results fitting on Xtrain, k is the Khat
    cp_hat <- c(0, tauMat[k, 1:k], n)
    eB <- tauMat[k, 1:k]
    haus <- hausdorff_dist_vector(B, eB, n)
    kdif <- k - length(B)  # Khat - Ktrue
    obj_fit <- cv_objfun(Xtrain, Xtrain, W, cp_hat)  # fit on Xtrain
    obj_cv <- cv_objfun(Xtrain, Xtest, W, cp_hat)  # fit on Xtrain, test on Xtest
    if (obj_fit < 1e-4) {obj_fit <- 0;}  # IMPORTANT!!!
    if (obj_cv < 1e-4) {obj_cv <- 0;}  # IMPORTANT!!!
    errMat[k, h_names] <- c(obj_fit, obj_cv, k, kdif, haus$OE, haus$UE)
    errMat[k, t_names[1:k]] <- tauMat[k, 1:k]
  }
  return(errMat)
}



copss_dp_ErrMat_hausdorff <- function(X, W, tau2n, Kmax = 15) {
  ### copss_dp_ErrMat_hausdorff is the same as copss_dp_hausdorff.
  ## normalized: 
  ## for X data itself, when calculating cost_matrix(X), normalized nor not.
  ## W: vector with p-length as the global normalizer.
  ## mat_e: CV error matrix on 'even' data, which uses 'odd' data to estimate.
  ## mat_o: CV error matrix on 'odd' data, which uses 'even' data to estimate.
  ## err_e matrix, CV error based on Xe (test data), where Xo is train data
  ## err_o matrix, CV error based on Xo (test data), where Xe is train data
  ## Xo, Xe, Xodd, Xeven
  Xo <- X[seq(1, nrow(X), 2), ]  # max(seq(1, nrow(X), 2))
  Xe <- X[seq(2, nrow(X), 2), ]  # max(seq(2, nrow(X), 2))
  if (!is.matrix(Xo)) {Xo = matrix(Xo, length(Xo), 1)}
  if (!is.matrix(Xe)) {Xe = matrix(Xe, length(Xe), 1)}
  
  Bo <- floor(tau2n * nrow(Xo))
  Be <- floor(tau2n * nrow(Xe))
  
  ## W: vector with p-length as the global normalizer.
  Kmax <- min(Kmax, nrow(X) / 3)   ### tau_mat is Kmax*Kmax matrix
  tauMat_o <- cp_hdmean_dp(Xo, W, Kmax = Kmax)
  tauMat_e <- cp_hdmean_dp(Xe, W, Kmax = Kmax)
  err_o <- dpTauMat2cvErrMat(Xo, Xe, W, Bo, tauMat_o)
  err_e <- dpTauMat2cvErrMat(Xe, Xo, W, Be, tauMat_e)
  ## err_o matrix, where Xo is used to train, Xe for CV test error.
  ## err_e matrix, where Xe is used to train, Xo for CV test error.
  
  # plot(err_e[, "K_hat"], err_e[, "objfun"])
  # plot(err_e[, "K_hat"], err_e[, "objfit"])
  # plot(err_o[, "K_hat"], err_o[, "objfun"])
  # plot(err_o[, "K_hat"], err_o[, "objfit"])
  res <- list(mat_e = err_e, mat_o = err_o)
  return(res)
}





