############################## ZouSettings ##############################
##### CP(A) - Both the number and the locations are fixed. #####
# Kn <- 11
# tau2n <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81);
# n <- 2048

##### CP(B) - Both the number and the locations can vary #####
##### with the sample size n.


# PS: TO DEAL with the numerical small values!
# if (obj_o < 1e-4) {obj_o <- 0;}  # IMPORTANT!!!


cpb_Kn <- function(n) {floor(log(n)^1.01)}


#' Change-point model B -- generate the true change-point locations
#' 
#' Change-point model B -- generate the true change-point locations
#' 
#' @param n sample size
#' @return the locations of the change-points, without including 0 and n.
#' @examples 
#'   tau <- cpb_true_tau(2048)
#'   n <- 400
#'   tau_true <- cpb_true_tau(n)
#'   tau_vec <- c(0, sort(tau_true), n)
#'   min(diff(tau_vec))
cpb_true_tau <- function(n) {
  Kn <- cpb_Kn(n)
  a <- floor(n^(1/4))
  tau_true <- c()
  for (j in 1:Kn) {  # j = 1
    tau_true[j] = j * floor(n/(Kn + 1)) + sample(-a:a, 1)
  }
  tau_true <- sort(tau_true)
  return(tau_true)
}


#' Generate piecewise constant vector 
#' 
#' Generate piecewise constant vector in 1D with n-length
#' 
#' Generate piecewise constant data with tau_true as the true change-points.
#' 
#' Generate piecewise constant vector with length n, 
#' which starts with 0, such as [0, a, a, a, 0, 0,...]. 
#' The number of the repeated a's or repeated 0's are determine by 
#' \eqn{\tau_{j+1} - \tau_{j}}.
#' 
#' @param n the length of the generated vector, same as the sample size
#' @param tau_true true change-points, \eqn{\tau}'s not the tau2n's, 
#'   without the 0 nor n included.
#'   
#' @param umax the max value of each element in the \eqn{\theta}-vector.
#'   \eqn{\theta} is further repeated to generate the result vector. 
#'   In other words, umax is the max value of the generated vecor.
#' @return a n-length vetor
#' 
#' @examples 
#' n <- 200
#' tau_true <- c(50, 100)
#' x <- gen_piecewise_const_vector(n, tau_true, umax=3)
gen_piecewise_const_vector <- function(n, tau_true, umax=1) {
  ### generate piecewise constant data with tau_vec as true change-points.
  Kn <- length(tau_true)
  ### (mean_true <- rlnorm(Kn + 1, meanlog = 0, sdlog = log2(10)/2))
  ### generating the parameter of interest, theta_true
  if ((Kn + 1) %% 2 == 0) {
    theta_true <- rep(c(0, umax), floor((Kn + 1) / 2))  # [0, 1, 0, 1, 0, 1]
  } else {
    theta_true <- c(rep(c(0, umax), floor((Kn + 1)/2)), 0)  # [0, 1, 0, 1, 0]
  }
  ### length(theta_true) == Kn + 1
  tau_vec <- c(0, tau_true, n)  # 0 and n are included.
  
  xsignal <- c()
  for (i in 1:(Kn+1)) {
    ni <- tau_vec[i+1] - tau_vec[i]
    xsignal <- c(xsignal, rep(theta_true[i], ni))
  }
  return(xsignal)
}



#' Generate theta vector
#' 
#' Generate theta vecor with (Kn+1)-length
#' 
#' We need to set the velues of n, p, r, umax, which finally can generate proper 
#' \eqn{\delta_{n,p}}, thus the signal matrix, \eqn{\mu} matrix. 
#' In order to generate proper \eqn{\delta_{n,p}}, we need \eqn{\theta}. 
#' The \eqn{\theta}'s vector takes value of 0, umax alternatively.
#' Then \eqn{\lambda_{n}} defines the minimum distance between true change-points,
#' \deqn{\lambda_{n} = min(\tau_{j+1} - \tau_{j}) . }
#' The \eqn{\delta_{n,p}} is defined as following, 
#' \deqn{\delta_{n,p} = min \| D^{-1/2}(\mu_{j+1} - \mu_{j}) \|^2 .} 
#' At the same time, approximately we have \eqn{\delta_{n,p}} of magnitude 
#'   \deqn{\delta_{n,p} = p * r * umax^2 , }
#' with each element being umax in \eqn{\bm{\mu}} .
#' 
#' \deqn{D = diagonal(\sigma_{j}^2)}
#' \deqn{signal = \lambda_{n} \delta_{n,p}}
#' 
#' Notes: \eqn{\lambda_{n}} is only about the positions, wheras,
#' \eqn{\delta_{n,p}} is about both the positions and the \eqn{\mu}'s.
#' 
#' Notes: e.g. tau2n = c(0.2, 0.4, 0.6), then Kn=3 not 5.
#' \itemize{
#' \item{If (Kn+1) is even, the generated vector will be
#'    [0, umax, 0, umax, ..., 0, umax] which ends with umax;}
#' \item{If (Kn+1) is odd, the generated vector will be
#'    [0, umax, 0, umax, ..., 0, umax, 0] which ends with 0.}
#' }
#' 
#' @param Kn  the number of change-points
#' @param umax the constant (or the max value) of each element of
#'   the generated vector
#' @return a (Kn+1)-length vector, which can be viewed as the base vector,
#'   which would further generate the \eqn{\mu}'s matrix.
#' 
#' @examples 
#'   tau2n <- c(0.2, 0.4, 0.6)  # then Kn=3 not 5.
#'   Kn <- length(tau2n)
#'   theta_true <- genconst_Kn_theta(Kn, umax=4)
genconst_Kn_theta <- function(Kn, umax=1){
  # Kn is the number of change-points
  # umax is the constant or max constant of each element
  # returns a (Kn+1)-length vector, which can be viewed as the mu's matrix.
  # e.g. tau2n = c(0.2, 0.4, 0.6), then Kn=3 not 5.
  # [0, a, 0, a, ..., 0, a] if (Kn+1) is even
  # [0, a, 0, a, ..., 0, a, 0] if (Kn+1) is odd
  if ((Kn + 1) %% 2 == 0) {
    theta_true <- rep(c(0, umax), floor((Kn + 1)/2))  # [0, 1, ... 0, 1]
  } else {
    theta_true <- c(rep(c(0, umax), floor((Kn + 1)/2)), 0)  # [0, 1, 0, 1, 0]
  }
  return(theta_true)
}


#' Generate (Kn+1)*p true signal matrix \eqn{\mu}
#' 
#' Generate (Kn+1)*p true signal matrix \eqn{\mu} in \deqn{X = \mu + \epsilon}
#' 
#' This function generate (Kn+1)*p true signal matrix \eqn{\mu} by extending 
#' theta_true[i] to a p-len vector with \eqn{nr = p * r} components being 
#' theta_true[i] (make sure nr is at least one).
#' 
#' @inheritParams genconst_Kn_theta
#' @param p dimension
#' @param r the ratio where floor(p*r) dims in total have a change.
#' @param diff_signal logical.
#'   Allow different signals along different dimensions or not.
#'   There are three options: TRUE, FALSE, "permutation", which is corresponding
#' @return a (Kn+1)*p matrix, with each row as a mean vector in each
#'   segmentation, for example, row-vector generated by 
#'   \code{\link{genconst_Kn_theta}}
genconst_mat_signal <- function(Kn, umax, p, r, diff_signal){
  ## This function generate (Kn+1)*p true mu matrix by extending 
  ## theta_true[i] to a p-len vector with nr components being theta_true[i].
  ## to make sure nr is at least one.
  ## diff_signal = "permutation"
  nr <- floor(p * r)
  if (nr == 0) {nr <- 1}
  
  ### Type2: generate piecewise constant mean data, [0, mu, 0, mu, ...]
  ### mu = umax * (1, 1/2, ..., 1/nr, 0, ..., 0)^{1/2}
  ### only the first nr components/(or columns) are non-zeroes 
  if (diff_signal == TRUE) {
      cj <- 1 /sqrt(1:nr)
      hmat <- matrix(0, Kn+1, p)
      for (i in 1:(Kn+1)) {hmat[i, 1:nr] <- cj * umax }
      m <- ceiling((Kn + 1) / 2)
      hmat[2*c(1:m) - 1, ] <- 0
      return(hmat)
  }
  
  ## sample(1:10), a random permutation, permutation can change the seed
  if (r == 1) {J <- c(1:p)} else {J <- sample(p, nr)}
  ## allow different signals along different dimensions. permutation once.
  if (diff_signal == "permutation") { cj <- sample(nr) / nr } else {cj <- 1}
  
  ### Type1: generate piecewise constant mean data, [0, mu, 0, mu, ...]
  ### mu = umax * I_p
  ## delta_np = nr * umax^2  # with each element in \mu being umax
  theta_true <- genconst_Kn_theta(Kn, umax)
  hmat <- matrix(0, Kn+1, p)
  for (i in 1:(Kn+1)) {hmat[i, J] <- cj * theta_true[i]}
  # muJ = hmat[, J]
  
  return(hmat)
}


#' Generate IND residual
#' 
#' Generate independent residual data
#' 
#' @param n integer, sample size
#' @param p integer, dimension
#' @param residualform string, the residual form, which can take value of 
#'   "normal", "t", "chisq", "uniform", 
#'   or "quasinormal" but only when p = 1.
#'   All the variance of these residual data are normalized to 1.
#'   So we can only focus on different shapes of different distrubution.
#' @return A random matrix of n*p.
gen_IND_residual <- function(n, p, residualform="normal") {
  ### generate the independent data based on different residual form.
  ### We need to normalize the variance s.t. variance = 1. 
  ### So we can only focus on different shapes.
  if (residualform == "normal") {
    eps <- matrix(rnorm(n * p, 0, 1), n, p)
  } else if (residualform == "t"){
    ### t_{v}, the df = v, then the mean is 0, and the variance is v / (v-2)
    tmp <- rt(n * p, df = 5) / sqrt(5/3)
    eps <- matrix(tmp, n, p)
  } else if (residualform == "chisq"){
    ### chisq_{v}, the df = v, then the mean is v, and the variance is 2*v
    tmp <- (rchisq(n * p, df = 5, ncp = 0) - 3) / sqrt(6)
    eps <- matrix(tmp, n, p)
  } else if (residualform == "uniform") {
    ### Uniform(a, b), mean (a+b)/2, variance (b-a)^2/12
    a <- (-sqrt(3))
    b <- sqrt(3)
    s <- sqrt((b - a)^2 / 12)  # Here s = 1
    tmp <- runif(n * p, min = a, max = b) / s
    eps <- matrix(tmp, n, p)
  } else if (residualform == "quasinormal" & p == 1) {
    tmp <- sin(2*pi * c(1:n)) / sum(sin(2*pi * c(1:n))) * rnorm(n, 0, 1)
    eps <- matrix(tmp, n, 1)
  } else {
    eps <- matrix(rnorm(n * p, 0, 1), n, p)
  }
  return(eps)
}


#' Generate 1 dim data
#' 
#' CP(A) - MODEL I: univariate mean change-point model.
#' 
#' @section Details:
#' CP(A) - Both the number and locations of changepoints are fixed.
#' 
#' CP(B) - Both the number and locations of changepoints can vary with the 
#' sample size n.
#' 
#' MODEL I: univariate mean change-point model.
#' 
#' We adopt the blocks setting which is widely used in the literature.
#' 
#' The signals \eqn{theta_{j}}'s (the mean) are chosen as a piecewise constant 
#' function with Kn and the scale parameter \code{sig_const = 0.5}.
#' \deqn{x_{i} = \theta_{j} + sig_const * \epsilon_{i}}.
#' \eqn{\theta_{j}} can take values alternatively between 0 and 1.
#' 
#' Specifically, \eqn{K_n = 11} and 
#' \deqn{ \tau^{*}_{K_{n}} / n = 
#'   (0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81 ). }
#' \eqn{K_n = 11} is the number of the true change-points \eqn{\tau^{*}}.
#' 
#' @inheritParams gen_IND_residual
#' @param sig_const scalar, default 0.5, the constant before the residual
#' @param verbose bool, default FALSE, whether to return a list of 
#'   data with (X, mu, eps) or X.
#' @return A random vector of length n.
#' 
#' @examples 
#' x <- generate_1D_data(n=2048, residualform="normal", sig_const=0.5)
generate_1D_data <- function(
    n=2048, residualform="normal", sig_const=0.5, verbose=FALSE) {
  ## n=2048; residualform="normal"; sig_const=0.5;
  # Kn: the number of change-points.
  # tau_true: stores the change-points locations, including the 0 and n.
  # theta_true: stores the true parameter of interest, theta_j.
  tau2n <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
  Kn <- length(tau2n)  # Kn <- 11
  
  ### generating the parameter of interest, theta_true
  # (mean_true <- rlnorm(Kn + 1, meanlog = 0, sdlog = log2(10)/2))
  if ((Kn + 1) %% 2 == 0) {
    theta_true <- rep(c(0, 1), floor((Kn + 1) / 2))  # [0, 1, 0, 1, 0, 1]
  } else {
    theta_true <- c(rep(c(0, 1), floor((Kn + 1)/2)), 0)  # [0, 1, 0, 1, 0]
  }
  # length(theta_true) == Kn + 1
  
  tau_true <- floor(tau2n * n)
  tau_vec <- c(0, tau_true, n)  # 0 and n are included.
  xsignal <- c()
  for (i in 1:(Kn+1)) {
    ni <- tau_vec[i+1] - tau_vec[i]
    xsignal <- c(xsignal, rep(theta_true[i], ni))
  }
  
  # generate the residual
  eps <- gen_IND_residual(n, p = 1, residualform = residualform)
  X <- xsignal + sig_const * eps  # X is n*1
  # plot(X, type='p')
  # for (i in 2:Kn) abline(v=tau_true[i], col='red')
  if (verbose){
    return(list(X=X, Xu=xsignal, eps=eps))
  } else {
    return(X)
  }
  
}


############################## WangSettings ##############################
### The mean-change model: X_i = mu_k + Sigma_matrix^{1/2} * z_i
### Generating error term eps_i = Sigma_matrix^{1/2} * z_i from an ARMA model,
### eps_{ij} - beta_{1} eps_{i, j-1} = e_{ij} + beta_{2} e_{i,j-1}
### (correlations between the j-1 and the j-th dimensions)
### where e_{ij} is the white noise with zero mean.
##### Four correlation structures #####
### (1) beta_{1}, beta_{2} = 0, 0   the independent(IND) case;
### (2) beta_{1}, beta_{2} = 3/4, 0   the autoregressive(AR) case;
### (3) beta_{1}, beta_{2} = 0, 3/4   the moving-average(MA) case;
### (4) beta_{1}, beta_{2} = 1/2, 1/2   the ARMA case;
##### Three scenarios for the innovation e_{ij}, residualform #####
### (a) N(0, 1), normal 
### (b) t(5), Student's t distribution with five degree of freedom.
### (c) chi-squared(3), chi-squared distribution with three degree of freedom.
### ALL simulation results are based on 1000 replications.
##### Two different data generation processes: #####
##### Model(1) - same signals along all the changed dimensions. #####
### We randomly assign p*r components to be changed from mu_k to mu_{k+1}
#####  Model(2) - different signals along different dimensions #####
### allow the signal strengths along different dimensions to be different.
### 



#' Generate ARMA residual
#' 
#' Generate ARMA residual data
#' 
#' The mean-change model: 
#' \deqn{X_i = \mu_{k} + \Sigma^{1/2} * z_i}
#' 
#' Generating error term 
#' \eqn{\epsilon_{i} = \Sigma^{1/2} * z_i} from an ARMA model,
#' which specifies the correlations between the j-1 and the j-th dimensions, 
#' i.e. correlations between \eqn{\epsilon_{i, j-1}} and \eqn{\epsilon_{i, j}}.
#' \deqn{\epsilon_{ij} - \beta_{1} \epsilon_{i, j-1} = e_{ij} + \beta_{2} e_{i,j-1}}
#' where \eqn{e_{ij}} and \eqn{e{i, j-1}} are the white noise with zero mean.
#' 
#' Four correlation structures:
#' \itemize{
#'   \item{(1) \eqn{\beta_1, \beta_{2} = 0, 0}   the independent(IND) case;}
#'   \item{(2) \eqn{\beta_1, \beta_{2} = 3/4, 0}   the autoregressive(AR) case;}
#'   \item{(3) \eqn{\beta_1, \beta_{2} = 0, 3/4}   the moving-average(MA) case;}
#'   \item{(4) \eqn{\beta_1, \beta_{2} = 1/2, 1/2}   the ARMA case;}
#'  }
#'  
#' @inheritParams gen_IND_residual
#' @param b1 \eqn{\beta_1}
#' @param b2 \eqn{\beta_2}
#' @return A random matrix n*p
gen_ARMA_residual <- function(n, p, b1=0, b2=0, residualform="normal") {
  ### X_{ij} - beta_{1} * X_{i, j-1} = Z_{ij} + beta_{2} Z_{i,j-1}
  ### X_{t} - phi * X_{t-1} = Z_{t} + theta * Z_{t-1}
  ### Z_{1:t} are white noise with mean zero and variance of sigma^2
  ### usually to make it easier we take sigma = 1.
  ### Here t is along the dimension p,
  ### generate ARMA data based on different residual form.
  Z <- gen_IND_residual(n, p, residualform)
  X <- matrix(NA, n, p)
  X[, 1] <- Z[, 1]
  for (i in 1:n){
    for (j in 2:p){
      X[i, j] <- b1 * X[i, j-1] + Z[i, j] + b2 * Z[i, j-1]
    }
  }
  return(X)
}


#' Generate independent data in mcp model
#' 
#' Generate data, independent case, multiple change-point model (mean change)
#' 
#' IND case: {X_i = u_k + sig_const * eps_i}. 
#'   X_i is an obsevation vector with p-dim. 
#'   Independent case has 0 correlation between components of eps_{i}, 
#'   namely 0 correlation between eps_{ij} and eps_{i,j-1}.
#' 
#' @inheritParams gen_IND_residual
#' @inheritParams genconst_mat_signal
#' @param n Sample size
#' @param p Dimension
#' @param r Sparsity ratio. There are p*r dims that have change-points.
#'   Set r=1 when it is in low dimension, i.e. p is small.
#' @param tau2n The locations of change-points relative to the sample size,
#'   without 0 and 1 included, just the pure change-points.
#' @param diff_signal logical.
#'   Allow different signals along different dimensions or not.
#' @param sig_const scalar, default 1, the constant before the residual.
#' @param amplify_signal scalar, default 1, the constant before the signal.
#' @param verbose bool, default FALSE, whether to return a list of 
#'   data with (X, mu, eps) or X.
#' @return X data matrix n*p
#' @examples 
#' n <- 200; p <- 500; r <- 0.05; umax <- 0.9
#' tau2n <- c(0.2, 0.40, 0.60, 0.80)
#' residualform="normal"; diff_signal=FALSE; sig_const=1;  amplify_signal=1
#' set.seed(12345)
#' X <- gen_data_IND_mcp(n, p, r, umax, tau2n, residualform="normal", 
#'   diff_signal=FALSE, sig_const=1, amplify_signal=1)
gen_data_IND_mcp <- function(
    n, p, r, umax, tau2n, residualform="normal", 
    diff_signal=FALSE, sig_const=1, amplify_signal=1, verbose=FALSE) {
  ### generate the mu_mat, which is (Kn+1)*p matrix.
  Kn <- length(tau2n)
  hmat <- genconst_mat_signal(Kn, umax, p, r, diff_signal)
  
  ### generate the signal -- Xu
  Theta <- matrix(0, n, 0)
  tau_true <- floor(tau2n * n)
  taut <- c(0, tau_true, n)
  for (j in 1:(Kn + 1)){  # j = K0+1; j = 1
    temp <- c(rep(0, taut[j]), rep(1, taut[j+1] - taut[j]), rep(0, n - taut[j+1]))
    Theta <- cbind(Theta, temp)
  }
  Xu <- Theta %*% hmat
  
  ### generate the IND residual
  eps <- gen_IND_residual(n, p, residualform)
  Xmat <- amplify_signal * Xu + sig_const * eps
  if (verbose){
    return(list(Xmat=Xmat, Xu=Xu, eps=eps))
    # A = list(Xmat=Xmat, Xu=Xu, eps=eps)
  } else {
    return(Xmat)
  }
}



#' Generate ARMA data in mcp model
#' 
#' Generate data, ARMA case, multiple change-point model (mean change)
#' 
#' @inheritParams gen_ARMA_residual
#' @inheritParams gen_data_IND_mcp
#' @inheritParams genconst_mat_signal
#' @return X data matrix n*p
#' @examples
#' n <- 200; p <- 100; r <- 0.1; umax <- 0.8
#' n <- 200; p <- 10; r <- 1; umax <- 0.8
#' tau2n <- c(0.2, 0.40, 0.60, 0.80)
#' residualform="normal"; diff_signal=FALSE; sig_const=1;  amplify_signal=1
#' set.seed(12345)
#' b1 <- 3/4; b2 <- 0
#' X1 <- gen_data_ARMA_mcp(n, p, r, umax, tau2n, b1, b2, 
#'   residualform="normal", diff_signal=FALSE, amplify_signal=1)
gen_data_ARMA_mcp <- function(
    n, p, r, umax, tau2n, b1, b2, residualform="normal", 
    diff_signal=FALSE, sig_const=1, amplify_signal=1, verbose=FALSE) {
  ### generate the mu_mat, which is (Kn+1)*p matrix.
  Kn <- length(tau2n)
  hmat <- genconst_mat_signal(Kn, umax, p, r, diff_signal)
  
  tau_true <- floor(tau2n * n)
  taut <- c(0, tau_true, n)
  Theta <- matrix(0, n, 0)
  for (j in 1:(Kn + 1)){  # j = K0+1; j = 1
    temp <- c(rep(0, taut[j]), rep(1, taut[j+1] - taut[j]), rep(0, n - taut[j+1]))
    Theta <- cbind(Theta, temp)
  }
  Xu <- Theta %*% hmat
  
  ## b1=0; b2=0; residualform="normal"; diff_signal=FALSE; amplify_signal=1
  ### generate the ARMA residual
  eps <- gen_ARMA_residual(n, p, b1, b2, residualform)
  Xmat <- amplify_signal * Xu + sig_const * eps
  if (verbose){
    return(list(Xmat=Xmat, Xu=Xu, eps=eps))
    # A = list(Xmat=Xmat, Xu=Xu, eps=eps)
  } else {
    return(Xmat)
  }
}


#' Generate high-dimension CV ARMA data
#' 
#' Generate ARMA data for high-dimension CV 
#' 
#' Assume diff_signal = FALSE forever; and amplify_signal=1 forever
#' 
#' @inheritParams gen_IND_residual
#' @inheritParams gen_ARMA_residual
#' @inheritParams genconst_mat_signal
#' @inheritParams gen_data_ARMA_mcp
#' 
#' @param umax maximum element in the \eqn{\mu = (0, a, 0, a, \cdots)}
#' @param cor_structure "IND", "AR", "MA", "ARMA"
#' @return X data matrix n*p
#' @examples 
#'   tau2n <- c(0.2, 0.4, 0.6, 0.8)
#'   n <- 40; p <- 100; r <- 1; umax <- 1;
#'   dt <- gen_hcv_ARMA_data(
#'       n, p, r, umax, tau2n, cor_structure="IND", 
#'       residualform="normal", 
#'       diff_signal=TRUE, verbose=TRUE)
#'   Xu <- dt$Xu
gen_hcv_ARMA_data <- function(
    n, p, r, umax, tau2n, cor_structure="IND", residualform="normal", 
    diff_signal=FALSE, verbose=FALSE){
  # umax means the max value of each element in the mu-vector
  # Assume D = I_p
  # delta_np = min \|D^{-1/2}(\mu_{j+1} - \mu_{j})\|^2
  # for the current data generating system 
  # delta_np = umax^2 * p * r
  # n = 400; p = 100; r = 0.1; umax=1; diff_signal = TRUE
  Kn <- length(tau2n)
  hmat <- genconst_mat_signal(Kn, umax, p, r, diff_signal)
  
  # ## generate the signal -- Xu matrix
  # tau_true <- floor(tau2n * n)
  # tau_vec <- c(0, tau_true, n)  # 0 and n are included.
  # Xu <- vector()
  # for (i in 1:(Kn + 1)) {
  #   u_i <- hmat[i, ]
  #   ni <- tau_vec[i+1] - tau_vec[i]
  #   mu_i <- matrix(rep(u_i, ni), ni, p, byrow = TRUE)
  #   Xu <- rbind(Xu, mu_i)
  # }
  
  ## generate the signal -- Xu matrix
  tau_true <- floor(tau2n * n)
  taut <- c(0, tau_true, n)
  Theta <- matrix(0, n, 0)
  Kn <- length(tau_true)
  for (j in 1:(Kn + 1)){
    temp <- c(rep(0, taut[j]), rep(1, taut[j+1] - taut[j]), rep(0, n - taut[j+1]))
    Theta <- cbind(Theta, temp)
  }
  Xu <- Theta %*% hmat
  
  ### matching the proper correlation structure ### 
  if (cor_structure == "IND") {
    b1 <- 0; b2 <- 0
  } else if (cor_structure == "AR") {
    b1 <- 3/4; b2 <- 0
  } else if (cor_structure == "MA") {
    b1 <- 0; b2 <- 3/4
  } else if (cor_structure == "ARMA") {
    b1 <- 1/2; b2 <- 1/2
  }
  
  eps <- gen_ARMA_residual(n, p, b1, b2, residualform)
  Xmat <- Xu + eps
  if (verbose){
    return(list(Xmat=Xmat, Xu=Xu, eps=eps))
    # A = list(Xmat=Xmat, Xu=Xu, eps=eps)
  } else {
    return(Xmat)
  }
}


#' Generate CP(B) data
#' 
#' Generate change-point Model B data
#' 
#' @inheritParams gen_hcv_ARMA_data
#' @inheritParams genconst_mat_signal
#' 
#' @return X data matrix n*p
#' @examples 
#'   tau2n <- c(0.2, 0.4, 0.6, 0.8)
#'   Kn <- length(tau2n)
#'   n <- 40; p <- 100; r <- 1; umax <- 1
#'   dt <- gen_CPB_data(
#'       n, p, r, umax, cor_structure="IND", residualform="normal", 
#'       diff_signal=TRUE, verbose=TRUE)
#'   Xu <- dt$Xu
gen_CPB_data <- function(
    n, p, r, umax, cor_structure="IND", residualform="normal", 
    diff_signal=FALSE, verbose=FALSE){
  # umax means the max value of each element in the mu-vector
  # Assume D = I_p
  # delta_np = min \|D^{-1/2}(\mu_{j+1} - \mu_{j})\|^2
  # for the current data generating system 
  # delta_np = umax^2 * p * r
  ##### Under CP(B), the locations of change-points vary from replication 
  ##### to replication. 
  ##### generate the tau2n of CPB #####
  amplify_signal <- 1
  tau_true <- cpb_true_tau(n)
  tau_vec <- c(0, tau_true, n)  # 0 and n are included.
  tau2n <- tau_true / n
  # tau_true - floor(tau2n * n)  # there maybe some difference
  Kn <- length(tau2n)
  hmat <- genconst_mat_signal(Kn, umax, p, r, diff_signal)
  
  # ## generate the signal -- Xu matrix
  # Xu <- vector()
  # for (i in 1:(Kn + 1)) {
  #   u_i <- hmat[i, ]
  #   ni <- tau_vec[i+1] - tau_vec[i]
  #   mu_i <- matrix(rep(u_i, ni), ni, p, byrow = TRUE)
  #   Xu <- rbind(Xu, mu_i)
  # }
  
  ## generate the signal -- Xu matrix
  tau_true <- floor(tau2n * n)
  taut <- c(0, tau_true, n)
  Theta <- matrix(0, n, 0)
  Kn <- length(tau_true)
  for (j in 1:(Kn + 1)){
    temp <- c(rep(0, taut[j]), rep(1, taut[j+1] - taut[j]), rep(0, n - taut[j+1]))
    Theta <- cbind(Theta, temp)
  }
  Xu <- Theta %*% hmat
  
  ### matching the proper correlation structure ### 
  if (cor_structure == "IND") {
    b1 <- 0; b2 <- 0
  } else if (cor_structure == "AR") {
    b1 <- 3/4; b2 <- 0
  } else if (cor_structure == "MA") {
    b1 <- 0; b2 <- 3/4
  } else if (cor_structure == "ARMA") {
    b1 <- 1/2; b2 <- 1/2
  }
  
  eps <- gen_ARMA_residual(n, p, b1, b2, residualform)
  Xmat <- Xu + eps
  if (verbose){
    return(list(Xmat=Xmat, Xu=Xu, eps=eps, tau2n=tau2n))
    # A = list(Xmat=Xmat, Xu=Xu, eps=eps)
  } else {
    return(Xmat)
  }
}


#' Generate residual from a mutivariate normal distribution
#' 
#' Simulate from a zero-mean Multivariate Normal Distribution
#' 
#' Generate residual data from a Multivariate normal distribution.
#' 
#' @param n sample size
#' @param p dimension
#' @param rho \eqn{S_{p*p} = \rho^{|i-j|}} the element in the covariance matrix,
#'   default \eqn{\rho=1}
#' @return n*p residual matrix
gen_MultNorm_residual <- function(n, p, rho=1) {
  # library(MASS)
  # Sigma <- matrix(c(10,3,3,2),2,2)
  # Sigma
  # var(mvrnorm(n = 1000, rep(0, 2), Sigma))
  # var(mvrnorm(n = 1000, rep(0, 2), Sigma, empirical = TRUE))
  S <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      S[i, j] <- rho^{abs(i - j)}
    }
  }
  eps <- MASS::mvrnorm(n, rep(0, p), Sigma=S)
  return(eps)
}



#' Generate high-dimension CV Multivariate Normal data
#' 
#' Generate Multivariate Normal data for high-dimension CV 
#' 
#' Assume diff_signal = FALSE forever; and amplify_signal=1 forever
#' 
#' @inheritParams gen_hcv_ARMA_data
#' @inheritParams gen_MultNorm_residual
#' @examples 
#'   tau2n <- c(0.2, 0.4, 0.6, 0.8)
#'   n <- 400; p <- 100; r <- 0.1; umax <- 1; rho <- 1
#'   dt <- gen_hcv_MvNorm_data(
#'       n, p, r, umax, tau2n, rho, diff_signal=TRUE, verbose=TRUE)
#'   Xu <- dt$Xu
gen_hcv_MvNorm_data <- function(
    n, p, r, umax, tau2n, rho, diff_signal=FALSE, verbose=FALSE){
  amplify_signal <- 1
  Kn <- length(tau2n)
  hmat <- genconst_mat_signal(Kn, umax, p, r, diff_signal)
  
  # ## generate the signal -- Xu matrix
  # tau_true <- floor(tau2n * n)
  # tau_vec <- c(0, tau_true, n)  # 0 and n are included.
  # Xu <- vector()
  # for (i in 1:(Kn + 1)) {
  #   u_i <- hmat[i, ]
  #   ni <- tau_vec[i+1] - tau_vec[i]
  #   mu_i <- matrix(rep(u_i, ni), ni, p, byrow = TRUE)
  #   Xu <- rbind(Xu, mu_i)
  # }
  
  ## generate the signal -- Xu matrix
  tau_true <- floor(tau2n * n)
  taut <- c(0, tau_true, n)  # 0 and n are included.
  Theta <- matrix(0, n, 0)
  Kn <- length(tau_true)
  for (j in 1:(Kn + 1)){
    temp <- c(rep(0, taut[j]), rep(1, taut[j+1] - taut[j]), rep(0, n - taut[j+1]))
    Theta <- cbind(Theta, temp)
  }
  Xu <- Theta %*% hmat
  
  eps <- gen_MultNorm_residual(n, p, rho=rho)
  Xmat <- Xu + eps
  if (verbose){
    return(list(Xmat=Xmat, Xu=Xu, eps=eps))
    # A = list(Xmat=Xmat, Xu=Xu, eps=eps)
  } else {
    return(Xmat)
  }
}
