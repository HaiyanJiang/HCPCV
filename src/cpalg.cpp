#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#define INFINITYNUMBER 1E7

// [[Rcpp::plugins(cpp11)]]
// // Add before some functions we want to use c++11

/* In the R file, before we install the package, we need to specify a compiler.
 require(Rcpp)
 Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
 system('g++ -v')
 Sys.which("gcc.exe")
 Sys.getenv()['PATH']
 */

// changepoint detection algprithms:
// DP (dynamic programming) algorithm, 
// OP (optimal partition),
// PELT (pruned exact linear time) method.


// [[Rcpp::export]]
void introduction_of_cpp_index() {
    std::string myintro = "\n                                                    \
// It is annoying to deal with the cpp_index and the normal R_index.\n           \
// \n                                                                            \
// One trick is to first make sure every cpp_index is in the proper range. \n    \
// A convenient way is first to check the bounds of the cpp_index, \n            \
// second to focus only on the total number of data points we considered\n       \
//   and ignore what the last index is.\n                                        \
// \n                                                                            \
// \n                                                                            \
// DP:\n                                                                         \
// cpp index:\n                                                                  \
//   H(Kmax+1, n) keep records of candidate taus, tau_0=0, tau_{Kmax+1}=n,\n     \
//      the j-th segment is (tau_{j-1}, tau_{j}], which means the left point\n   \
//      is not included, but the right one is. H(*, s), s=0,...,n-1.\n           \
// F(Kmax+1, n), keep records of optimal loss, proper s index is [0,..., n-1].\n \
//      F(*, s) means until the index of s, [0,...,s], so there are (s+1)\n      \
//      data points included, the number of data points is (s+1).\n              \
// \n                                                                            \
// Paper notation: F(k, s) = F(k-1, t) + cost(t+1, s); with changepoint t.\n     \
// Paper index: [1,...,s], [1,...,t], [t+1,...,s]\n                              \
// Code notation: F(k, s) = F(k-1, t) + cost(t+1, s); s = 0,...,n-1\n            \
// Code index: [0,...,s], [0,...,t], [t+1,...,s].\n                              \
// s can take values from [0,...,n-1], when s=0, means 1 data, s=n-1, means n.\n \
// \n                                                                            \
// The number of data points considered in the first part is (t+1), \n           \
// then (t+1) is some candidate changepoint, tau_{s} = t+1, so H[s] = t + 1.\n   \
// \n                                                                            \
// OP and PELT: \n                                                               \
//   initialized F(n+1) and H(n+1).\n                                            \
// \n                                                                            \
// F[s] means s data points included, [0,...,s-1]. s = [1,...,n]\n               \
// Paper notation: F[s] = F[t] + cost(t+1, s); with changepoint t.\n             \
// Paper index: [1,...,s], [1,...,t], [t+1,...,s]\n                              \
// Code notation: F[s] = F[t] + cost(t, s-1); s = 1, ..., n\n                    \
// Code index: [0,...,s-1], [0,...,t-1], [t,...,s-1]\n                           \
// \n                                                                            \
// The number of data points considered in the first part is (t), \n             \
// then (t) is some candidate changepoint, t is already the R index, \n          \
// tau_{s} = t, so H[s] = t.\n ";
  cout << myintro << endl;
}


// It is annoying to deal with the cpp_index and the normal R_index.
// 
// One trick is to first make sure every cpp_index is in the proper range. 
// A convenient way is first to check the bounds of the cpp_index, 
// second to focus only on the total number of data points we considered
//   and ignore what the last index is.
// 
// 
// DP:
// cpp index:
//   H(Kmax+1, n) keep records of candidate taus, tau_0=0, tau_{Kmax+1}=n,
//      the j-th segment is (tau_{j-1}, tau_{j}], which means the left point
//      is not included, but the right one is. H(*, s), s=0,...,n-1.
// F(Kmax+1, n), keep records of optimal loss, proper s index is [0,..., n-1].
//      F(*, s) means until the index of s, [0,...,s], so there are (s+1)
//      data points included, the number of data points is (s+1).
// 
// Paper notation: F(k, s) = F(k-1, t) + cost(t+1, s); with changepoint t.
// Paper index: [1,...,s], [1,...,t], [t+1,...,s]
// Code notation: F(k, s) = F(k-1, t) + cost(t+1, s); s = 0,...,n-1
// Code index: [0,...,s], [0,...,t], [t+1,...,s].
// s can take values from [0,...,n-1], when s=0, means 1 data, s=n-1, means n.
// 
// The number of data points considered in the first part is (t+1), 
// then (t+1) is some candidate changepoint, tau_{s} = t+1, so H[s] = t + 1.
// 
// OP and PELT: 
//   initialized F(n+1) and H(n+1).
// 
// F[s] means s data points included, [0,...,s-1]. s = [1,...,n]
// Paper notation: F[s] = F[t] + cost(t+1, s); with changepoint t.
// Paper index: [1,...,s], [1,...,t], [t+1,...,s]
// Code notation: F[s] = F[t] + cost(t, s-1); s = 1, ..., n
// Code index: [0,...,s-1], [0,...,t-1], [t,...,s-1]
// 
// The number of data points considered in the first part is (t), 
// then (t) is some candidate changepoint, t is already the R index, 
// tau_{s} = t, so H[s] = t.






// [[Rcpp::export]]
NumericVector estimated_sigma_square(NumericMatrix Y) {
  // Return the estimated variance sigma_hat^2 vector.
  // Y is an n*p matrix, return p-length vector nor(p).
  // nor[j] = 1/(2*(n - 1)) * sum_{i=1}^{n-1} (Y[i+1,j] - Y[i,j])^2.
  int n = Y.nrow();
  int p = Y.ncol();
  NumericVector nor(p);
  for (int j = 0; j < p; j++) {
    nor[j] = 0;
    for (int i = 0; i < n - 1; i++){
      nor[j] = nor[j] + pow(Y(i + 1,j) - Y(i,j), 2);
    }
    nor[j] = nor[j]/(2*(n - 1));
  }
  return nor;
}


// [[Rcpp::export]]
NumericMatrix standardize_matrix(NumericMatrix X, NumericVector D) {
  // X is a n*p matrix, D is a p-length vector.
  // Assume E(X_i) = 0, Var(X_i) = E[(X_i)^2] = D = Sigma
  // Y_i = D^(-1/2) * X_i. then E[(Y_i)^2] = I_{p}.
  // Y_{n*p} = X_{n*p} * diag(D)_{p*p}^{-1/2}.
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix Y(n, p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      Y(i, j) = X(i, j) / pow(D[j], 0.5); // the inverse of D sqrt.
    }
  }
  return Y;
}


// [[Rcpp::export]]
NumericMatrix normalize_matrix(NumericMatrix X, NumericVector W, bool normalized = false) {
  // X is a n*p matrix, D is a p-length vector.
  // Yi = D^(-1/2) * Xi. Y_{n*p} = X_{n*p} * diag(D)_{p*p}^{-1/2}.
  // W = D^(-1), then D^(-1/2) is W^(1/2)
  if (not normalized) {
    return X;
  }
  
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix Y(n, p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      Y(i, j) = X(i, j) * pow(W[j], 0.5); // the sqrt of D.
    }
  }
  return Y;
}


// [[Rcpp::export]]
double sum_knorm(NumericVector x, int k) {
  // sum_{i=1}^{n} xi^k;
  int n = x.size();
  double sk = 0;
  for (int i = 0; i < n; i++) {
    sk += pow(x[i], k);
  }
  return sk;
}


// [[Rcpp::export]]
NumericMatrix cost_matrix(NumericMatrix X, NumericVector W) {
  // cost_matrix and cost_matrix_recursion_v1 are the same functions.
  // The recursive version of calculating the cost matrix function. 
  // ONLY recursion, use the relationship between C(t,s) and C(t-1, s).
  // C(t,s) = C(t-1, s) - \frac{s-t+1}{s-t+2} ( Y_{t-1} - \bar{Y}_{[t,s]})^2
  // WITHOUT storing any cumulative sum MATRIX of the Y matrix.
  // ONE p-dim vector to store the cumulative sum of Y[0:s], cums[0:s]
  // ONE p-dim temp vector to store the sum from Y[t:s], cts
  // cts_{t, s}[j] = cums_{0:s}[j] - Y(t, j)
  // The t indicates the initial date, while s indicates the terminal date.
  // We need to calculate each value of C(t, s);
  
  NumericMatrix Y = normalize_matrix(X, W);
  
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix costmat(n, n);  // the return results.
  double sq = 0;  // which is used to record the sum_{i=0}^{s} ||Y_{i}||^2.
  NumericVector cums(p);  // to record cumx[j] = sum_{i=0}^{s} Y_{ij}
  for (int j = 0; j < p; j++) {
    cums[j] = 0;
  }
  for (int s = 0; s < n; s++) {
    for (int j = 0; j < p; j++) {
      sq += pow(Y(s, j), 2);
      cums[j] += Y(s, j);
    }
    
    costmat(0, s) = sq - sum_knorm(cums, 2) / (s + 1);  // (s+1) terms.
    
    NumericVector cts(p);  // cumulative sum of [t:s] after updating,
    // before updating is sum of [t+1:s], because updating is substraction.
    // need to be initialized inside each s loop.
    for (int j = 0; j < p; j++) {
      cts[j] = cums[j] - Y(0, j);
    }
    for (int t = 1; t <= s - 1; t++) {
      double tmp = 0;  // the temp squres
      for (int j = 0; j < p; j++) {
        tmp += pow(Y(t - 1, j) - cts[j] / (s - t + 1), 2);
        cts[j] -= Y(t, j);
      }
      costmat(t, s) = costmat(t - 1, s) - tmp * (s - t + 1) / (s - t + 2); 
    }
  }
  return costmat;
}




// [[Rcpp::export]]
NumericVector seek_dp_tau(int L, NumericMatrix H) {
  // L: given L changepoints, 1 <= L <= Kmax
  // H: (Kmax +1)*n matrix, the optimal locations of the estimated tau.
  // H(L, ) stores tauHat_{1}, tauHat_{L}, the locations
  // (L+1)-th row of H, with row index L, with L-change-points.
  int n = H.ncol();
  NumericVector tau(L);  // 0, ..., k-1
  
  tau[L - 1] = H(L, n - 1);
  for (int i = L - 1; i > 0; i--) {
    tau[i - 1] =  H(i, tau[i] - 1);  // KEY POINTS: tau[i] - 1!!!!
    // the previous index is i-1.
  }
  // for (int i = L - 2; i >= 0; i--){
  //   tau[i] = H(i+1, tau[i+1] - 1);
  // }
  return tau;
}


// NOTE: here we cannot use // [[Rcpp::export]] to exptort this function
// because, Rcpp do not have pointer* nor reference &, so this function can
// only be used inside other cpp functions, not in R;
void min_which(double *array, int n, double &minout, int &whichout) {
  // Function to find minimum of an array with n elements that is put in min 
  minout = *array;
  whichout = 0;
  int i;
  for(i = 1; i < n; i++) {
    if(*(array+i) < minout){
      minout = *(array+i); 
      whichout = i;
    }
  }
}


// [[Rcpp::export]]
void findMin(NumericVector arr, int n, double& minVal, int& minIdx){
  /* We are assigning the first array element to the minVal variable and
   * then we are comparing all the array elements with the minVal inside
   * loop and if the element is smaller than minVal
   * then the minVal value is replaced by that.
   * This way we always have the smallest value in minVal.
   * Finally we are reference minVal.
   */
  // location is (cpp_index + 1)
  // so i is the cpp index, and the location is (i+1)
  minVal = arr[0];
  minIdx = 0;
  for (int i = 0; i < n; i++) {
    if(minVal > arr[i]) {
      minVal = arr[i];
      minIdx = i;
    }
  }
}



//' The Dynamic Programming Method
//' 
//' Change-point detection algorithm for the high dimension mean change model
//' using the dynamic programming.
//' 
//' Dynamic Programming of changepoint detection in order to get all 
//' possible changepoints given a maximum number of changepoints, \code{Kmax}.
//' 
//' @param X a n*p data matrix
//' @param W the normalized scale of the data
//' @param Kmax scalar, the maximum number of changepoints
//' @examples 
//' n <- 100
//' tau2n <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
//' tau2n <- c(0.1, 0.40, 0.65, 0.8)
//' (tau_true <- floor(tau2n * n))
//' x <- gen_piecewise_const_vector(n, tau_true)
//' X <- matrix(x, n, 1)
//' # W <- estimated_sigma_square(X)
//' W <- rep(1, ncol(X))
//' Kmax <- length(tau_true)
//' tau_mat <- cp_hdmean_dp(X, W, Kmax)
//' tau_mat[Kmax, ]
// [[Rcpp::export]]
NumericMatrix cp_hdmean_dp(NumericMatrix X, NumericVector W, int Kmax) {
  // change-point, high dimension mean, dynamic programming, cpp codes
  // X: the n*p data matrix. 
  // Kmax: scalar, the possible maximum number of change-points.
  // Here, s and t are the proper cpp index taking value from [0,...,n-1].
  // Floss(Kmax+1, n): k = [0,..., Kmax], each row represents with k number
  // change-points, F(k, s) the optimal loss with k change-points of y[0:s].
  // H(Kmax + 1, n): the candidate tau's locations (location=cpp_index + 1).
  // tau_mat(Kmax, Kmax); row k is a k-length vector, with k changepoints.
  // Floss(k, s) = min_{0 <= t <= s-1} Floss(k-1, t) + cost(t+1, s), where 
  // s and t are the proper cpp index taking value from [0,...,n-1].
  // IMPORTANT: Floss, H, initialize with all 0's! 
  // The initialized value can be different values, other than 0's, because
  // We ONLY use the upper matrix of Floss and H, t is from k.
  int n = X.nrow();
  
  NumericMatrix Floss(Kmax + 1, n);
  NumericMatrix H(Kmax + 1, n);  // Floss, H, cost_matrix are with same n.
  // Floss, H, initialize with all 0's!
  
  // cost_matrix is n*n matrix.
  NumericMatrix cost = cost_matrix(X, W);  // index [0,...,n-1]
  
  for (int j = 0; j < n; j++) {
    Floss(0, j) = cost(0, j);
  }
  // i means i changepoint
  // j is the cpp index, maximum (n-1)
  // t is from i-1, is the cpp index.
  for (int i = 1; i < Kmax + 1; i++) {
    for (int j = i; j < n; j++) {
      NumericVector arr(j-i+1);
      for (int t = i - 1; t <= j - 1; t++) {
        arr[t-i+1] = Floss(i - 1, t) + cost(t + 1, j);
      }
      double fmin;
      int tmin;
      findMin(arr, j-i+1, fmin, tmin);
      
      Floss(i, j) = fmin;
      H(i, j) = tmin + i;
      // tmin + (i-1) + 1; (i-1) is that t starts from i-1; 
      // 1 is the gap between cpp_index and location
      // location is (cpp_index + 1)
      // so tmin is the cpp index, so the location is (tmin+1)
    }
  }
  
  // backwards searching for the locations of the changepoints.
  // tau_mat(i, ) is i-len vecotr with cp locations - tau's.
  NumericMatrix tauMat(Kmax, Kmax);
  for (int i = 1; i <= Kmax; i++) {
    NumericVector tau = seek_dp_tau(i, H);  // here i. i-len vector.
    for (int j = 0; j < i; j++) {
      tauMat(i - 1, j) = tau[j];  // the index of tauMat is i-1 not i.
      // because i starts with 1 and ends with Kmax.
    }
  }
  List A = List::create(
    Named("Floss") = wrap(Floss), Named("Hindex") = wrap(H), 
    Named("cptaumat") = wrap(tauMat));
  return tauMat;
}





// // [[Rcpp::export]]
NumericMatrix remove_cpp_rows_of_matrix(NumericMatrix X, NumericVector a) {
  // a: the cpp row index that we want to remove from X, [0, ..., n1-1]
  int p = X.ncol();
  int n1 = X.nrow();
  set<int> s;
  for (int i = 0; i < a.size(); i++) {
    s.insert(a[i]);
  }
  int n2 = s.size();
  int n = n1 - n2;
  NumericMatrix Y(n, p);
  int nless = 0;
  for (int i = 0; i < n; i++) {
    int nx = i + nless;
    while (s.count(nx)) {
      nless += 1;
      nx = i + nless;
    }
    
    for (int j = 0; j < p; j ++) {
      Y(i, j) = X(nx, j);
    }
  }
  return Y;
}


// [[Rcpp::export]]
NumericMatrix remove_rows_of_matrix(NumericMatrix X, NumericVector a) {
  // a: the normal row index (not cpp index) that we want to remove from X.
  // a can have values from [1, ..., n1]
  int p = X.ncol();
  int n1 = X.nrow();
  set<int> s;
  for (int i = 0; i < a.size(); i++) {
    s.insert(a[i] - 1);  // here we need to transform it to the cpp index.
  }
  int n2 = s.size();
  int n = n1 - n2;
  NumericMatrix Y(n, p);
  int nless = 0;
  for (int i = 0; i < n; i++) {
    int nx = i + nless;
    while (s.count(nx)) {
      nless += 1;
      nx = i + nless;
    }
    
    for (int j = 0; j < p; j ++) {
      Y(i, j) = X(nx, j);
    }
  }
  return Y;
}


// [[Rcpp::export]]
NumericVector init_vector_with_constant(int len, int c) {
  // [c, c+1, ..., c+len-1]
  NumericVector vec(len); 
  for (int i = 0; i < len; i++) {
    vec[i] = c + i;
  }
  return vec;
}

//' Calculate the trace of \eqn{R^2} using Rcpp codes
//' 
//' Calculate the trace of \eqn{R^2}, \eqn{R^2} is correlation matrix of data X
//' 
//' It is slower compared to the same one but using R to calculate the value.
//' If we use the penalty from Yunlong Wang paper, which is calculated
//' \deqn{ penalty = p + c_0 * (trace(R^2))^{1/2} * log(n)^{1+\alpha}}
//' With 
//' \deqn{
//' trace(R^2) = 1 /(4(n-3)) \sum_{i=1}^{n-3}
//'   Z_{i}^T * (D_{-(i,i+1,i+2,i+3)})^(-1) * Z_{i+2}
//' }
//' where \eqn{ Z_{i} = X_{[i,]} - X_{[i+1,]} } }
//' 
//' @param X n*p data matrix
//' @return the value of trace(\eqn{R^2})
//' 
//' @examples
//' X = matrix(1:100, 20, 5)
//' trace_R_square_cpp(X)
//' 
// [[Rcpp::export]]
double trace_R_square_cpp(NumericMatrix X) {
  // function used to calculate the penalty from Yunlong Wang's paper.
  // pen = p + 1.5 * trace(R^2)^{1/2} * log(n)^{1.1}
  // traceR2 = 1/(4*(n-3)) \sum_{i=1}^{n-3}
  // (X_[i,] - X_[i+1,]) * (D_{-(i,i+1,i+2,i+3)})^(-1) * (X_[i+2,] - X_[i+3,])
  int n = X.nrow();
  int p = X.ncol();
  double s = 0;
  int len = 4;
  
  for (int i = 0; i < n - 3; i++) {
    NumericVector a = init_vector_with_constant(len, i);
    NumericMatrix Xtemp = remove_cpp_rows_of_matrix(X, a);
    NumericVector D = estimated_sigma_square(Xtemp);
    double sij = 0;
    for (int j = 0; j < p; j++) {
      sij += ((X(i, j) - X(i + 1, j)) * (X(i + 2, j) - X(i + 3, j)) / D[j]);
    }
    s += pow(sij, 2);
  }
  return (s / (4 * (n - 3)));
}


//' The Optimal Partitioning Method of High Dimension Mean Change Model
//' 
//' Optimal Partitioning method of changepoint detection.
//' 
//' Change-point detection algorithm for the high dimension mean change model
//' using Optimal Partitioning.
//' 
//' @inheritParams cp_hdmean_dp
//' @param pen the penalty of the cost function
//' @return the locations of estimated change-points
//' 
//' @examples
//' n <- 100
//' tau2n <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
//' tau2n <- c(0.1, 0.40, 0.65, 0.8)
//' (tau_true <- floor(tau2n * n))
//' x <- gen_piecewise_const_vector(n, tau_true)
//' X <- matrix(x, n, 1)
//' W <- estimated_sigma_square(X)
//' # W <- rep(1, ncol(X))
//' pen <- log(n)^1.1
//' op_list <- cp_hdmean_op(X, W, pen)
//' (op_list$cptau)
//' 
// [[Rcpp::export]]
List cp_hdmean_op(NumericMatrix X, NumericVector W, double pen) {
  // change-point, high dimension mean, optimal partitioning, cpp codes.
  // X: the n*p data matrix.
  // pen: the penalty of the objective function.
  // Floss(n+1): vector, F[s] the optimal loss with s data points,
  // with cpp index, it would be y[0:(s-1)],
  // since Floss[0] is occupied by the negative penalty. 
  // Floss[s], with length s, means s data points of y_i's, [0,..,s-1].
  // H(n+1): H[s], tau's index (location=index+1) corresponding to Floss[s].
  // tau_hat(ncp); vector, with ncp locations (index+1) of the changepoints.
  // NOTES: Floss, H with (n+1) elements, while cost_matrix with n.
  int n = X.nrow();
  NumericVector Floss(n + 1);  // Floss (n+1) elements, cost_matrix with n.
  NumericVector H(n + 1);  // (n+1) elements, H[s] to be accordant to Floss[s].
  for (int i = 0; i < n + 1; i++) {H[i] = -6;}  // initialize to be (-6)
  // H[i] = 0 means there is no change-point for all possible cuts of data
  // y[0:(i-1)].
  
  NumericMatrix cost = cost_matrix(X, W);  // index [0,...,n-1]
  
  Floss[0] = - pen;
  
  for (int s = 1; s < n + 1; s++) { // s maximum n, outside cost_matrix index
    double fmin = INFINITYNUMBER;
    int tmin = - INFINITYNUMBER;
    for (int t = 0; t <= s - 1; t++) {
      double fcur = Floss[t] + cost(t, s - 1) + pen;
      // Theoretically, F(s) = min_{0<=t<=s-1} F(t) + cost(t+1, s).
      // F(s) means with total s points, the optimal loss value.
      // s maximum n, so s is one ahead compared to the cpp index pf cost(t, s).
      // Both t and s are one index ahead. If t is not ahead, it would be
      // Floss[t] + cost(t+1, s-1), when s=n and t=s-1, then t+1=s=n, 
      // index out of bounds, index of cost(i, j), i, j is [0, ..., n-1].
      // Floss[t] is related to data y[0:(t-1)], with t data points.
      // R: [t+1:s], cpp: [t:s-1] because s starts with 1 ends with n.
      // attention to the index of cost matrix, here is (t, s-1) not (t+1,s)
      if (fcur < fmin) {
        fmin = fcur;
        tmin = t;  // t is the cpp index, not the location which is (t+1)
        // check the bounds: when s=n and t=s-1, then t+1=n, which means 
        // (t+1) is out of cpp index bounds, so (t+1) is the location.
        // t can take values from [0,...,n-1].
      }
    }
    Floss[s] = fmin;
    H[s] = tmin;  // H, with (n + 1) elements.
    
  }
  
  // backwards searching for the locations of the changepoints. tau_candidate
  NumericVector tau_cand(n);  // initialized value are zeros.
  // for (int i = 0; i < n; i++) { tau_cand[i] = -8;}
  int ncp = 0;  // the number of changepoints. H is [0,..., n-1, n], 
  int t_cur = H[n];  // H with (n+1) elements, the last proper index is n.
  while (t_cur > 0) {
    // cout << t_cur << ' ' << endl;  // used to check if the codes right.
    tau_cand[ncp] = t_cur;
    ncp += 1;
    t_cur = H[t_cur - 1];  // the index moves one step ahead.
    // Moving one index ahead is right, but I do not know why.
    // Why should the index move one step ahead?
    // For example: n = 100, H[100] = 60; H[60] = 10; H[10] = 0;
    // F[100] = F[60] + cost[61:100], H[100] = 60
    // F[s] = min_{0 <= t <= s-1} F[t] + cost[t+1:s]
    // 60 is the index, actually it is already the changepoint.
    // cost(t, s), the t is included in the last part.
    // (1, 10], (10, 60], (60, 100]
  }
  NumericVector tau_hat(ncp);  // only keep the useful part.
  for (int i = 0; i < ncp; i++) {
    tau_hat[i] = tau_cand[ncp - 1 - i];
  }
  
  List A = List::create(
    Named("Floss") = wrap(Floss), Named("Hindex") = wrap(H), 
    Named("cpnum") = ncp, Named("cptau") = wrap(tau_hat));
  return A;
  // return tau_hat;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List cp_hdmean_pelt(NumericMatrix X, NumericVector W, double pen) {
  // change-point, high dimension mean, PELT, cpp codes
  // X: the n*p data matrix. 
  // pen: the penalty of the objective function.
  // Floss(n+1): vector, F[s] the optimal loss with s data points y[0:(s-1)],
  // since Floss[0] is occupied by the negative penalty. 
  // Floss[s], with length s, means s data points of y_i's, [0,..,s-1].
  // H(n+1): H[s], the tau's index (location=index+1) corresponding to Floss[s].
  // tau_hat(ncp); vector, with ncp locations (index+1) of the changepoints.
  // R: a set, the remaining set of index (not location but index).
  
  int n = X.nrow();
  NumericVector Floss(n + 1);  // Floss, H with (n+1), cost_matrix with n.
  NumericVector H(n + 1);  // (n+1) elements, H[s] to be accordant to Floss[s].
  for (int i = 0; i < n + 1; i++) {H[i] = -6;}  // initialize to be (-6)
  // H[i] = 0 means there is no changepoint for all possible cut of the data.
  set<int> R = {0};    // the remaining set of index;
  int nerase = 0;
  
  NumericMatrix cost = cost_matrix(X, W);  // index [0,...,n-1]
  
  Floss[0] = - pen;
  double K = 0;  // K is used to check if the current t is pruned or not.
  
  for (int s = 1; s < n + 1; s++) { // s maximum n, outside cost_matrix index
    
    int ns = R.size();
    NumericVector Ftemp(ns);
    NumericVector Htemp(ns);
    int i = 0;
    
    double fmin = INFINITYNUMBER;
    int tmin = -5;
    for (set<int>::iterator it = R.begin(); it != R.end(); it++) {
      // if (s == 2) { cout << ' ' << *it << endl; }
      double loss = Floss[*it] + cost(*it, s - 1);
      // s is one ahead compared to the cpp index of cost(t, s). 
      // Both t and s are one index ahead. If t is not ahead, it would be
      // Floss[t] + cost(t+1, s-1), when s=n, t=s-1, then t+1=n, 
      // index out of bounds, index of cost(i, j), i, j is [0, ..., n-1]
      // Floss[t] is related to data y[0:(t-1)]
      // R: [t+1:s], cpp: [t:s-1]
      // attention to the index of cost matrix, here is (t, s-1) not (t,s)
      double fcur = loss + pen;
      if (fcur < fmin) {
        fmin = fcur;
        tmin = *it;  // keep records of the cpp index of the data.
      }
      // keep records of the vector related to the pruning.
      Ftemp[i] = loss + K;
      Htemp[i] = *it;  // cpp index of the data.
      i += 1;
    }
    Floss[s] = fmin;
    H[s] = tmin;  // H, with (n+1) elements. 
    // H[s] with s data points, index is [0:s-1]
    
    // keep records of the Ftemp and Htemp, not re-compute
    for (int i = 0; i < ns; i++) {
      /*
       if (s == 2) {
       cout << "The loss: " << Ftemp[i] << endl;
       cout << "The index: " << Htemp[i] << endl;
       cout << "The optimal loss: " << Floss[s] << endl;
       if (Ftemp[i] >= Floss[s]) {
       cout << "pruning: " << Htemp[i] << endl;
       }
       }
       */
      if (Ftemp[i] >= Floss[s]) {
        R.erase(Htemp[i]);
        nerase += 1;
      }
    }
    R.insert(s);
  }
  
  // backwards searching for the locations of the changepoints. tau_candidate
  NumericVector tau_cand(n);  // initialized value are zeros.
  // for (int i = 0; i < n; i++) { tau_cand[i] = -8;}
  int ncp = 0;  // the number of changepoints. H is [0,..., n-1, n], 
  int t_cur = H[n];  // H with (n+1) elements, the last proper index is n.
  while (t_cur > 0) {
    // cout << t_cur << ' ' << endl;  // used to check if the codes right.
    tau_cand[ncp] = t_cur;
    ncp += 1;
    t_cur = H[t_cur - 1];  // the index moves one step ahead.
    // Moving one index ahead is right, but I do not know why.
    // Why should the index move one step ahead?
    // For example: n = 100, H[100] = 60; H[60] = 10; H[10] = 0;
    // F[100] = F[60] + cost[61:100], H[100] = 60
    // F[s] = min_{0 <= t <= s-1} F[t] + cost[t+1:s]
    // 60 is the index, actually it is already the changepoint.
    // cost(t, s), the t is included in the last part.
    // (1, 10], (10, 60], (60, 100]
  }
  
  NumericVector tau_hat(ncp);  // only keep the useful part.
  for (int i = 0; i < ncp; i++) {
    tau_hat[i] = tau_cand[ncp - 1 - i];
  }
  
  List A = List::create(
    Named("Floss") = wrap(Floss), Named("Hindex") = wrap(H),
    Named("cpnum") = ncp, Named("cptau") = wrap(tau_hat),
    Named("nonprunedset") = wrap(R), Named("rsize") = nerase);
  return A;
  // return tau_hat;
}


//' \eqn{C(\tau^{train}; Xtest)}
//' 
//' The function of definition of \eqn{C(\tau^{train}; Xtest)}
//' 
//' \deqn{C(\tau^{train}; Xtest) = \sum_{l=0}^{L}
//'   \sum_{i=\tau_{l}+1}^{\tau{l+1}}
//'   (Xtrain_i - \bar(Xtest)_{\tau_{l}, \tau{l+1}})^{2} / W }
//' 
//' @param Xtrain, matrix with odd index, which is used to estimate tau.
//' @param Xtest, matrix with even index, used to construct the prediction error.
//' @param W, matrix with p*p, as the global normalizer, by dividide W[j]
//' @param tau, the vector of change-points, including 0, and n. 
//'   Usually tau is obtained by Xtrian, and this cv_objfun is used to 
//'   calculate the cross-validation squared loss.
//' @return mean of the cross validation error
// [[Rcpp::export]]
double cv_objfun(
    NumericMatrix Xtrain, NumericMatrix Xtest, 
    NumericVector W, NumericVector tau) {
  // High dimension CV objective function, mean change-point model, cpp codes.
  // Xtrain: X matrix with odd index, which is used to train/estimate/get tau.
  // Xtest: X matrix with even index, used to construct the prediction error.
  // W: vector with p-length as the global normalizer.
  // tau: tau_0=0, tau_{1}, tau_{2}, ..., tau_{L}, tau_{L+1}=n.
  // sum_{k=1}^{L+1} \sum_{tau_{k-1}+1}^{tau_{k}}
  // sum_{k=0}^{L} \sum_{tau_{k}+1}^{tau_{k+1}}
  // sum_{j=1}^{p}(Xtest(i, j) - mo[j])^2 / W[j]
  int p = Xtrain.ncol();
  int L = tau.size() - 2;  // ncp = tau.size(); ncp = L + 2
  // It also applies when L = 0
  
  double obj = 0;
  for (int k = 0; k < L + 1; k++) {
    int a = tau[k] + 1;
    int b = tau[k + 1];
    // R index: [tau_{k} + 1, tau_{k+1}]
    // cpp index: [tau_{k}, tau_{k+1} - 1], a-1, b-1
    // To get the mean of the Xtrain in each segment, with (b-a+1) length
    NumericVector mo(p);  // mean of Xtrain, p-len, keep track of segment-mean
    for (int j = 0; j < p; j++) {
      mo[j] = 0;
      for (int i = a - 1; i < b; i++) {mo[j] = mo[j] + Xtrain(i, j);}
      mo[j] /= (b - a + 1);
    }
    
    double s = 0;
    for (int i = a - 1; i < b; i++) {
      for (int j = 0; j < p; j++) {
        s += pow(((Xtest(i, j) - mo[j]) / pow(W[j], 0.5)), 2);
      }
    }
    obj += s;
  }
  // we do not add nor minus L*p, obj is only the squared loss.
  return obj;
}


//' \eqn{S_{xy}}, the function of \eqn{S_{xy}}
//' 
//' The function of definition of \eqn{S_{xy}}
//' 
//' \deqn{ S_{xy}(\tau(L); W) = \sum_{l=0}^{L}
//'   \sum_{i=\tau_{l}+1}^{\tau{l+1}}
//'   (x_i - \bar(x)_{\tau_{l}, \tau{l+1}})^{\top} W 
//'   (y_i - \bar(y)_{\tau_{l}, \tau{l+1}}) }
//' 
//' @inheritParams cv_objfun
//' @param X data matrix X
//' @param Y data matrix Y
//' 
//' @return the value of Sxy
//' 
// [[Rcpp::export]]
double Sxy_fun(
  NumericMatrix X, NumericMatrix Y, NumericVector W, NumericVector tau) {
  // High dimension CV objective function, mean change-point model, cpp codes.
  // X: X matrix with n*p
  // Y: Y matrix with n*p
  // W: vector with p-length as the global normalizer.
  // tau: some candidate locations of change-points 
  //   tau_0=0, tau_{1}, tau_{2}, ..., tau_{L}, tau_{L+1}=n
  // sum_{k=1}^{L+1} \sum_{tau_{k-1}+1}^{tau_{k}}
  // sum_{k=0}^{L} \sum_{tau_{k}+1}^{tau_{k+1}}
  // sum_{j=1}^{p}(Xe(i, j) - mo[j])^2 / W[j]
  int p = X.ncol();
  int L = tau.size() - 2;  // ncp = tau.size(); ncp = L + 2
  // It also applies when L = 0
  
  double obj = 0;
  for (int k = 0; k < L + 1; k++) {
    int a = tau[k] + 1;
    int b = tau[k + 1];
    // R index: [tau_{k} + 1, tau_{k+1}]
    // cpp index: [tau_{k}, tau_{k+1} - 1], a-1, b-1
    // To get the mean of the Xo in each segment, with (b-a+1) length
    NumericVector xm(p);  // mean of X, p-len, keep track of segment-mean
    NumericVector ym(p);  // mean of Y, p-len, keep track of segment-mean
    for (int j = 0; j < p; j++) {
      xm[j] = 0;
      ym[j] = 0;
      for (int i = a - 1; i < b; i++) {
        xm[j] = xm[j] + X(i, j);
        ym[j] = ym[j] + Y(i, j);
      }
      xm[j] /= (b - a + 1);
      ym[j] /= (b - a + 1);
    }
    
    double s = 0;
    for (int i = a - 1; i < b; i++) {
      for (int j = 0; j < p; j++) {
        s += ((X(i, j) - xm[j]) * (Y(i, j) - ym[j]) / W[j]);
      }
    }
    obj += s;
  }
  // we do not add nor minus L*p, obj is only the squared loss.
  return obj;
}



// [[Rcpp::export]]
NumericMatrix test_matrix_plus_scalar(NumericMatrix X, double a) {
  return X + a;
}



// [[Rcpp::export]]
void test_pen_char_param(NumericMatrix X, std::string penaltyform = "BIC") {
  // BIC: nor*(log(n))^1.1
  // Zou: nor*(log(n))^2.1/2;
  int n = X.nrow();
  NumericVector D = estimated_sigma_square(X);
  int p = D.size();
  double pen_default = log(n);
  double pen = 0;
  if (penaltyform == "BIC") {
    for (int j = 0; j < p; j++) {
      pen += pow(pen_default, 1.1) * D[j];
    }
    cout << "We can use the 'string' type with double quots!" << endl;
    cout << "BIC penalty: " << pen << endl;
  } else if (penaltyform == "ZOU"){
    for (int j = 0; j < p; j++) {
      pen += pow(pen_default, 2.1) * D[j] / 2;
    }
    cout << "ZOU penalty: " << pen << endl;
    cout << "lost!" << endl;
  } else if (penaltyform == "WANG") {
    // the penalty from Yunlong Wang's paper.
    // pen = p + 1.5 * trace(R^2)^{1/2} * log(n)^{1.1}
    pen = p + 1.5 * pow(trace_R_square_cpp(X), 0.5) * pow(log(n), 1.1);
    cout << "Wang penalty: " << pen << endl;
  } else {
    for (int j = 0; j < p; j++) {
      pen += pen_default * D[j];
    }
    cout << "Default penalty: " << pen << endl;
    cout << "Default!" << endl;
  }
  
}


// [[Rcpp::export]]
void test_bool_param(NumericMatrix X, bool normalized = true) {
  
  // We create a new numeric vector of length n with a constructor: 
  // NumericVector out(n). 
  // Another useful way of making a vector is to copy an existing one: 
  // NumericVector zs = clone(ys).
  int n = X.nrow();
  NumericVector W(n);
  for (int i = 0; i < n; i++) {
    W[i] = 1;
  }
  if (normalized) {
    NumericMatrix Y = normalize_matrix(X, W);
    cout << "To normalize!" << endl;
  } else {
    NumericMatrix Y = clone(X);
    cout << "No normalize!" << endl;
  }
  if (normalized) {
    cout << "To normalize!" << endl;
  } else {
    cout << "No normalize!" << endl;
  }
  // cout << Y(0, 0) << endl;  // error
  // cannot be Y.nrow() here, still don't know why.
  int p = X.ncol();
  cout << "nrow of Y: " <<  n << endl;
  cout << "ncol of Y: " <<  p << endl;
}


// [[Rcpp::export]]
NumericVector test_vector_bounds(NumericVector H) {
  int n = H.size() - 1;
  NumericVector tau_hat(n + 1);
  for (int i = 0; i < n + 1; i++) {
    tau_hat[i] = -8;
  }
  int ncp = 0;
  int t_last = n + 1;  // H is [0, ..., n]
  int t_cur = H[t_last - 1];
  while (t_cur > 1) {
    tau_hat[ncp] = t_cur;
    ncp += 1;
    t_last = t_cur; 
    cout << t_cur << endl;
    t_cur = H[t_last - 1];
  }
  cout << "Changepoint number: " << ncp << endl;
  cout << "The last changepoint: " << t_last << endl;
  return tau_hat;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
void test_insert_set() {
  set<int> R = {6, 0, 2, 12, 34, -23, -100};  // the remaining set index;
  R.insert(2);  // insert some element
  R.insert(1);
  cout << "myset contains: " << endl;
  for (set<int>::iterator it = R.begin(); it != R.end(); it ++) {
    cout << ' ' << *it;
    if (*it < 0)
      R.erase(*it); // delete some element
  }
  cout << endl;
  
  cout << "erased myset contains: " << endl;
  for (set<int>::iterator it = R.begin(); it != R.end(); it ++) {
    cout << ' ' << *it;
  }
  cout << endl;
  
  set<int> R2 = R;
  R2.erase(1);  // delete some element
  cout << "after erase myset contains: " << endl;
  for (set<int>::iterator it = R2.begin(); it != R2.end(); it ++) {
    cout << ' ' << *it;
    
  }
  
  
}


// [[Rcpp::export]]
void test_init_vector_with_constant(int len = 4, int c = 6) {
  NumericVector vec(len); 
  for (int i = 0; i < len; i++) {
    vec[i] = c + i;
  }
  for (int i = 0; i < len; i++) {
    cout << ' ' << vec[i];
  }
  
}


