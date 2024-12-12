#include <Rcpp.h>
using namespace Rcpp;

// helper function for matrix vector multiplication
NumericVector matrix_vector_multiply(NumericMatrix mat, NumericVector vec) {
  int nrow = mat.nrow(), ncol = mat.ncol();
  NumericVector out(nrow);
  
  for (int i = 0; i < nrow; i++) {
    double sum = 0;
    for (int j = 0; j < ncol; j++) {
      sum += mat(i, j) * vec[j];
    }
    out[i] = sum;
  }
  return out;
}

// [[Rcpp::export]]
int sample_z_k_cpp(NumericVector beta, IntegerVector z, NumericMatrix X, NumericVector Y, int k, double prob) {
  /**
   * Sample z_k from p(z_k | z_{-k}, beta, X, Y)
   */
  NumericVector X_k = X(_, k-1); // Adjust for 0-indexing in C++
  NumericMatrix X_mk = X;
  // X_mk(_, k-1) = NumericVector(X.nrow(), 0); // Set column k to 0
  
  double beta_k = beta[k-1]; // Adjust for 0-indexing
  NumericVector beta_mk = clone(beta);
  beta_mk[k-1] = 0; // Set element k to 0
  
  IntegerVector z_mk = clone(z);
  z_mk[k-1] = 0; // Set element k to 0
  
  // Convert z_mk to a NumericVector for compatible operations
  NumericVector z_mk_numeric = as<NumericVector>(z_mk); 
  // compute masked beta_mk
  NumericVector masked_beta_mk = z_mk_numeric * beta_mk;
  NumericVector res_k = Y - matrix_vector_multiply(X_mk, masked_beta_mk);
  // NumericVector res_k = Y - (X_mk * (z_mk * beta_mk));
  
  double log_un_prob_1 = -0.5 * sum(pow(beta_k * X_k, 2)) + sum(res_k * beta_k * X_k) + log(prob);
  double log_un_prob_0 = log(1 - prob);
  
  double prob_1 = 1 / (1 + exp(log_un_prob_0 - log_un_prob_1));
  
  // Generate uniform random number
  double u = R::runif(0, 1);
  
  return (u < prob_1) ? 1 : 0;
}

// [[Rcpp::export]]
IntegerVector sample_z_cpp(NumericVector beta, IntegerVector z, NumericMatrix X, NumericVector Y, double prob) {
  /**
   * Sample z_1, ..., z_k iteratively (using updated z_{1, ..., i-1} when computing z_{i]}).
   */
  int p = z.size();
  IntegerVector z_new = clone(z); // Make a copy to avoid modifying the original z
  
  for(int k = 0; k < p; ++k) {
    z_new[k] = sample_z_k_cpp(beta, z_new, X, Y, k+1, prob); // k+1 due to 0-indexing in C++ and 1-indexing in R
  }
  
  return z_new;
}

// [[Rcpp::export]]
double lupost_beta_k_cpp(double beta_k, NumericVector beta_mk, IntegerVector z, NumericMatrix X, NumericVector Y, int k, double lambda) {
  /**
   * Compute unnormalised log posterior density of beta_k | beta_{-k}, z, X, Y
   */
  NumericVector X_k = X(_, k-1); // Adjust for 0-indexing
  NumericMatrix X_mk = X;
  // X_mk(_, k-1) = NumericVector(X.nrow(), 0); // Exclude column k
  // removing this line as does not need to be done (already masked by setting z_mk[k-1] = 0 below).
  
  IntegerVector z_mk = clone(z);
  z_mk[k-1] = 0; // Exclude element k
  
  // Compute residuals without the effect of the k-th predictor
  // Convert z_mk to a NumericVector for compatible operations
  NumericVector z_mk_numeric = as<NumericVector>(z_mk); 
  // compute masked beta_mk
  NumericVector masked_beta_mk = z_mk_numeric * beta_mk;
  NumericVector res_k = Y - matrix_vector_multiply(X_mk, masked_beta_mk);
  // NumericVector res_k = Y - (X_mk * (z_mk * beta_mk));
  
  double z_k = z[k-1]; // Current inclusion indicator for beta_k
  double penalty = lambda * std::abs(beta_k); // Laplace prior penalty
  
  // Log posterior computation
  return -0.5 * sum(pow(z_k * beta_k * X_k - res_k, 2)) - penalty;
}

// [[Rcpp::export]]
double sample_proposal_cpp(double beta_old, double b) {
  /**
   * Generate a sample from a Laplace distribution using the inverse transform method
   */
  double u = R::runif(0, 1) - 0.5; // Shift uniform distribution to center at 0
  double laplace_sample;
  if (u < 0) {
    laplace_sample = beta_old + b * log(1 + 2 * u);
  } else {
    laplace_sample = beta_old - b * log(1 - 2 * u);
  }
  return laplace_sample;
}

// [[Rcpp::export]]
double sample_beta_k_cpp(NumericVector beta, IntegerVector z, NumericMatrix X, NumericVector Y, int k, double lambda) {
  /**
   * Sample beta_k from p(beta_k | beta_{-k}, z, X, Y)
   */
  double beta_k_old = beta[k-1]; // Adjust for 0-indexing
  
  NumericVector beta_mk = clone(beta);
  beta_mk[k-1] = 0; // Set element k to 0
  
  double beta_k_prop = sample_proposal_cpp(beta_k_old, 1); // Assuming b=1 for simplicity
  
  double lupost_prop = lupost_beta_k_cpp(beta_k_prop, beta_mk, z, X, Y, k, lambda);
  double lupost_old = lupost_beta_k_cpp(beta_k_old, beta_mk, z, X, Y, k, lambda);
  
  double A = std::min(1.0, exp(lupost_prop - lupost_old));
  
  double u = R::runif(0, 1);
  return (u < A) ? beta_k_prop : beta_k_old;
}

// [[Rcpp::export]]
NumericVector sample_beta_cpp(NumericVector beta, IntegerVector z, NumericMatrix X, NumericVector Y, double lambda) {
  /**
   * Sample beta_1, ..., beta_p iteratively (using updated beta_{1, ..., i-1} when computing beta_{i]}).
   */
  int p = beta.size();
  NumericVector beta_new = clone(beta); // Make a copy to avoid modifying the original beta
  
  for(int k = 0; k < p; ++k) {
    beta_new[k] = sample_beta_k_cpp(beta_new, z, X, Y, k+1, lambda); // k+1 due to 0-indexing in C++ and 1-indexing in R
  }
  
  return beta_new;
}

// [[Rcpp::export]]
List sample_mc_cpp(int n_samples, NumericVector beta_init, IntegerVector z_init, NumericMatrix X, NumericVector Y, double lambda = 2.0, double prob = 0.05, bool verbose = false) {
  /**
   * Sample a markov chain whose limiting distribution is the posterior.
   */
  int p = beta_init.size();
  NumericMatrix beta_samples(n_samples, p);
  IntegerMatrix z_samples(n_samples, p);
  
  // Initialize the first row of samples
  for(int j = 0; j < p; ++j) {
    beta_samples(0, j) = beta_init[j];
    z_samples(0, j) = z_init[j];
  }
  
  for(int i = 1; i < n_samples; ++i) {
    // Convert current row to NumericVector and IntegerVector for compatibility
    NumericVector beta_current = beta_samples(i-1, _);
    IntegerVector z_current = z_samples(i-1, _);
    
    IntegerVector z_new = sample_z_cpp(beta_current, z_current, X, Y, prob);
    NumericVector beta_new = sample_beta_cpp(beta_current, z_new, X, Y, lambda);
    
    // Update the matrices with new samples
    for(int j = 0; j < p; ++j) {
      beta_samples(i, j) = beta_new[j];
      z_samples(i, j) = z_new[j];
    }
    
    if(verbose && (i % 100 == 0)) {
      Rcout << "Progress: " << (i / (double)n_samples) * 100 << "%. \n";
    }
  }
  
  return List::create(Named("z") = z_samples, Named("beta") = beta_samples);
}