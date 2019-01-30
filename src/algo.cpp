// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
using namespace arma;



// 00) Prelim
void nan_to_minus_infinity(double & x) {
  // turn NaN x to -Inf
  if (isnan(x)) {
    x = -INFINITY;
  }
}

const bool update(double & par_curr, const double par_prop, double & lpost_curr, const double lpost_prop, double & s, const int i, const int burnin, const double factor = 30.0) {
  // M-H update
  const bool accept_reject = log(runif(1)[0]) < lpost_prop - lpost_curr;
  par_curr = accept_reject ? par_prop : par_curr;
  lpost_curr = accept_reject ? lpost_prop : lpost_curr;
  if (i < burnin) {
    s = sqrt(s * s + (accept_reject ? 3.0 : (-1.0)) * s * s / factor / sqrt(i + 1.0));
  }
  return accept_reject;
}

DataFrame df_scalars(const double alpha_init, const double s_alpha_init, const double s_alpha_final, const double a, const double b, const int n_swap, const int N, const int thin, const int burnin, const int print_freq) {
  DataFrame scalars =
    DataFrame::create(Named("alpha_init") = NumericVector::create(alpha_init),
                      Named("s_alpha_init") = NumericVector::create(s_alpha_init),
                      Named("s_alpha_final") = NumericVector::create(s_alpha_final),
                      Named("a") = NumericVector::create(a),
                      Named("b") = NumericVector::create(b),
                      Named("n_swap") = IntegerVector::create(n_swap),
                      Named("N") = IntegerVector::create(N),
                      Named("thin") = IntegerVector::create(thin),
                      Named("burnin") = IntegerVector::create(burnin),
                      Named("print_freq") = IntegerVector::create(print_freq));
  return scalars;
}

const arma::mat table(const arma::imat X, const int K) {
  // matrix of counts: element (p, i) is count of i in X's pth row
  if (K <= 0) {
    stop("table: K has to be positive.");
  }
  const int n = X.n_rows;
  arma::mat counts(n, K);
  arma::uvec indices;
  for (int p = 0; p < n; p++) { 
    for (int i = 0; i < K; i++) {
      indices = find(X.row(p) == i);
      counts(p, i) = (double) indices.size();
    }
  }
  return counts;
}

template <class T>
T reorder(T Y, const arma::uvec o) {
  // reorder rows & columns of Y simultaneously
  return Y.submat(o, o);
}

//' Reorder rows and columns of a square matrix
//'
//' @param Y a square matrix
//' @param o a permutation vector of 0,1,...,(n-1), where n is the number of rows/columns of Y
//' @return the reordered square matrix of same size as Y
//' @export
// [[Rcpp::export]]
const arma::mat reorder_once(const arma::mat Y, const arma::uvec o) {
  return reorder(Y, o);
}

const arma::mat reorder_twice(const arma::mat Y, const arma::uvec o) {
  // expect to return Y
  return reorder(reorder(Y, o), sort_index(o));
}

const arma::mat checks(const arma::mat Y,
                       const arma::uvec o_init,
                       const double alpha_init,
                       const double s_alpha_init,
                       const double a,
                       const double b,
                       const arma::mat A,
                       const arma::mat B
                       ) {
  // generic checks at the beginning of an MCMC algorithm
  // initialises Y_star = reorder(Y, o_init) if everything passes
  const int n = Y.n_cols, K = A.n_cols;
  if (!Y.is_square()) {
    stop("Y has to be a square matrix.");
  }
  if (o_init.size() != n) {
    stop("Length of o_init has to be n.");
  }
  const IntegerVector seq_n = seq_len(n) - 1,
    o_sorted = wrap(sort(o_init));
  if (is_true(any(seq_n != o_sorted))) {
    stop("o_init has to be a permutation of {0,1,...,n-1}.");
  }
  const arma::mat Y_star = reorder(Y, o_init),
    Y_lower = trimatl(Y_star); // including main diagonal
  if (any(vectorise(Y_lower))) {
    stop("o_init has to be such that reorder(Y, o_init) is upper triangular.");
  }
  const vec vals_alpha = {alpha_init, s_alpha_init, a, b};
  if (any(vals_alpha <= 0.0)) {
    stop("Initial value, proposed standard deviation & hyperparameters for alpha have to be positive.");
  }
  const arma::uvec dims = {A.n_rows, B.n_cols, B.n_rows};
  if (any(dims != K)) {
    stop("A & B have to be square matrices of same dimensions.");
  }
  if (any(vectorise(A) <= 0.0) || any(vectorise(B) <= 0.0)) {
    stop("All elements of A & B have to be positive.");
  }
  return Y_star;
}

const arma::uvec indices_ij(const int i, const int j, const arma::imat Z_upper, const arma::imat Z_lower) {
  // find overlapping indices
  // Z_upper is upper tri., Z_lower is lower tri.
  const arma::uvec ui = find(Z_upper == i),
    uj = find(Z_lower.t() == j),
    u = intersect(ui, uj);
  return u;
}



// 01) simulate functions - random
const int sample_1(const IntegerVector seq, const NumericVector l = NumericVector::create()) {
  // sample 1 value from seq w/ weights exp(l) (if non-empty)
  IntegerVector seq0 = clone(seq);
  NumericVector p = clone(l);
  if (l.size() != 0) {
    NumericVector l0 = clone(l);
    while (max(l0) - min(l0) > 700.0) { // to prevent overflow (due to underflow)
      seq0 = seq0[l0 != min(l0)];
      l0 = l0[l0 != min(l0)];
    }
    p = exp(l0 - min(l0));
  }
  return Rcpp::RcppArmadillo::sample(seq0, 1, true, p)[0];
}

const arma::mat sim_C(const int K, const arma::mat A, const arma::mat B, const arma::imat Z_upper, const arma::imat Z_lower, const arma::mat Y_star) {
  // simulate C (for both prior & cond. posterior)
  arma::mat C(K, K);
  arma::uvec u;
  double a0, b0;
  const arma::mat X_star = trimatu(1.0 - Y_star, 1); // omit main diag
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      u = indices_ij(i, j, Z_upper, Z_lower);
      a0 = A(i, j);
      b0 = B(i, j);
      if (u.size() != 0)  {
        a0 += sum(Y_star(u));
        b0 += sum(X_star(u));
      }
      C(i, j) = rbeta(1, a0, b0)[0];
    }
  }
  return C;
}

const arma::mat sim_C_undir(const int K, const arma::mat A, const arma::mat B, const arma::imat Z_upper, const arma::imat Z_lower, const arma::mat Y) {
  // simulate C (for both prior & cond. posterior)
  arma::mat C(K, K);
  arma::uvec u;
  double a0, b0;
  const arma::mat X = trimatu(1.0 - Y, 1);
  for (int i = 0; i < K; i++) {
    for (int j = i; j < K; j++) {
      u = indices_ij(i, j, Z_upper, Z_lower);
      a0 = A(i, j);
      b0 = B(i, j);
      if (u.size() != 0) {
        a0 += sum(Y(u));
        b0 += sum(X(u));
      }
      if (j != i) {
        u = indices_ij(j, i, Z_upper, Z_lower);
        if (u.size() != 0) {
          a0 += sum(Y(u));
          b0 += sum(X(u));
        }
      }
      C(i, j) = rbeta(1, a0, b0)[0];
      C(j, i) = C(i, j);
    }
  }
  return C;
}

const arma::mat sim_D(const int n, const int K, const arma::mat alphas) {
  // simulate D (for both prior & cond. posterior)
  // alphas assumed to be nxK matrix
  NumericVector d(K);
  arma::mat D(n, K);
  for (int p = 0; p < n; p++) {
    for (int i = 0; i < K; i++) {
      d[i] = rgamma(1, alphas(p, i), 1.0)[0];
    }
    d = d / sum(d); // normalising gives Dirichlet
    D.row(p) = as<rowvec>(d); // each row sums to 1.0
  }
  return D;
}

const arma::imat sim_Z(const IntegerVector seq_K, const arma::mat D) {
  // simulate Z (for initialisation)
  // D assumed to be nxK matrix
  const int n = D.n_rows;
  NumericVector d;
  IntegerVector z(n);
  arma::imat Z(n, n);
  for (int p = 0; p < n; p++) {
    d = wrap(D.row(p));
    z = Rcpp::RcppArmadillo::sample(seq_K, n, true, d);
    Z.row(p) = as<irowvec>(z);
  }
  Z.diag().fill(-1); // Zpp = -1 & never updated
  return Z;
}

const irowvec sim_Z_partial(const IntegerVector seq_K, const arma::mat L) {
  // simulate Z (Gibbs step)
  // L assumed to be mxK matrix of log probabilities
  const int m = L.n_rows, K = seq_K.size();
  rowvec l_arma(K);
  NumericVector l_rcpp(K);
  irowvec z(m);
  for (int q = 0; q < m; q++) {
    l_arma = L.row(q);
    l_rcpp = wrap(l_arma); // will fail if wrap(L.row(q)) - why?
    z[q] = sample_1(seq_K, l_rcpp);
  }
  return z;
}



// 02) Regular Gibbs sampler (RGS) of MMSBM for DAG
const double lpost_rgs_alpha(const double alpha, const int n, const int K, const arma::mat D, const double a, const double b) {
  // log conditional posterior of alpha in RGS
  double lpost;
  if (alpha <= 0.0) {
    lpost = -INFINITY;
  }
  else {
    const double n0 = (double) n, K0 = (double) K;
    const NumericVector alpha0 = NumericVector::create(alpha);
    lpost = (alpha - 1.0) * accu(log(D)) + n0 * lgamma(K0 * alpha) - n0 * K0 * lgamma(alpha) + dgamma(alpha0, a, 1.0 / b, true)[0];
  }
  nan_to_minus_infinity(lpost);
  return lpost;
}

//' Run the regular Gibbs sampler of MMSBM
//'
//' @param Y the adjacency matrix representing the graph
//' @param o_init initial value of the ordering vector
//' @param alpha_init initial value of alpha
//' @param s_alpha_init initial value of proposal standard deviation for alpha's Metropolis step
//' @param a,b hyperparameters for alpha's prior
//' @param A,B square matrices of the hyperparameters for group-to-group probabilities' prior
//' @param n_swap number of swaps in order per iteration
//' @param N,thin,burnin,print_freq MCMC quantities
//' @return a list of data frames and matrices representing the MCMC output
//' @export
// [[Rcpp::export]]
List rgs_mmsbm(const arma::mat Y,
               const arma::uvec o_init,
               const double alpha_init,
               const double s_alpha_init,
               const double a,
               const double b,
               const arma::mat A,
               const arma::mat B,
               const int n_swap = 1,
               const int N = 1000,
               const int thin = 1,
               const int burnin = 100,
               const int print_freq = 100
               ) {
  // a) initialisation: data & latent variables
  const int n = Y.n_cols, K = A.n_cols;
  const IntegerVector seq_n = seq_len(n) - 1,
    seq_K = seq_len(K) - 1;
  arma::mat Y_star = checks(Y, o_init, alpha_init, s_alpha_init, a, b, A, B),
    alphas(n, K);
  alphas.fill(alpha_init);
  arma::imat Z_upper(n, n), Z_lower(n, n);
  Z_upper.fill(-1);
  Z_lower.fill(-1);
  arma::mat C = sim_C(K, A, B, Z_upper, Z_lower, Y_star),
    D = sim_D(n, K, alphas), D_star(n, K);
  arma::imat Z = sim_Z(seq_K, D), Z_star(n, n);
  // b) initialisation: for updating
  arma::mat counts(n, K); // for D
  double alpha_curr = alpha_init, alpha_prop, // for alpha
    s_alpha = s_alpha_init, // for alpha
    lpost_curr, lpost_prop; // for alpha
  bool aor; // for alpha
  ivec Z_temp; // for Z
  arma::mat C_temp, Y_temp, D_temp, L_temp; // for Z
  arma::uvec o = o_init; // for o
  IntegerVector seq_temp; // for o
  double accept_prob; // for o
  int o_temp; // for o
  // c) initialisation: for saving
  vec alpha_vec(N);
  arma::mat C_mat(N, K * K), D_mat(N, n * K);
  arma::umat o_mat(N, n), i_mat(N, n);
  // d) run
  int i, j, k, p, q, t;
  for (t = 0; t < N * thin + burnin; t++) {
    // update all star matrices 1st
    Z_star = reorder(Z, o);
    Y_star = reorder(Y, o);
    D_star = D.rows(o);
    // update C & D: conjugate, Gibbs
    Z_upper = trimatu(1 + Z_star) - 1; // set lower tri. to -1
    Z_lower = trimatl(1 + Z_star) - 1; // set upper tri. to -1
    C = sim_C(K, A, B, Z_upper, Z_lower, Y_star);
    counts = table(Z, K);
    D = sim_D(n, K, alpha_curr + counts);
    // update alpha: Metropolis
    alpha_prop = rnorm(1, alpha_curr, s_alpha)[0];
    lpost_prop = lpost_rgs_alpha(alpha_prop, n, K, D, a, b);
    lpost_curr = lpost_rgs_alpha(alpha_curr, n, K, D, a, b);
    aor = update(alpha_curr, alpha_prop, lpost_curr, lpost_prop, s_alpha, t, burnin);
    // update Z_star & Z
    for (p = 0; p < n-1; p++) {
      Z_temp = Z_star(span(p+1, n-1), p); // length (n-1-p)
      // all (n-1-p)xK matrices below
      C_temp = C.cols(conv_to<uvec>::from(Z_temp)).t();
      Y_temp = repmat(Y_star(p, span(p+1, n-1)).t(), 1, K);
      D_temp = repmat(D_star.row(p), (n-1-p), 1);
      L_temp = Y_temp % log(C_temp) + (1.0 - Y_temp) % log(1.0 - C_temp) + log(D_temp);
      Z_star(p, span(p+1, n-1)) = sim_Z_partial(seq_K, L_temp);
    }
    for (q = 1; q < n; q++) {
      Z_temp = Z_star(span(0, q-1), q); // length q
      // all qxK matrices below
      C_temp = C.rows(conv_to<uvec>::from(Z_temp));
      Y_temp = repmat(Y_star(span(0, q-1), q), 1, K);
      D_temp = repmat(D_star.row(q), q, 1);
      L_temp = Y_temp % log(C_temp) + (1.0 - Y_temp) % log(1.0 - C_temp) + log(D_temp);
      Z_star(q, span(0, q-1)) = sim_Z_partial(seq_K, L_temp);
    }
    Z = reorder(Z_star, sort_index(o));
    // update o
    for (k = 0; k < n_swap; k++) {
      p = sample_1(seq_n);
      if (p == 0) {
        q = 1;
      }
      else if (p == n-1) {
        q = n-2;
      }
      else {
        q = (runif(1)[0] < 0.5) ? p-1 : p+1;
      }
      if (Y_star(min(p,q), max(p,q)) != 1) { // st. reject if equals 1
        i = Z_star(p, q);
        j = Z_star(q, p);
        accept_prob = (1.0 - C(j, i)) / (1.0 - C(i, j));
        if (p == 0 || p == n-1) {
          accept_prob *= 0.5;
        }
        else if (q == 0 ||  q == n-1) {
          accept_prob *= 2.0;
        }
        if (runif(1)[0] < accept_prob) {
          o_temp = o[p];
          o[p] = o[q];
          o[q] = o_temp;
          Z_star.swap_rows(p, q);
          Z_star.swap_cols(p, q);
          Y_star.swap_rows(p, q);
          Y_star.swap_cols(p, q);
          D_star.swap_rows(p, q);
        }
      }
    }
    // print & save
    if ((t+1) % print_freq == 0) {
      Rcout << "Iteration " << t+1 << endl;
      Rcout << "alpha (curr): " << alpha_curr << endl;
      if (t < burnin) {
        Rcout << "s_alpha: " << s_alpha << endl;
      }
      Rcout << "C (curr): " << C << endl;
    }
    if (t >= burnin && (t - burnin + 1) % thin == 0) {
      k = (t - burnin + 1) / thin - 1;
      alpha_vec[k] = alpha_curr;
      C_mat.row(k) = vectorise(C, 1);
      D_mat.row(k) = vectorise(D, 1);
      o_mat.row(k) = conv_to<urowvec>::from(o);
      i_mat.row(k) = conv_to<urowvec>::from(sort_index(o));
    }
  }
  // e) save
  DataFrame scalars = df_scalars(alpha_init, s_alpha_init, s_alpha, a, b, n_swap, N, thin, burnin, print_freq);
  List L = List::create(Named("alpha") = alpha_vec,
                        Named("C") = C_mat,
                        Named("D") = D_mat,
                        Named("o") = o_mat,
                        Named("i") = i_mat,
                        Named("scalars") = scalars,
                        Named("Y") = Y,
                        Named("o_init") = o_init,
                        Named("A") = A,
                        Named("B") = B);
  return L;
}



// 03) all other possible samplers
// List rgs_mmsbm_dag() {} // same as rgs_mmsbm()

// List rgs_mmsbm_und() {}

// List mwg_mmsbm_dag() {}

// List mwg_mmsbm_und() {}

// List rgs_a_sbm_dag() {}

// List rgs_a_sbm_und() {}

// List mwg_a_sbm_dag() {}

// List mwg_a_sbm_und() {}

// List rgs_amsbm_dag() {}

// List rgs_amsbm_und() {}

// List mwg_amsbm_dag() {}

// List mwg_amsbm_und() {}
