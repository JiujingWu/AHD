// [[Rcpp::depends(Rcpp, RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using RcppParallel::RMatrix;
using RcppParallel::parallelFor;
using RcppParallel::RVector;

/*
 Helper: sort the i0-th row of the *permuted* Dx (denoted Dx*)
 without materializing Dx*. We apply the permutation `perm` to both
 the row index and the column indices and return 1-based column
 order indices in `ord`.

 Dx: pairwise distance matrix for X (n x n)
 perm: length-n permutation in 0-based indices
 i0: 0-based row index
 ord: output; length n; values in 1..n (1-based columns after sorting)
 */
static inline void row_order_perm_one(const RMatrix<double>& Dx,
                                      const std::vector<int>& perm,
                                      int i0,
                                      std::vector<int>& ord)
{
  const int n = Dx.nrow();
  std::vector<std::pair<double,int>> v(n);
  const int ii = perm[i0];
  for (int j = 0; j < n; ++j) {
    const int jj = perm[j];
    // store (value, 1-based column index in the *reordered* layout)
    v[j] = std::make_pair(Dx(ii, jj), j+1);
  }
  std::stable_sort(v.begin(), v.end(),
                   [](const std::pair<double,int>& a, const std::pair<double,int>& b){
                     if (a.first != b.first) return a.first < b.first;
                     return a.second < b.second;
                   });
  for (int j = 0; j < n; ++j) ord[j] = v[j].second;
}

/*
 Rank with ties.method = "max" for one row of Dy after columns are
 reordered according to `ord1based` (which is the 1-based order
 returned from the Dx* row sort). We do *not* permute Dy itself.

 Dy: pairwise distance matrix for Y (n x n)
 ord1based: 1..n indices of the *reordered* columns for this row
 i0: 0-based row index in Dy
 rk: output ranks (1..n) with max-rank for ties, aligned to the
 *reordered* column positions (i.e., the j-th rank corresponds
 to ord1based[j]).
 */
static inline void rank_max_row_from_reordered_cols(const RMatrix<double>& Dy,
                                                    const std::vector<int>& ord1based,
                                                    int i0,
                                                    std::vector<int>& rk)
{
  const int n = Dy.nrow();
  std::vector<std::pair<double,int>> v(n);
  for (int j = 0; j < n; ++j) {
    int col = ord1based[j]-1;
    v[j] = std::make_pair(Dy(i0, col), j);
  }
  std::stable_sort(v.begin(), v.end(),
                   [](const std::pair<double,int>& a, const std::pair<double,int>& b){
                     if (a.first != b.first) return a.first < b.first;
                     return a.second < b.second;
                   });
  int pos = 0;
  while (pos < n) {
    int start = pos;
    double val = v[pos].first;
    while (pos+1 < n && v[pos+1].first == val) ++pos;
    int end = pos;
    int r = end + 1;  // max-rank for the tie block
    for (int t = start; t <= end; ++t) rk[v[t].second] = r;
    ++pos;
  }
}

/*
 A single-permutation AHD computation.
 Returns AHD and its four quadrant components.
 */
struct AHDObs {
  double AHD, AHD11, AHD12, AHD21, AHD22;
};

static inline AHDObs compute_AHD_one_perm(const RMatrix<double>& Dx,
                                          const RMatrix<double>& Dy,
                                          const std::vector<int>& perm)
{
  const int n = Dx.nrow();
  const double denom = (double)n - 2.0;

  std::vector<int> ord(n), rk(n);

  // accumulate four block contributions
  double S11=0.0, S12=0.0, S21=0.0, S22=0.0;

  std::vector<double> A1_col(n);
  for (int j = 0; j < n; ++j)
    A1_col[j] = ((j+1) > 2 ? ((double)(j+1) - 2.0)/denom : 0.0);

  for (int i = 0; i < n; ++i) {
    // 1) sort the i-th row of Dx* (via indices)
    row_order_perm_one(Dx, perm, i, ord);

    // 2) reorder i-th row of Dy according to ord and compute max-rank
    rank_max_row_from_reordered_cols(Dy, ord, i, rk);

    // 3) row-wise A1(i, ·) from ranks
    std::vector<double> A1_row(n);
    for (int j = 0; j < n; ++j) {
      int r = rk[j];
      A1_row[j] = (r > 2 ? ((double)r - 2.0)/denom : 0.0);
    }

    // 4) for each column j, compute A11/A22 counts, derive A12/A21,
    //    form Hellinger-based contributions, and accumulate
    for (int j = 0; j < n; ++j) {
      const double thr = A1_row[j];

      int c11 = 0;
      for (int k = 0; k <= j; ++k) if (A1_row[k] <= thr) ++c11;
      c11 -= 2; if (c11 < 0) c11 = 0;

      int c22 = 0;
      for (int k = j; k < n; ++k) if (A1_row[k] > thr) ++c22;

      const double A11 = c11 / denom;
      const double A22 = c22 / denom;
      const double A1_ = A1_col[j];       // column marginal
      const double A_1 = thr;             // row marginal
      const double A2_ = 1.0 - A1_;
      const double A_2 = 1.0 - A_1;

      const double A12 = std::max(0.0, A1_ - A11);
      const double A21 = std::max(0.0, A_1 - A11);

      auto sq = [](double z){ return z*z; };
      const double t11 = std::sqrt(std::max(0.0, A11)) - std::sqrt(std::max(0.0, A_1 * A1_));
      const double t12 = std::sqrt(std::max(0.0, A12)) - std::sqrt(std::max(0.0, A1_ * A_2));
      const double t21 = std::sqrt(std::max(0.0, A21)) - std::sqrt(std::max(0.0, A2_ * A_1));
      const double t22 = std::sqrt(std::max(0.0, A22)) - std::sqrt(std::max(0.0, A_2 * A2_));
      S11 += sq(t11); S12 += sq(t12); S21 += sq(t21); S22 += sq(t22);
    }
  }

  const double cst = 4.0 * ((double)n - 2.0);
  AHDObs out;
  out.AHD11 = cst * S11;
  out.AHD12 = cst * S12;
  out.AHD21 = cst * S21;
  out.AHD22 = cst * S22;
  out.AHD   = out.AHD11 + out.AHD12 + out.AHD21 + out.AHD22;
  return out;
}

/*
 Exported R function: AHD(Dx, Dy)
 Computes the observed AHD statistic and its four components without permutation.
 */
// [[Rcpp::export]]
Rcpp::List AHD(const Rcpp::NumericMatrix& Dx,
               const Rcpp::NumericMatrix& Dy) {
  const int n = Dx.nrow();
  if (n < 3) Rcpp::stop("n must be >= 3");
  if (Dx.ncol()!=n || Dy.nrow()!=n || Dy.ncol()!=n)
    Rcpp::stop("Dx and Dy must be n x n");

  // identity "permutation" 0...n-1
  std::vector<int> id(n);
  for (int i = 0; i < n; ++i) id[i] = i;

  AHDObs s = compute_AHD_one_perm(RcppParallel::RMatrix<double>(Dx),
                                  RcppParallel::RMatrix<double>(Dy),
                                  id);

  return Rcpp::List::create(
    Rcpp::_["AHD"]   = s.AHD,
    Rcpp::_["AHD11"] = s.AHD11,
    Rcpp::_["AHD12"] = s.AHD12,
    Rcpp::_["AHD21"] = s.AHD21,
    Rcpp::_["AHD22"] = s.AHD22
  );
}

/*
 Parallel worker: each thread processes a batch of permutations and
 writes a binary indicator (>= observed) into per-permutation slots.
 */
struct AHDPermWorker : public RcppParallel::Worker {
  const RMatrix<double> Dx;
  const RMatrix<double> Dy;
  const int n;
  const int B;
  const int* perms;  // flattened (B x n) permutation table (0-based)

  AHDObs obs;

  std::vector<long long> count_all, count_11, count_12, count_21, count_22;

  AHDPermWorker(const NumericMatrix& Dx_,
                const NumericMatrix& Dy_,
                const IntegerMatrix& perms_,
                const AHDObs& obs_)
    : Dx(Dx_), Dy(Dy_), n(Dx_.nrow()), B(perms_.nrow()), obs(obs_)
  {
    // copy permutations to contiguous memory (faster reads; avoid SEXP contention)
    perms = (int*)malloc(sizeof(int)*B*n);
    for (int b = 0; b < B; ++b)
      for (int i = 0; i < n; ++i)
        const_cast<int*>(perms)[b*n + i] = perms_(b,i) - 1; // store 0-based

    count_all.resize(B, 0);
    count_11.resize(B, 0);
    count_12.resize(B, 0);
    count_21.resize(B, 0);
    count_22.resize(B, 0);
  }

  ~AHDPermWorker() {
    if (perms) free((void*)perms);
  }

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> perm(n);
    for (std::size_t b = begin; b < end; ++b) {
      // read permutation b
      const int* row = perms + b*n;
      for (int i = 0; i < n; ++i) perm[i] = row[i];

      // compute statistics for this permutation
      AHDObs st = compute_AHD_one_perm(Dx, Dy, perm);

      // write per-permutation results (no sharing = no races)
      if (st.AHD   >= obs.AHD)   count_all[b] = 1;
      if (st.AHD11 >= obs.AHD11) count_11[b]  = 1;
      if (st.AHD12 >= obs.AHD12) count_12[b]  = 1;
      if (st.AHD21 >= obs.AHD21) count_21[b]  = 1;
      if (st.AHD22 >= obs.AHD22) count_22[b]  = 1;
    }
  }
};

/*
 Exported R function: AHD_test(Dx, Dy, perms)
 Permutation test with B permutations (perms is B x n, 1..n).
 Returns p-values for AHD and each quadrant component using (count+1)/(B+2).
 */
// [[Rcpp::export]]
List AHD_test(const NumericMatrix& Dx,
              const NumericMatrix& Dy,
              const IntegerMatrix& perms) {
  const int n = Dx.nrow();
  if (n < 3) stop("n must be >= 3");
  if (Dx.ncol()!=n || Dy.nrow()!=n || Dy.ncol()!=n)
    stop("Dx and Dy must be n x n");
  const int B = perms.nrow();
  if (perms.ncol() != n) stop("perms must be B x n (1..n)");

  // observed (no parallel)
  std::vector<int> id(n); for (int i=0;i<n;++i) id[i]=i;
  AHDObs obs = compute_AHD_one_perm(RMatrix<double>(Dx), RMatrix<double>(Dy), id);

  // parallel over permutations
  AHDPermWorker worker(Dx, Dy, perms, obs);
  parallelFor(0, B, worker);

  // reduce
  long long c_all=0, c11=0, c12=0, c21=0, c22=0;
  for (int b=0;b<B;++b) {
    c_all += worker.count_all[b];
    c11   += worker.count_11[b];
    c12   += worker.count_12[b];
    c21   += worker.count_21[b];
    c22   += worker.count_22[b];
  }

  // unbiased permutation p-values: (k+1)/(B+2)
  double pv   = (c_all + 1.0) / (B + 2.0);
  double pv11 = (c11   + 1.0) / (B + 2.0);
  double pv12 = (c12   + 1.0) / (B + 2.0);
  double pv21 = (c21   + 1.0) / (B + 2.0);
  double pv22 = (c22   + 1.0) / (B + 2.0);

  return List::create(
    _["pvalue"]   = pv,
    _["pvalue11"] = pv11,
    _["pvalue12"] = pv12,
    _["pvalue21"] = pv21,
    _["pvalue22"] = pv22
  );
}

