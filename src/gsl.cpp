#include "gsl.h"

#define GSL_SUCCESS 0
#define GSL_EZERODIV 1
#define GSL_EBADLEN 2

namespace GSL {
/* plain gauss elimination, only not bothering with the zeroes
 *
 *       diag[0]  abovediag[0]             0   .....
 *  belowdiag[0]       diag[1]  abovediag[1]   .....
 *             0  belowdiag[1]       diag[2]
 *             0             0  belowdiag[2]   .....
 */
int solve_tridiag_nonsym(const std::vector<double>& diag, const std::vector<double>& abovediag,
                         const std::vector<double>& belowdiag, const std::vector<double>& rhs, std::vector<double>& x,
                         size_t N) {
  int status = GSL_SUCCESS;
  std::vector<double> alpha(N);
  std::vector<double> z(N);

  size_t i, j;

  /* Bidiagonalization (eliminating belowdiag)
   & rhs update
   diag' = alpha
   rhs' = z
   */
  alpha[0] = diag[0];
  z[0] = rhs[0];

  if (alpha[0] == 0) {
    status = GSL_EZERODIV;
  }

  for (i = 1; i < N; i++) {
    const double t = belowdiag[i - 1] / alpha[i - 1];
    alpha[i] = diag[i] - t * abovediag[i - 1];
    z[i] = rhs[i] - t * z[i - 1];
    if (alpha[i] == 0) {
      status = GSL_EZERODIV;
    }
  }

  /* backsubstitution */
  x[N - 1] = z[N - 1] / alpha[N - 1];
  if (N >= 2) {
    for (i = N - 2, j = 0; j <= N - 2; j++, i--) {
      x[i] = (z[i] - abovediag[i] * x[i + 1]) / alpha[i];
    }
  }

  if (status == GSL_EZERODIV) {
    std::cout << "Error : matrix must be positive definite!"
              << "\n";
  }

  return status;
}

int gsl_linalg_solve_tridiag(const std::vector<double>& diag, const std::vector<double>& abovediag,
                             const std::vector<double>& belowdiag, const std::vector<double>& rhs,
                             std::vector<double>& solution) {
  if (diag.size() != rhs.size()) {
    std::cout << "size of diag must match rhs" << std::endl;
    exit(GSL_EBADLEN);
  } else if (abovediag.size() != rhs.size() - 1) {
    std::cout << "size of abovediag must match rhs-1" << std::endl;
    exit(GSL_EBADLEN);
  } else if (belowdiag.size() != rhs.size() - 1) {
    std::cout << "size of belowdiag must match rhs-1" << std::endl;
    exit(GSL_EBADLEN);
  } else if (solution.size() != rhs.size()) {
    std::cout << "size of solution must match rhs" << std::endl;
    exit(GSL_EBADLEN);
  } else {
    return solve_tridiag_nonsym(diag, abovediag, belowdiag, rhs, solution, diag.size());
  }
  return 0;
}

}  // namespace GSL
