#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

#include "timer.h"

template <typename T>
std::vector<T> generateRandomVector(size_t size, T minValue, T maxValue) {
  std::random_device rnd_device;
  std::mt19937 mersenne_engine{rnd_device()};
  std::uniform_real_distribution<T> dist{minValue, maxValue};
  auto gen = [&dist, &mersenne_engine]() { return dist(mersenne_engine); };
  std::vector<T> data(size);
  std::generate(data.begin(), data.end(), gen);
  return data;
}

struct indGenerator {
  int current_;
  indGenerator(int start) : current_(start) {}
  int operator()() { return current_++; }
};

template <typename T>
std::vector<T> generateIncrementalVector(size_t size) {
  indGenerator g(1);
  std::vector<T> data(size);
  std::generate(data.begin(), data.end(), g);
  return data;
}

void saveGslVector(std::string filename, gsl_vector* x, int size) {
  std::ofstream outfile(filename.c_str());
  if (outfile.is_open()) {
    outfile << std::scientific << std::setprecision(4);
    for (size_t i = 0; i < size; ++i) {
      outfile << gsl_vector_get(x, i) << "\n";
    }
    outfile.close();
  }
}

class TestMatrix {
 public:
  TestMatrix(int _N) : N(_N) { initVectors(); }

  void initVectors() {
    lower_diag = generateRandomVector<double>(N - 1, 0., 1.);
    upper_diag = generateRandomVector<double>(N - 1, 0., 1.);
    central_diag = generateRandomVector<double>(N, 1., 10.);
    rhs = generateIncrementalVector<double>(N);
  }

  void computeLU() {
    gsl_vector* b = gsl_vector_alloc(N);
    for (size_t i = 0; i < N; ++i) {
      gsl_vector_set(b, i, rhs.at(i));
    }

    gsl_matrix* m = gsl_matrix_alloc(N, N);
    for (size_t i = 0; i < N; ++i)
      for (size_t j = 0; j < N; ++j) {
        if (i == j)
          gsl_matrix_set(m, i, j, central_diag[i]);
        else if (i == j - 1)
          gsl_matrix_set(m, i, j, upper_diag[i]);
        else if (i == j + 1)
          gsl_matrix_set(m, i, j, lower_diag[j]);
        else
          gsl_matrix_set(m, i, j, 0.);
      }

    gsl_vector* x = gsl_vector_alloc(N);
    int s;
    gsl_permutation* p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(m, p, &s);

    {
      NM::Timer timer;
      for (int i = 0; i < 10000; ++i) gsl_linalg_LU_solve(m, p, b, x);
    }

    saveGslVector("LU.txt", x, N);

    gsl_permutation_free(p);
    gsl_vector_free(x);

    gsl_vector_free(b);
    gsl_matrix_free(m);
  }

  void tridiag() {
    gsl_vector* x = gsl_vector_alloc(N);

    auto gsl_central_diag = gsl_vector_const_view_array(&central_diag[0], N);
    auto gsl_lower_diag = gsl_vector_const_view_array(&lower_diag[0], N - 1);
    auto gsl_upper_diag = gsl_vector_const_view_array(&upper_diag[0], N - 1);
    auto gsl_rhs = gsl_vector_const_view_array(&rhs[0], N);

    {
      NM::Timer timer;
      for (int i = 0; i < 10000; ++i)
        gsl_linalg_solve_tridiag(&gsl_central_diag.vector, &gsl_upper_diag.vector, &gsl_lower_diag.vector,
                                 &gsl_rhs.vector, x);
    }

    saveGslVector("tridiag.txt", x, N);

    gsl_vector_free(x);
  }

 protected:
  int N = 0;
  std::vector<double> central_diag;
  std::vector<double> lower_diag;
  std::vector<double> upper_diag;
  std::vector<double> rhs;
};

int main(void) {
  TestMatrix matrix(300);

  matrix.computeLU();
  matrix.tridiag();

  return 0;
}
