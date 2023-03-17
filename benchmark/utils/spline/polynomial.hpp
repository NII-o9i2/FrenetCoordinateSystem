#include "array"
#include "string"
#pragma once

#include "fplus/fplus.hpp"

namespace msquare {
namespace pnc {
constexpr int N_DEG = 5;
constexpr int N_DIM = 2;
constexpr int N_CONTINUOUS = 4;  // continuity up to jerk should be enough
constexpr int REGULATOR_TOP_DEGREE = 3;
constexpr std::array<int, 9> FACTORIALS = {1, 1, 2, 6, 24, 120, 720, 5040, 40320};

template <int N_DEG>
class Polynomial {
  /**
   * @brief Polynomial class, coeffs are stored in reversed order
   * The parameterization is given by
   * f(s) = coeffs_(n) + coeffs_(n-1)/1!^s + ... + coeffs_(0)/(n!)*s^n
   * f(s) = coeffs_normal_order_(0) + coeffs_normal_orders_(1)^s + ... +
   * coeffs_normal_order_(n)*s^n
   * coeffs are scaled to have straightforward physical
   * meaning, also improve the efficiency when evaluating the derivatives.
   */
 public:
  using Coeffs = std::array<double, N_DEG + 1>;
  Polynomial() {
    coeffs_.fill(0.0);
    coeffs_normal_order_.fill(0.0);
  }
  explicit Polynomial(const Coeffs& coeffs) : coeffs_(coeffs) {
    for (int i = 0; i < coeffs_.size(); i++) {
      coeffs_normal_order_[i] = coeffs_[N_DEG - i] / FACTORIALS[i];
    }
  }
  inline double Evaluate(const double s, const int degree) const {
    double ret = coeffs_[0] / FACTORIALS[N_DEG - degree];
    for (int i = 1; i + degree <= N_DEG; i++) {
      ret = ret * s + coeffs_[i] / FACTORIALS[N_DEG - degree - i];
    }
    return ret;
  }
  inline double Evaluate(const double s) const {
    double ret = coeffs_normal_order_[N_DEG];
    for (int i = 1; i <= N_DEG; i++) {
      ret = ret * s + coeffs_normal_order_[N_DEG - i];
    }
    return ret;
  }
  std::string DebugString() const {
    std::stringstream ss;
    ss << "Polynomial_" << N_DEG;
    ss << "/n coeffs: ";
    ss << fplus::show_cont(coeffs_);
    ss << "/n coeffs_normal_order: ";
    ss << fplus::show_cont(coeffs_normal_order_);
    return ss.str();
  }

 private:
  Coeffs coeffs_;
  Coeffs coeffs_normal_order_;
};

template <int N_DEG>
class Polynomial2D {
 public:
  using Coeffs = std::array<double, N_DEG + 1>;
  Polynomial2D() = default;
  Polynomial2D(const Coeffs& coeffs_x, const Coeffs& coeffs_y)
      : polys_({Polynomial<N_DEG>(coeffs_x), Polynomial<N_DEG>(coeffs_y)}) {}
  inline std::pair<double, double> Evaluate(const double s,
                                            const int degree) const {
    return std::make_pair(polys_[0].Evaluate(s, degree),
                          polys_[1].Evaluate(s, degree));
  }
  inline std::pair<double, double> Evaluate(const double s) const {
    return {polys_[0].Evaluate(s), polys_[1].Evaluate(s)};
  }
  std::string DebugString() const {
    std::stringstream ss;
    ss << "Polynomial2D_" << N_DEG;
    ss << "/n x: ";
    ss << polys_[0].DebugString();
    ss << "/n y: ";
    ss << polys_[1].DebugString();
    return ss.str();
  }
  Polynomial<N_DEG>& operator[](const unsigned char j) {
    assert(j < 2);
    return polys_[j];
  }

 private:
  std::array<Polynomial<N_DEG>, 2> polys_;
};
using QuinticPoly = Polynomial<N_DEG>;
using QuinticPoly2D = Polynomial2D<N_DEG>;

}  // namespace pnc
}  // namespace msquare
