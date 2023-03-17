//
// Created by ted.tai on 22-11-18.
//

#pragma once
#include <utility>
#include "polynomial.hpp"


namespace msquare {
namespace pnc {
template <int N_DEG>
class Spline2D {
 public:
  enum class EvaluateStatus {
    OK = 0,
    OUT_OF_RANGE = 1,
  };

  using PII = std::pair<double, double>;
  using PolynomialType = Polynomial2D<N_DEG>;
  struct Circle {
    PII center;
    double radius;
    Circle(PII center, double r) : center(std::move(center)), radius(r) {}
  };
  Spline2D() = default;
  bool SetDomainsAndPolys(const std::vector<double>& domains,
                          const std::vector<PolynomialType>& polys) {
    if (polys.empty()) {
      return false;
    }
    if (domains.size() != polys.size() + 1) {
      return false;
    }
    if (!std::is_sorted(domains.begin(), domains.end())) {
      return false;
    }
    accumulate_domains_ = domains;
    polys_ = polys;
    circles_.clear();
    circles_.reserve(polys_.size());
    for (int i = 0; i < polys_.size(); i++) {
      double l = accumulate_domains_[i + 1] - accumulate_domains_[i];
      const auto& seg = polys_[i];
      const PII pt_center = seg.Evaluate(0.5 * l);
      const auto radius = l / 2 + 5.0;
      circles_.emplace_back(pt_center, radius);
    }
    return true;
  }
  int num_segments() const { return static_cast<int>(polys_.size()); }
  const std::vector<double>& accumulate_domains() const {
    return accumulate_domains_;
  }
  double start_domain() const { return accumulate_domains_.front(); }
  double end_domain() const { return accumulate_domains_.back(); }
  const std::vector<PolynomialType>& polys() const { return polys_; }
  std::pair<EvaluateStatus, PII> Evaluate(const double s,
                                          const int degree) const {
    if (s < accumulate_domains_.front()) {
      return {EvaluateStatus::OUT_OF_RANGE, polys_.front().Evaluate(s, degree)};
    }
    if (s >= accumulate_domains_.back()) {
      return {
          EvaluateStatus::OUT_OF_RANGE,
          polys_.back().Evaluate(
              s - accumulate_domains_[accumulate_domains_.size() - 2], degree)};
    }
    auto idx = fplus::find_last_idx_by(
        [s](const double& domain) { return s >= domain; }, accumulate_domains_);
    return {EvaluateStatus::OK,
            polys_[idx.unsafe_get_just()].Evaluate(
                s - accumulate_domains_[idx.unsafe_get_just()], degree)};
  }
  std::pair<EvaluateStatus, PII> Evaluate(const double s) const {
    if (s < accumulate_domains_.front()) {
      return {EvaluateStatus::OUT_OF_RANGE, polys_.front().Evaluate(0.0)};
    }
    if (s >= accumulate_domains_.back()) {
      return {EvaluateStatus::OUT_OF_RANGE,
              polys_.back().Evaluate(
                  s - accumulate_domains_[accumulate_domains_.size() - 2])};
    }
    auto idx = fplus::find_last_idx_by(
        [s](const double& domain) { return s >= domain; }, accumulate_domains_);
    return {EvaluateStatus::OK,
            polys_[idx.unsafe_get_just()].Evaluate(
                s - accumulate_domains_[idx.unsafe_get_just()])};
  }
  Polynomial<N_DEG>& operator()(const uint8_t idx, uint8_t dim) {
    assert(idx < polys_.size());
    return polys_[idx](dim);
  }
  const PolynomialType& PolyAt(const uint8_t idx) const {
    assert(idx < polys_.size());
    return polys_[idx];
  }
  const Circle& CircleAt(const uint8_t idx) const {
    assert(idx < polys_.size());
    return circles_[idx];
  }
  std::string DebugString() const {
    std::stringstream ss;
    ss << "Spline2D_" << N_DEG;
    ss << "/n accumulate_domains: ";
    ss << fplus::show_cont(accumulate_domains_);
    ss << "/n polys: ";
    for (const auto& poly : polys_) {
      ss << poly.DebugString();
    }
    return ss.str();
  }

 private:
  std::vector<PolynomialType> polys_;
  std::vector<double> accumulate_domains_;
  std::vector<Circle> circles_;
};
}  // namespace pnc
}  // namespace npp