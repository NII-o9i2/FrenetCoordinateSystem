//
// Created by jojo-feng on 23-2-1.
//

#ifndef BENCHMARK_FRENET_COORDINATE_SYSTEM_FRENET_COORDINATE_SYSTEM_NEW_H_
#define BENCHMARK_FRENET_COORDINATE_SYSTEM_FRENET_COORDINATE_SYSTEM_NEW_H_

#include "lane.hpp"
#include "frenet_coordinate_system.h"
#include "spline/spline_generator.hpp"
#include <cassert>
#include <vector>

namespace msquare {
class FrenetCoordinateSystemNew {
 public:
  FrenetCoordinateSystemNew() = default;
  ~FrenetCoordinateSystemNew() = default;
  explicit FrenetCoordinateSystemNew(
      const std::vector<Point2D> &vec_pts,
      const FrenetCoordinateSystemParameters &fcs_params);

  bool Available() const { return available_; }
  double GetLength() const { return lane_.length(); }
  double GetSlength() const { return lane_.length(); }
  Point2D GetRefCurvePoint(double s) const {
    assert(s <= GetLength() && s >= 0);
    return lane_.GetPositionUnsafe(s);
  }
  double GetRefCurveHeading(double s) const {
    assert(s <= GetLength() && s >= 0);
    return lane_.GetHeadingUnsafe(s);
  }
  double GetRefCurveCurvature(double s) const {
    assert(s <= GetLength() && s >= 0);
    return lane_.GetCurvatureUnsafe(s);
  }
  double GetRefCurveDCurvature(double s) const {
    assert(s <= GetLength() && s >= 0);
    return lane_.GetDCurvatureUnsafe(s);
  }

  TRANSFORM_STATUS CartCoord2FrenetCoord(const Point2D &cart, Point2D &frenet,
                                         bool has_heuristics = false,
                                         double s_begin = 0.0,
                                         double s_end = 120.0) const;
  TRANSFORM_STATUS FrenetCoord2CartCoord(const Point2D &frenet,
                                         Point2D &cart) const;
  TRANSFORM_STATUS CartState2FrenetState(const CartesianState &cart_state,
                                         FrenetState &frenet_state) const;
  TRANSFORM_STATUS FrenetState2CartState(const FrenetState &frenet_state,
                                         CartesianState &cart_state) const;

 private:
  FrenetCoordinateSystemParameters fcs_params_{}; // this value is no use
  bool available_{false};
  pnc::Lane lane_;
  pnc::SplineFitting::SplineFittingParams fitting_params_;
};
}  // namespace msquare

#endif //BENCHMARK_FRENET_COORDINATE_SYSTEM_FRENET_COORDINATE_SYSTEM_NEW_H_
