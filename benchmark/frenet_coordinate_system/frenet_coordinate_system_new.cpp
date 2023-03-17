//
// Created by jojo-feng on 23-2-1.
//

#include "frenet_coordinate_system_new.h"
#include "lane.hpp"
#include "cartesian_frenet_conversion.h"

namespace msquare {
using planning_math::CartesianFrenetConverter;
using pnc::SplineFitting;
using std::array;
using std::vector;
FrenetCoordinateSystemNew::FrenetCoordinateSystemNew(
    const vector<Point2D> &vec_pts,
    const FrenetCoordinateSystemParameters &fcs_params)
    : fcs_params_(fcs_params) {
  fitting_params_.has_coarse_projection_method = true;
  const auto accumulate_s = SplineFitting::GenAccumulateS(vec_pts);
  const auto fitting_res_with_time =
      SplineFitting::GenPiecewiseQuinticSplineWithTimeCost(
          vec_pts, accumulate_s, fitting_params_, &lane_);
  const auto &time_cost_ns = fitting_res_with_time.first;
//  NLOGE("[frenet construct time cost]: %lu(ns)", time_cost_ns);
  const bool is_fitting_success = fitting_res_with_time.second.first;
  if (!is_fitting_success) {
    const auto &fail_msg = fitting_res_with_time.second.second;
//    NLOGE("[frenet construct failed]: %s", fail_msg.c_str());
    return;
  }
  available_ = true;
}

TRANSFORM_STATUS FrenetCoordinateSystemNew::CartCoord2FrenetCoord(
    const Point2D &cart, Point2D &frenet, bool has_heuristics, double s_begin,
    double s_end) const {
  if (!Available()) {
    return TRANSFORM_FAILED;
  }
  auto to_point_2d = [](const pnc::Lane::ProjectionPt &pt) -> Point2D {
    return {pt.t, pt.d};
  };
  if (has_heuristics) {
    frenet = to_point_2d(lane_.GetProjectionAtDomain(cart, s_begin, s_end));
  } else {
    frenet = to_point_2d(lane_.GetProjectionCoarse(cart));
  }
  return TRANSFORM_SUCCESS;
}

TRANSFORM_STATUS
FrenetCoordinateSystemNew::FrenetCoord2CartCoord(const Point2D &frenet,
                                                 Point2D &cart) const {
  const auto s = frenet.x, l = frenet.y;
  if (!Available() || s < 0 || s > GetLength()) {
    return TRANSFORM_FAILED;
  }
  const auto ref_pt = lane_.GetPositionUnsafe(s);
  const auto ref_theta = lane_.GetHeadingUnsafe(s);
  cart.x = ref_pt.x - l * std::sin(ref_theta);
  cart.y = ref_pt.y + l * std::cos(ref_theta);
  return TRANSFORM_SUCCESS;
}

TRANSFORM_STATUS
FrenetCoordinateSystemNew::CartState2FrenetState(
    const CartesianState &cart_state, FrenetState &frenet_state) const {
  if (!Available()) {
    return TRANSFORM_FAILED;
  }
  const auto sl_pt = lane_.GetProjectionCoarse({cart_state.x, cart_state.y});
  const auto ref_pt = lane_.GetPositionUnsafe(sl_pt.t);
  const auto ref_theta_kappa_dkappa =
      lane_.GetHeadingCurvatureDCurvatureUnsafe(sl_pt.t);
  array<double, 3> s_condition{}, d_condition{};
  CartesianFrenetConverter::cartesian_to_frenet(
      sl_pt.t, ref_pt.x, ref_pt.y, ref_theta_kappa_dkappa[0],
      ref_theta_kappa_dkappa[1], ref_theta_kappa_dkappa[2], cart_state.x,
      cart_state.y, cart_state.speed, cart_state.acceleration, cart_state.yaw,
      cart_state.curvature, &s_condition, &d_condition);
  frenet_state.s = sl_pt.t;
  frenet_state.r = sl_pt.d;
  frenet_state.ds = s_condition[1];
  frenet_state.dds = s_condition[2];
  frenet_state.dr_ds = d_condition[1];
  frenet_state.ddr_dsds = d_condition[2];
  return TRANSFORM_SUCCESS;
}

TRANSFORM_STATUS FrenetCoordinateSystemNew::FrenetState2CartState(
    const FrenetState &frenet_state, CartesianState &cart_state) const {
  if (!Available()) {
    return TRANSFORM_FAILED;
  }
  const auto ref_pt = lane_.GetPositionUnsafe(frenet_state.s);
  const auto ref_theta_kappa_dkappa =
      lane_.GetHeadingCurvatureDCurvatureUnsafe(frenet_state.s);
  const array<double, 3> s_condition{frenet_state.s, frenet_state.ds,
                                     frenet_state.dds};
  const array<double, 3> d_condition{frenet_state.r, frenet_state.dr_ds,
                                     frenet_state.ddr_dsds};
  CartesianFrenetConverter::frenet_to_cartesian(
      frenet_state.s, ref_pt.x, ref_pt.y, ref_theta_kappa_dkappa[0],
      ref_theta_kappa_dkappa[1], ref_theta_kappa_dkappa[2], s_condition,
      d_condition, &cart_state.x, &cart_state.y, &cart_state.yaw,
      &cart_state.curvature, &cart_state.speed, &cart_state.acceleration);
  return TRANSFORM_SUCCESS;
}
}  // namespace msquare