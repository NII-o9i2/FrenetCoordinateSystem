//
// Created by jojo-feng on 23-1-29.
//

#include "frenet_coordinate_system.h"
#include "eigen3/Eigen/Eigen"
#include "iostream"
#define POLYFIT_GAIN 1000.0

double PointsSquareDistance(const Point2D &p1, const Point2D &p2) {
  return pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
}
struct Segment2D {
  Point2D point[2];
};

double CrossProduct(const Segment2D &s1, const Segment2D &s2) {
  double product;
  double xx1 = s1.point[1].x - s1.point[0].x;
  double yy1 = s1.point[1].y - s1.point[0].y;
  double xx2 = s2.point[1].x - s2.point[0].x;
  double yy2 = s2.point[1].y - s2.point[0].y;
  product = xx1 * yy2 - yy1 * xx2;
  return product;
}

double DotProduct(const Segment2D &s1, const Segment2D &s2) {
  double product;
  double xx1 = s1.point[1].x - s1.point[0].x;
  double yy1 = s1.point[1].y - s1.point[0].y;
  double xx2 = s2.point[1].x - s2.point[0].x;
  double yy2 = s2.point[1].y - s2.point[0].y;
  product = xx1 * xx2 + yy1 * yy2;
  return product;
}

double Norm(const Segment2D &s1) {
  //    double product;
  double xx1 = s1.point[1].x - s1.point[0].x;
  double yy1 = s1.point[1].y - s1.point[0].y;

  return sqrt(xx1 * xx1 + yy1 * yy1);
}
const double PI = 3.1415926;
namespace msquare {

FrenetCoordinateSystem::FrenetCoordinateSystem(
    double length, const spline &x_s, const spline &y_s,
    const FrenetCoordinateSystemParameters &fcs_params)
    : m_fcs_params(fcs_params), m_length(length), m_x_s(x_s), m_y_s(y_s) {
  // If the curvature profile is not given, then compute the curvature profile
  double tmp_s = 0;
  vector<double> vec_s, vec_k;
  while (tmp_s <= length) {
    double dx = x_s.deriv(1, tmp_s);
    double dy = y_s.deriv(1, tmp_s);
    double ddx = x_s.deriv(2, tmp_s);
    double ddy = y_s.deriv(2, tmp_s);
    double tmp_k = (dx * ddy - dy * ddx) / pow(dx * dx + dy * dy, 1.5);
    if (tmp_k > 0.5) {
      tmp_k = 0.5;
    } else if (tmp_k < -0.5) {
      tmp_k = -0.5;
    }
    vec_k.push_back(tmp_k);
    vec_s.push_back(tmp_s);
    tmp_s += m_fcs_params.step_s;
  }
  m_k_s.set_points(vec_s, vec_k);
}

FrenetCoordinateSystem::FrenetCoordinateSystem(
    double length, const spline &x_s, const spline &y_s, const spline &k_s,
    const FrenetCoordinateSystemParameters &fcs_params)
    : m_fcs_params(fcs_params), m_length(length), m_x_s(x_s), m_y_s(y_s),
      m_k_s(k_s) {}

FrenetCoordinateSystem::FrenetCoordinateSystem(
    const vector<Point2D> &vec_pts,
    const FrenetCoordinateSystemParameters &fcs_params)
    : m_fcs_params(fcs_params) {
  vector<double> vec_x, vec_y;
  for (auto &pt : vec_pts) {
    vec_x.push_back(pt.x);
    vec_y.push_back(pt.y);
  }
  ConstructFrenetCoordinateSystem(vec_x, vec_y);
}
FrenetCoordinateSystem::FrenetCoordinateSystem(
    const vector<double> &vec_x, const vector<double> &vec_y,
    const FrenetCoordinateSystemParameters &fcs_params)
    : m_fcs_params(fcs_params) {
  ConstructFrenetCoordinateSystem(vec_x, vec_y);
  //    int i=0;
}
FrenetCoordinateSystem::FrenetCoordinateSystem(
    const vector<double> &vec_x, const vector<double> &vec_y,
    const FrenetCoordinateSystemParameters &fcs_params, double s_start,
    double s_end)
    : m_fcs_params(fcs_params) {
  ConstructFrenetCoordinateSystem(vec_x, vec_y, s_start, s_end);
}

void FrenetCoordinateSystem::ConstructFrenetCoordinateSystem(
    const vector<double> &vec_x, const vector<double> &vec_y,
    double s_polyfit_start, double s_polyfit_end) {

  cached_vec_s.clear();
  cached_vec_s.reserve(vec_x.size() + 1);
  cached_vec_c.clear();

  cached_vec_s.push_back(0);
  for (int i = 1; i < (int)vec_x.size(); i++) {
    cached_vec_s.push_back(cached_vec_s[i - 1] +
        sqrt(pow(vec_x[i] - vec_x[i - 1], 2) +
            pow(vec_y[i] - vec_y[i - 1], 2)));
  }
  m_x_s.set_points(cached_vec_s, vec_x);
  m_y_s.set_points(cached_vec_s, vec_y);
  m_length = cached_vec_s[(int)vec_x.size() - 1];
  vector<double> vec_k;

  s_polyfit_start = std::max(0.0, s_polyfit_start);
  double s_end = std::min(m_length, s_polyfit_end);

  for (int i = 0; i < (int)cached_vec_s.size(); i++) {
    double tmp_s = cached_vec_s[i];
    double dx = m_x_s.deriv(1, tmp_s);
    double dy = m_y_s.deriv(1, tmp_s);
    double ddx = m_x_s.deriv(2, tmp_s);
    double ddy = m_y_s.deriv(2, tmp_s);
    double tmp_k = (dx * ddy - dy * ddx) / pow(dx * dx + dy * dy, 1.5);

    if (tmp_k > 0.5) {
      tmp_k = 0.5;
    } else if (tmp_k < -0.5) {
      tmp_k = -0.5;
    }
    vec_k.push_back(tmp_k);

    if (tmp_s >= s_polyfit_start && tmp_s <= s_end) {
      if (tmp_s > cached_vec_s[0] + 30.0) {
        cached_vec_c.push_back(std::make_pair(0.0, tmp_s));
      } else {
        cached_vec_c.push_back(std::make_pair(tmp_k * POLYFIT_GAIN, tmp_s));
      }
    }
  }
  m_k_s.set_points(cached_vec_s, vec_k);
  // polyfit curvature of frenet coordinate
  curvature_fit_parameters =
      polyfit(cached_vec_c, 4); // Fitting 3rd degree polynomial:4 coefficients
  curvature_fit_parameters.push_back(POLYFIT_GAIN);
}

// Squared distance to a point on the curve (represented by the offset s)
double FrenetCoordinateSystem::SquaredDistanceToOffSetOnCurve(
    double s, const Point2D &cart0) const {
  Point2D cur_point = {m_x_s(s), m_y_s(s)};
  Point2D cart = cart0;
  return PointsSquareDistance(cart, cur_point);
}

// Gradient of the previous function
double FrenetCoordinateSystem::GradientSquaredDistanceToOffSetOnCurve(
    double s, const Point2D &cart0) const {
  return 2.0 * m_x_s.deriv(1, s) * (m_x_s(s) - cart0.x) +
      2.0 * m_y_s.deriv(1, s) * (m_y_s(s) - cart0.y);
}

// Hessian of the distance function
double FrenetCoordinateSystem::HessianSquaredDistanceToPointOnCurve(
    double s, const Point2D &cart0) const {
  return 2.0 * m_x_s.deriv(2, s) * (m_x_s(s) - cart0.x) +
      2.0 * m_y_s.deriv(2, s) * (m_y_s(s) - cart0.y) +
      2.0 * pow(m_x_s.deriv(1, s), 2) + 2.0 * pow(m_y_s.deriv(1, s), 2);
}

// Find the offset that has the minimal distance to the given point
double FrenetCoordinateSystem::FindMinDistanceOffsetOnCurve(
    const Point2D &cart0, double yaw, bool has_yaw, bool has_heuristics,
    double s_begin, double s_end) const {
  double squared_distance = std::numeric_limits<double>::infinity();
  double shortest_s = 0.0;
  double check_s_begin{0.0}, check_s_end{m_length};
  if (has_heuristics) {
    check_s_begin = std::max(0.0, s_begin);
    check_s_end = std::min(s_end, m_length);
  }
  double cur_s = check_s_begin;
  // Coarsed linear search to counter possible local minimum
  // Check if the yaw difference between the road heading and the vehicle
  // heading is too large
  // If larger than pi/2, means that the mapping is not correct
  while (cur_s <= check_s_end) {
    double temp = SquaredDistanceToOffSetOnCurve(cur_s, cart0);
    if (has_yaw) {
      if (squared_distance > temp &&
          cos(yaw - GetRefCurveHeading(cur_s)) > 0.0) {
        squared_distance = temp;
        shortest_s = cur_s;
      }
    } else {
      if (squared_distance > temp) {
        squared_distance = temp;
        shortest_s = cur_s;
      }
    }
    cur_s += m_fcs_params.coarse_step_s;
  }

  // Gradient descend to find the minimum
  double begin_s = std::max(shortest_s - m_fcs_params.coarse_step_s, 0.0);
  double end_s = std::min(shortest_s + m_fcs_params.coarse_step_s, m_length);
  double s_old = std::numeric_limits<double>::max();
  double s_new = (begin_s + end_s) / 2.0;
  int count = 0;
  while (!(std::abs(s_new - s_old) <= m_fcs_params.coord_transform_precision ||
      s_new - begin_s <= 0 || s_new - end_s >= 0 ||
      count > m_fcs_params.max_iter)) {
    s_old = s_new;
    double step_val = -(m_fcs_params.optimization_gamma /
        HessianSquaredDistanceToPointOnCurve(s_old, cart0)) *
        GradientSquaredDistanceToOffSetOnCurve(s_old, cart0);
    step_val = std::max<double>(-m_fcs_params.coarse_step_s / 10, step_val);
    step_val = std::min<double>(m_fcs_params.coarse_step_s / 10, step_val);

    s_new += step_val;
    count++;
  }
  return s_old;
}

// Cartesian coordinate to frenet coordinate
TRANSFORM_STATUS FrenetCoordinateSystem::CartCoord2FrenetCoord(
    const Point2D &cart, Point2D &frenet, bool has_heuristics, double s_begin,
    double s_end) const {
  double start_time;
  double min_point_s;
  if (has_heuristics) {
    min_point_s = FindMinDistanceOffsetOnCurve(cart, 0.0, false, has_heuristics,
                                               s_begin, s_end);
  } else {
    min_point_s = FindMinDistanceOffsetOnCurve(cart);
  }
  double distance =
      std::sqrt(SquaredDistanceToOffSetOnCurve(min_point_s, cart));
  Point2D nearest_point = {m_x_s(min_point_s), m_y_s(min_point_s)};
  Segment2D roadHeading = {
      nearest_point,
      {nearest_point.x + m_fcs_params.step_s * m_x_s.deriv(1, min_point_s),
       nearest_point.y + m_fcs_params.step_s * m_y_s.deriv(1, min_point_s)}};

  Segment2D tangent = {nearest_point, cart};

  double cross_p = CrossProduct(roadHeading, tangent);
  frenet.x = min_point_s;
  double cos_yaw_diff = std::fabs(DotProduct(roadHeading, tangent) / (Norm(roadHeading) * Norm(tangent)));
//  std::cout << "distance is "<< distance << std::endl;
//  std::cout << "cos_yaw_diff is "<< cos_yaw_diff << std::endl;
  if (distance > 1.0                  /*far away from the reference line*/
      && std::abs(cos_yaw_diff) > 0.1 /*tolerance for not completely ortho*/
      ) {
    // MSD_LOG(INFO,
    // "Frenet_Check1:Dist:%3.2f,COS:%3.4f",distance,cos_yaw_diff);
    // MSD_LOG(INFO, "Frenet_Check1:min_point_s:%3.2f",min_point_s);
    // MSD_LOG(INFO, "Frenet_Check1:NP_dx:%3.2f,NP_dy:%3.2f",
    //         nearest_point.x - cart.x,nearest_point.y - cart.y);
    return TRANSFORM_FAILED;
  }
  if (cross_p > 0)
    frenet.y = distance;
  else
    frenet.y = -distance;
  return TRANSFORM_SUCCESS;


}
/*
std::vector<double> FrenetCoordinateSystem::
CartCoord2FrenetCoord1(const Point2D& cart, Point2D& frenet) const {
    std::vector<double> res;
    double min_point_s = FindMinDistanceOffsetOnCurve(cart);
    double distance = std::sqrt(SquaredDistanceToOffSetOnCurve(min_point_s,
cart)); Point2D nearest_point = { m_x_s(min_point_s),m_y_s(min_point_s) };
    Segment2D roadHeading = { nearest_point,
                                { nearest_point.x + m_fcs_params.step_s *
m_x_s.deriv(1, min_point_s), nearest_point.y + m_fcs_params.step_s *
m_y_s.deriv(1, min_point_s)
                                }
                            };

    Segment2D tangent = { nearest_point, cart };

    double cross_p = CrossProduct(roadHeading, tangent);
    frenet.x = min_point_s;
    res.push_back(frenet.x);
    res.push_back(cross_p);
    res.push_back(m_fcs_params.step_s * m_x_s.deriv(1, min_point_s));
    res.push_back(m_fcs_params.step_s * m_y_s.deriv(1, min_point_s));
    res.push_back(cart.x - nearest_point.x);
    res.push_back(cart.y -nearest_point.y);
    res.push_back(std::abs<double>(
                DotProduct(roadHeading,tangent)
                /(Norm(roadHeading)*Norm(tangent))));

    if (cross_p > 0)
        frenet.y = distance;
    else
        frenet.y = -distance;
    return res;
}*/

// Cartesian coordinate to frenet coordinate, with yaw information
TRANSFORM_STATUS FrenetCoordinateSystem::CartCoord2FrenetCoord(
    const Point2D &cart, Point2D &frenet, double yaw) const {
  double min_point_s = FindMinDistanceOffsetOnCurve(cart, yaw, true);
  double distance = sqrt(SquaredDistanceToOffSetOnCurve(min_point_s, cart));

  Point2D nearest_point = {m_x_s(min_point_s), m_y_s(min_point_s)};
  Segment2D roadHeading = {
      nearest_point,
      {nearest_point.x + m_fcs_params.step_s * m_x_s.deriv(1, min_point_s),
       nearest_point.y + m_fcs_params.step_s * m_y_s.deriv(1, min_point_s)}};

  Segment2D tangent = {nearest_point, cart};

  double cross_p = CrossProduct(roadHeading, tangent);
  frenet.x = min_point_s;
  // Frenet coordinate ill defined if closing to the end of the curve
  if (frenet.x > m_length - m_fcs_params.step_s) {
    return TRANSFORM_FAILED;
  }
  double cos_yaw_diff = std::fabs(DotProduct(roadHeading, tangent) /
      (Norm(roadHeading) * Norm(tangent)));
  if (distance > 1.0 /*far away from the reference line*/ &&
      std::abs(cos_yaw_diff) > 0.1 /*tolerance for not completely ortho*/)
    return TRANSFORM_FAILED;
  if (cross_p > 0)
    frenet.y = distance;
  else
    frenet.y = -distance;
  return TRANSFORM_SUCCESS;
}
// Frenet coordinate to cartesian coordinate
TRANSFORM_STATUS
FrenetCoordinateSystem::FrenetCoord2CartCoord(const Point2D &frenet,
                                              Point2D &cart) const {
  if (frenet.x < 0 || frenet.x > m_length)
    return TRANSFORM_FAILED;
  Point2D ref_point;
  ref_point = {m_x_s(frenet.x), m_y_s(frenet.x)};
  // int index = FindClosestIndex(frenet0.x, 0, m_vec_s.size(), m_vec_s);
  // ref_point = { m_vec_x[index],m_vec_y[index] };

  double angle = atan2(m_y_s.deriv(1, frenet.x), m_x_s.deriv(1, frenet.x));
  double x1 = ref_point.x + cos(angle + PI / 2) * frenet.y;
  double y1 = ref_point.y + sin(angle + PI / 2) * frenet.y;
  cart.x = x1;
  cart.y = y1;
  return TRANSFORM_SUCCESS;
}

// Cartesian state to coupled frenet state
TRANSFORM_STATUS
FrenetCoordinateSystem::CartState2FrenetState(const CartesianState &cart_state,
                                              FrenetState &frenet_state) const {
  Point2D cart_coord;
  cart_coord.x = cart_state.x;
  cart_coord.y = cart_state.y;
  Point2D frenet_coord;
  TRANSFORM_STATUS stat =
      CartCoord2FrenetCoord(cart_coord, frenet_coord, cart_state.yaw);
  if (stat == TRANSFORM_FAILED)
    return TRANSFORM_FAILED;
  frenet_state.s = frenet_coord.x;
  frenet_state.r = frenet_coord.y;
  double ref_yaw = GetRefCurveHeading(frenet_state.s);
  double ref_curvature = GetRefCurveCurvature(frenet_state.s);
  double yaw_diff = cart_state.yaw - ref_yaw;
  // To get dds and ddr, we should first compute dr_ds and ddr_ds_ds
  double dr_ds = (1 - ref_curvature * frenet_state.r) * tan(yaw_diff);
  frenet_state.dr_ds = dr_ds;
  // dr
  frenet_state.dr = cart_state.speed * sin(yaw_diff);
  double dcurvature = GetRefCurveDCurvature(frenet_state.s);
  double term_a = (dcurvature * frenet_state.r + ref_curvature * dr_ds);
  double term_b = (cart_state.curvature * (1 - ref_curvature * frenet_state.r) /
      cos(yaw_diff) -
      ref_curvature);
  double term_c = (1 - ref_curvature * frenet_state.r);
  double ddr_dsds = -term_a * tan(yaw_diff) +
      term_b * term_c / (cos(yaw_diff) * cos(yaw_diff));
  frenet_state.ddr_dsds = ddr_dsds;
  frenet_state.ds =
      cart_state.speed * cos(yaw_diff) / (1 - ref_curvature * frenet_state.r);
  frenet_state.dds =
      (cart_state.acceleration -
          pow(frenet_state.ds, 2) * (term_b * term_c * tan(yaw_diff) - term_a) /
              cos(yaw_diff)) /
          (term_c / cos(yaw_diff));

  // Get ddr
  frenet_state.ddr =
      frenet_state.dds * dr_ds + pow(frenet_state.ds, 2) * ddr_dsds;
  return TRANSFORM_SUCCESS;
}

// Coupled frenet state to cartesian state
TRANSFORM_STATUS FrenetCoordinateSystem::FrenetState2CartState(
    const FrenetState &frenet_state, CartesianState &cart_state) const {
  Point2D frenet_coord;
  frenet_coord.x = frenet_state.s;
  frenet_coord.y = frenet_state.r;
  Point2D cart_coord;
  TRANSFORM_STATUS stat = FrenetCoord2CartCoord(frenet_coord, cart_coord);
  if (stat == TRANSFORM_FAILED)
    return TRANSFORM_FAILED;
  // Cartesian position
  cart_state.x = cart_coord.x;
  cart_state.y = cart_coord.y;
  double ref_yaw = GetRefCurveHeading(frenet_state.s);
  double ref_curvature = GetRefCurveCurvature(frenet_state.s);
  double dcurvature = GetRefCurveDCurvature(frenet_state.s);

  cart_state.speed =
      sqrt(pow((1 - ref_curvature * frenet_state.r) * frenet_state.ds, 2) +
          pow(frenet_state.dr, 2));
  double yaw_diff =
      atan(frenet_state.dr_ds / (1 - ref_curvature * frenet_state.r));
  cart_state.yaw = ref_yaw + yaw_diff;
  double term_a =
      dcurvature * frenet_state.r + ref_curvature * frenet_state.dr_ds;
  double term_b = 1 - ref_curvature * frenet_state.r;
  cart_state.curvature = ((frenet_state.ddr_dsds + term_a * tan(yaw_diff)) *
      (cos(yaw_diff) * cos(yaw_diff)) / term_b +
      ref_curvature) *
      cos(yaw_diff) / term_b;
  cart_state.acceleration =
      frenet_state.dds * term_b / cos(yaw_diff) +
          (pow(frenet_state.ds, 2) / cos(yaw_diff)) *
              (term_b * tan(yaw_diff) *
                  (cart_state.curvature * term_b / cos(yaw_diff) - ref_curvature) -
                  term_a);
  return TRANSFORM_SUCCESS;
}
//-----------------------------------------------------------------//
//--- Return fit coefficient of polynomial, least squares fit   ---//
//--- array: value, pair.first to pair.second                   ---//
//--- n: the number of coefficients                             ---//
//-----------------------------------------------------------------//
vector<double>
FrenetCoordinateSystem::polyfit(vector<std::pair<double, double>> &array,
                                int n) {
  int num = array.size();
  double a[num * n];
  double b[num];
  double c[n];
  std::memset(c, 0, sizeof(c));
  for (int i = 0; i < num; i++) {
    for (int j = 0; j < n; j++) {
      a[n * i + j] = pow(array[i].second, n - 1 - j);
    }
    b[i] = array[i].first;
  }
  // CvMat A = cvMat(num, n, CV_64FC1, a);
  // CvMat B = cvMat(num, 1, CV_64FC1, b);
  // CvMat X = cvMat(n, 1, CV_64FC1, c);

  // cvSolve(&A,   &B,   &X,   CV_LU);

  // vector<double> res;
  // for(int i = 0; i < n; i++){
  //     res.push_back(X.data.db[i]);
  // }

  // Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > A;
  // A = Eigen::MatrixXd::Zero( num, n);
  // Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > B;
  // B = Eigen::MatrixXd::Zero( num, 1);
  // Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > X;
  // X = Eigen::MatrixXd::Zero( n, 1);
  Eigen::MatrixXd A, B, X;
  A.resize(num, n);
  B.resize(num, 1);
  X.resize(n, 1);
  for (int i = 0; i < num; i++) {
    for (int j = 0; j < n; j++) {
      A(i, j) = a[i, j];
    }
  }
  for (int i = 0; i < num; i++) {
    B(i, 0) = b[i];
  }
  X = A.fullPivLu().solve(B);
  vector<double> res;
  for (int i = 0; i < n; i++) {
    res.push_back(X(i, 0));
  }
  return res;
}

void FrenetCoordinateSystem::VerifyCurveFit() {
  // std::ofstream outFile_;
  // outFile_.open("/home/ros/frenet_curvature_fit_output_data.csv",
  //               std::ios::app);
  // outFile_ << "s" << ',' << "cur_orign" << ',' << "cur_poly" << std::endl;

  double s = 0.0;
  double s_max = 80.0;
  double s_resolution = 0.5;
  while (s < s_max) {
    double cur = curvature_fit_parameters[0] * std::pow(s, 3) +
        curvature_fit_parameters[1] * std::pow(s, 2) +
        curvature_fit_parameters[2] * std::pow(s, 1) +
        curvature_fit_parameters[3];
    cur /= curvature_fit_parameters[4];
    double cur_orign = m_k_s.operator()(s);

    // outFile_ << s << ',' << cur_orign << ',' << cur << std::endl;
    s += s_resolution;
  }
  // std::cout << "****** curvature error max: " << cur_error_max << "percent:"
  // << cur_error_max * 100.0 / cur_base << "******" << std::endl;
  // outFile_ << -1 << ',' << -1 << ',' << -1 << ',' << -1 << std::endl;
  // outFile_.close();
}

void FrenetCoordinateSystem::Update(
    const vector<double> &vec_x, const vector<double> &vec_y,
    const FrenetCoordinateSystemParameters &fcs_params, double s_start,
    double s_end) {
  m_length = 0.0;
  m_x_s.clear();
  m_y_s.clear();
  m_k_s.clear();
  curvature_fit_parameters.clear();

  m_fcs_params = fcs_params;
  ConstructFrenetCoordinateSystem(vec_x, vec_y, s_start, s_end);
}
}