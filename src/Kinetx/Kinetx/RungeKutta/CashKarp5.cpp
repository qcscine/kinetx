/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/RungeKutta/CashKarp5.h"
#include <Eigen/Dense>

namespace Scine {
namespace Kinetx {

CashKarp5::CashKarp5(Network& net) : RungeKutta(net) {
}

void CashKarp5::propagateY(Eigen::VectorXd& y, double& t, double& dt) const {
  /* Data and coefficients taken from
   * https://arxiv.org/pdf/1309.2710.pdf
   */

  // Cash-Karp Constants
  /* Only required for time dependent gradients
  constexpr double a2 = 1.0/5.0;
  constexpr double a3 = 3.0/10.0;
  constexpr double a4 = 3.0/5.0;
  constexpr double a5 = 1.0;
  constexpr double a6 = 7.0/8.0;
  */
  /* Intermediate step coefficients */
  constexpr const double b21 = 1.0 / 5.0;
  constexpr const double b31 = 3.0 / 40.0;
  constexpr const double b32 = 9.0 / 40.0;
  constexpr const double b41 = 3.0 / 10.0;
  constexpr const double b42 = -9.0 / 10.0;
  constexpr const double b43 = 6.0 / 5.0;
  constexpr const double b51 = -11.0 / 54.0;
  constexpr const double b52 = 5.0 / 2.0;
  constexpr const double b53 = -70.0 / 27.0;
  constexpr const double b54 = 35.0 / 27.0;
  constexpr const double b61 = 1631.0 / 55296.0;
  constexpr const double b62 = 175.0 / 512.0;
  constexpr const double b63 = 575.0 / 13824.0;
  constexpr const double b64 = 44275.0 / 110592.0;
  constexpr const double b65 = 253.0 / 4096.0;

  /* 5th order coefficients */
  constexpr const double c1 = 37.0 / 378.0;
  // constexpr const double c2 = 0.0;
  constexpr const double c3 = 250.0 / 621.0;
  constexpr const double c4 = 125.0 / 594.0;
  // constexpr const double c5 = 0.0;
  constexpr const double c6 = 512.0 / 1771.0;

  /* embedded 4th order coefficients */
  constexpr const double d1 = 2825.0 / 27648.0;
  // constexpr const double d2 = 0.0;
  constexpr const double d3 = 18575.0 / 48384.0;
  constexpr const double d4 = 13525.0 / 55296.0;
  constexpr const double d5 = 277.0 / 14336.0;
  constexpr const double d6 = 1.0 / 4.0;

  Eigen::VectorXd tmp;
  const Eigen::VectorXd k1 = dt * this->g(y);
  tmp = y + b21 * k1;
  const Eigen::VectorXd k2 = dt * this->g(tmp);
  tmp = y + b31 * k1 + b32 * k2;
  const Eigen::VectorXd k3 = dt * this->g(tmp);
  tmp = y + b41 * k1 + b42 * k2 + b43 * k3;
  const Eigen::VectorXd k4 = dt * this->g(tmp);
  tmp = y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4;
  const Eigen::VectorXd k5 = dt * this->g(tmp);
  tmp = y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5;
  const Eigen::VectorXd k6 = dt * this->g(tmp);

  // Calculate desired accuracy
  const double tol = 1.0e-10;
  Eigen::VectorXd acc = Eigen::VectorXd::Zero(y.size());
  acc.array() = y.array().abs() + k1.array().abs() + 1.0e-30;
  acc.array() = tol * acc.array();

  /* original code:
   * y4th = y + d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6;
   * y5th = y + c1 * k1 + c2 * k2 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6;
   *
   * Dropped 0 multiplications:
   */
  auto y4th = y + d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6;
  auto y5th = y + c1 * k1 + c3 * k3 + c4 * k4 + c6 * k6;

  // Error estimate
  auto err = y5th - y4th;

  bool reject = false;
  for (unsigned int i = 0; i < acc.size(); i++) {
    if (fabs(err[i]) > acc[i]) {
      reject = true;
      break;
    }
  }
  if (!reject) {
    // const double maxval = (acc.array()/err.array()).abs().maxCoeff();
    // const double incr = 0.9 * pow(maxval, 1.0/5.0);
    const double minval = (acc.array() / err.array()).abs().minCoeff();
    const double incr = 0.9 * pow(minval, 1.0 / 5.0);
    t += dt;
    y = y5th;
    dt *= fmin(incr, 5.0);
  }
  else {
    // const double maxval = (acc.array()/err.array()).abs().maxCoeff();
    // const double decr = 0.9 * pow(maxval, 1.0/4.0);
    const double minval = (acc.array() / err.array()).abs().minCoeff();
    const double decr = 0.9 * pow(minval, 1.0 / 4.0);
    // const double decr = 0.9 * minval;
    dt *= fmax(decr, 0.1);
  }
}

} /* namespace Kinetx */
} /* namespace Scine */
