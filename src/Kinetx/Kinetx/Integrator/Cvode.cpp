/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/Integrator/Cvode.h"
#include "Kinetx/Integrator/RungeKutta.h"
#include <cvode/cvode.h>               // prototypes for CVODE fcts., consts.
#include <cvode/cvode_direct.h>        // access to CVDls interface
#include <nvector/nvector_serial.h>    // access to serial N_Vector
#include <sundials/sundials_math.h>    // contains the macros ABS, SUNSQR, EXP
#include <sundials/sundials_types.h>   // defs. of realtype, sunindextype
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <iomanip>
#include <iostream>

static int check_flag(void* flagvalue, const char* funcname, int opt) {
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return (1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int*)flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return 1;
  }

  return 0;
}

namespace Scine {
namespace Kinetx {

struct Cvode::Impl : public IntegratorBase {
  Impl(Network& net) : IntegratorBase(net) {
    y = N_VNew_Serial(net.nCompounds, context);
    cvode_mem = CVodeCreate(CV_BDF, context);
  }

  static int rhs_func(realtype /*t*/, N_Vector c, N_Vector g, void* user_data) {
    Network net = (*static_cast<Network*>(user_data));

    realtype* c_ptr = N_VGetArrayPointer(c); // pointer u vector data
    realtype* g_ptr = N_VGetArrayPointer(g); // pointer to udot vector data

    Eigen::Map<Eigen::VectorXd> concentrations(c_ptr, net.nCompounds);

    auto rf = std::get<0>(net.rateConstants);
    auto rb = std::get<1>(net.rateConstants);

    const auto sfT = std::get<0>(net.stoichiometryTransposed);
    const auto sbT = std::get<1>(net.stoichiometryTransposed);

    Eigen::VectorXd fRateBase = Eigen::VectorXd::Ones(net.nReactions);
    Eigen::VectorXd bRateBase = Eigen::VectorXd::Ones(net.nReactions);
    for (unsigned int iRxn = 0; iRxn < net.nReactions; ++iRxn) {
      for (Eigen::SparseMatrix<int>::InnerIterator it(sfT, iRxn); it; ++it) {
        fRateBase.data()[iRxn] *= fast_pow(concentrations[it.row()], it.value());
      }
      for (Eigen::SparseMatrix<int>::InnerIterator it(sbT, iRxn); it; ++it) {
        bRateBase.data()[iRxn] *= fast_pow(concentrations[it.row()], it.value());
      }
    }
    rf.array().colwise() *= fRateBase.array();
    rb.array().colwise() *= bRateBase.array();
    const auto totalRates = rf - rb;

    const Eigen::MatrixXd cdt = net.totalStoichiometry.cast<double>().transpose() * totalRates;
    const Eigen::VectorXd new_g = cdt.rowwise().sum();

    for (unsigned int i = 0; i < net.nCompounds; i++) {
      g_ptr[i] = new_g[i];
    }
    return 0;
  };

  void propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                 Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const {
    const Eigen::VectorXd yInitial = concentrations;
    this->propagateY(concentrations, t, dt);
    this->trackVertexAndEdgeFluxes(concentrations, yInitial, yFlux, rFlux, rForwardFlux, rBackwardFlux, dt);
  }

  void propagateY(Eigen::VectorXd& concentrations, double& tStart, double& dt) const {
    assert(_net.nCompounds == concentrations.size());
    realtype* y_ptr = N_VGetArrayPointer(y);
#pragma omp parallel for
    for (unsigned int i = 0; i < _net.nCompounds; i++) {
      y_ptr[i] = concentrations[i];
    }

    int flag;
    if (first_run) {
      flag = CVodeSetUserData(cvode_mem, &_net);

      flag = CVodeInit(cvode_mem, rhs_func, 0.0, y); // TODO
      if (check_flag(&flag, "CVodeSetUserData", 1))
        return;

      flag = CVodeSStolerances(cvode_mem, reltol, abstol);
      if (check_flag(&flag, "CVodeSStolerances", 1))
        return;

      SUNMatrix A = SUNDenseMatrix(_net.nCompounds, _net.nCompounds, context);
      if (check_flag((void*)A, "SUNDenseMatrix", 0))
        return;

      LS = SUNLinSol_Dense(y, A, context);
      if (check_flag((void*)LS, "SUNLinSol_Dense", 0))
        return;

      flag = CVodeSetLinearSolver(cvode_mem, LS, A);
      if (check_flag(&flag, "CVodeSetLinearSolver", 1))
        return;
      CVodeSetMaxStep(cvode_mem, 0.0);
      CVodeSetMaxNumSteps(cvode_mem, 10000);
      first_run = false;
    }

    flag = CVode(cvode_mem, tStart + dt, y, &tStart, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1))
      throw std::runtime_error("CVode failed.");
#pragma omp parallel for
    for (unsigned int i = 0; i < _net.nCompounds; i++) {
      concentrations[i] = y_ptr[i];
    }
  }

  ~Impl() {
    SUNLinSolFree(LS);
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
  };

  mutable bool first_run = true;
  mutable SUNLinearSolver LS;
  realtype abstol = 1e-21; // real tolerance of system
  realtype reltol = 1e-9;  // absolute tolerance of system
  sundials::Context context;
  N_Vector y;
  void* cvode_mem = NULL; // Problem dedicated memory.
};

Cvode::Cvode(Network& net) : _pimpl(new Impl(net)) {
}

Cvode::~Cvode() = default;

Eigen::MatrixXd Cvode::runIntegrationByTime(Eigen::VectorXd y, double t, double dt, Eigen::VectorXd& rFlux,
                                            Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, const double tMax,
                                            const unsigned int batchInterval, const double convergenceConcentrationChange) {
  return _pimpl->runIntegrationByTime(y, t, dt, rFlux, rForwardFlux, rBackwardFlux, tMax, batchInterval,
                                      convergenceConcentrationChange);
}

Eigen::MatrixXd Cvode::runIntegration(Eigen::VectorXd y, double t, double dt, Eigen::VectorXd& rFlux,
                                      Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux,
                                      const unsigned int batchInterval, const unsigned int nBatches,
                                      const double convergenceConcentrationChange) {
  return _pimpl->runIntegration(y, t, dt, rFlux, rForwardFlux, rBackwardFlux, batchInterval, nBatches,
                                convergenceConcentrationChange);
}
void Cvode::propagate(Eigen::VectorXd& concentrations, Eigen::VectorXd& yFlux, Eigen::VectorXd& rFlux,
                      Eigen::VectorXd& rForwardFlux, Eigen::VectorXd& rBackwardFlux, double& t, double& dt) const {
  _pimpl->propagate(concentrations, yFlux, rFlux, rForwardFlux, rBackwardFlux, t, dt);
}

} /* namespace Kinetx */
} /* namespace Scine */