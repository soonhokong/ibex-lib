/* ============================================================================
 * D Y N I B E X - Definition of the Runge-Kutta order 4 Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU
 * LGPL. See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex
 * Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT
 * AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation
 * ParisTech and FDO ENSTA"
 * ----------------------------------------------------------------------------
 */
#ifndef IBEX_SOL_RK4_H
#define IBEX_SOL_RK4_H

#include "ibex_solution.h"
#include <iomanip>
#include <stdlib.h>

namespace ibex {

class solution_j_rk4 : public solution_j {
public:
  // method to define

  IntervalVector picard(const IntervalVector& y0, ivp_ode *_ode, int ordre) {
    return picard_tayl(y0, _ode, ordre);
  }

  // the LTE
  Affine3Vector LTE(const IntervalVector& y0, ivp_ode *_ode, double h) {
    // computation of the RK4 error
    Affine3Vector err_aff = _ode->computeRK4derivative(Affine3Vector(y0));
    err_aff *= (std::pow(h, 5) / 120);

    return err_aff;
  }

  // the factor for the next stepsize computation
  double factor_stepsize(double test) {
    return std::min(1.8, std::max(0.4, 0.95 * std::pow(1.0 / test, 0.2)));
  }

  // compute the sharpest jn+1
  int calcul_jnh(ivp_ode *_ode) {
    // with RK4 and affine form
    *box_jnh_aff = remainder_rk4(_ode);
    return 1;
  };

  // constructor
  solution_j_rk4(const Affine3Vector& _box_jn, double tn, double h,
                 ivp_ode *_ode, double a, double fac)
      : solution_j(_box_jn, tn, h, _ode, a, fac) {}

  // destructor
  ~solution_j_rk4() {}

private:
  // rk4 with remainder
  Affine3Vector remainder_rk4(ivp_ode *_ode) {
    double h = time_j.diam();

    Affine3Vector boxj1(*box_jn_aff);

    Affine3Vector k1 = _ode->compute_derivatives_aff(1, boxj1);

    Affine3Vector boxj2(k1);
    boxj2 *= (0.5 * h);

    Affine3Vector k2 = _ode->compute_derivatives_aff(1, *box_jn_aff + boxj2);

    Affine3Vector boxj3(k2);
    boxj3 *= (0.5 * h);

    Affine3Vector k3 = _ode->compute_derivatives_aff(1, *box_jn_aff + boxj3);

    Affine3Vector boxj4(k3);
    boxj4 *= (h);

    Affine3Vector k4 = _ode->compute_derivatives_aff(1, *box_jn_aff + boxj4);

    k2 *= (2.0);
    k3 *= (2.0);

    Affine3Vector int_rk4 = k1 + k2 + k3 + k4;
    int_rk4 *= (h / 6.0);
    int_rk4 += *box_jn_aff;

    return int_rk4 + *box_err_aff; // df
  };
};
} // namespace ibex

#endif
