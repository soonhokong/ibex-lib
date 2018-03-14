/* ============================================================================
 * D Y N I B E X - Definition of the ODE structure
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
#ifndef IBEX_ODE_H
#define IBEX_ODE_H

#include "ibex_AffineVector.h"
#include "ibex_Ctc3BCid.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcHC4.h"

namespace ibex {

const int MAX_FCT = 60;
const double tol_ctc = 1e-15;

class ivp_ode {
public:
  ivp_ode(const Function &_ydot, double _t0, const IntervalVector &_yinit) {
    nbvar = _yinit.size();

    yinit = new IntervalVector(_yinit);
    yinit_aff = new Affine3Vector(_yinit);

    ydot = new Function(_ydot);
    t0 = _t0;
    embedded_ctc = NULL;
  };

  ivp_ode(const Function &_ydot, double _t0, const Affine3Vector &_yinit_aff) {
    nbvar = _yinit_aff.size();

    yinit = new IntervalVector(_yinit_aff.itv());
    yinit_aff = new Affine3Vector(_yinit_aff);

    ydot = new Function(_ydot);
    t0 = _t0;
    embedded_ctc = NULL;
  };

  // with additive constraints
  ivp_ode(const Function &_ydot, double _t0, const IntervalVector &_yinit,
          const Array<NumConstraint> &csp) {
    nbvar = _yinit.size();

    yinit = new IntervalVector(_yinit);
    yinit_aff = new Affine3Vector(_yinit);

    ydot = new Function(_ydot);
    t0 = _t0;
    Array<NumConstraint> csp_temp(csp);
    CtcHC4 *hc4 = new CtcHC4(csp_temp, tol_ctc);
    Ctc3BCid *cid = new Ctc3BCid(*hc4);
    embedded_ctc = new CtcFixPoint(*cid, tol_ctc);
  };

  ivp_ode(const Function &_ydot, double _t0, const Affine3Vector &_yinit_aff,
          const Array<NumConstraint> &csp) {
    nbvar = _yinit_aff.size();

    yinit = new IntervalVector(_yinit_aff.itv());
    yinit_aff = new Affine3Vector(_yinit_aff);

    ydot = new Function(_ydot);
    t0 = _t0;
    Array<NumConstraint> csp_temp(csp);
    CtcHC4 *hc4 = new CtcHC4(csp_temp, tol_ctc);
    Ctc3BCid *cid = new Ctc3BCid(*hc4);
    embedded_ctc = new CtcFixPoint(*cid, tol_ctc);
  };

  void frechet_precomputing(int order) {
    edtfr = new edtree_frechet(order, ydot, nbvar);
  }

  IntervalVector compute_derivatives(int ordre, const IntervalVector &yi) {
    assert(ordre < MAX_FCT);

    if (ordre == 1) {
      IntervalMatrix J(yi.size(), yi.size());
      ydot->hansen_matrix(yi, J);
      return ydot->eval_vector(yi.mid()) + J * (yi - yi.mid());
    }
    IntervalVector der(nbvar);

    for (int j = 0; j < nbvar; j++)
      der[j] = edtfr->get_derivatives(ordre, Affine3Vector(yi), j).itv();

    return der;
  };

  Affine3Vector computeRK4derivative(const Affine3Vector &y) {
    Affine3Vector rk4_deriv(y.size(), 0);

    for (int i = 0; i < y.size(); i++) {
      rk4_deriv[i] = edtfr->lteExplicitRK4(i, y);
    }
    return rk4_deriv;
  }

  Affine3Vector computeRADAU3derivative(const Affine3Vector &y) {
    Affine3Vector radau3_deriv(y.size(), 0);

    for (int i = 0; i < y.size(); i++) {
      radau3_deriv[i] = edtfr->lteImplicitRadau3(i, y);
    }
    return radau3_deriv;
  }

  Affine3Vector computeLC3derivative(const Affine3Vector &y) {
    Affine3Vector lc3_deriv(y.size(), 0);

    for (int i = 0; i < y.size(); i++) {
      lc3_deriv[i] = edtfr->lteImplicitLobbato3c4(i, y);
    }
    return lc3_deriv;
  }

  Affine3Vector computeHEUNderivative(const Affine3Vector &y) {
    Affine3Vector heun_deriv(y.size(), 0);

    for (int i = 0; i < y.size(); i++) {
      heun_deriv[i] = edtfr->lteExplicitHeun(i, y);
    }
    return heun_deriv;
  }

  Affine3Vector computeIEULERderivative(const Affine3Vector &y) {
    Affine3Vector ieuler_deriv(y.size(), 0);

    for (int i = 0; i < y.size(); i++) {
      ieuler_deriv[i] = edtfr->lteImplicitEuler(i, y);
    }
    return ieuler_deriv;
  }

  Affine3Vector computeIMIDPOINTderivative(const Affine3Vector &y) {
    Affine3Vector imidpoint_deriv(y.size(), 0);

    for (int i = 0; i < y.size(); i++) {
      imidpoint_deriv[i] = edtfr->lteImplicitMidpoint(i, y);
    }
    return imidpoint_deriv;
  }

  Affine3Vector computeLA3derivative_aff(const Affine3Vector &y) {
    Affine3Vector la3_deriv(y.size());

    for (int i = 0; i < y.size(); i++) {
      la3_deriv[i] = edtfr->lteImplicitLobbato3a4(i, y);
    }
    return la3_deriv;
  }

  // with complete affine form
  Affine3Vector compute_derivatives_aff(int ordre, const Affine3Vector &yi) {
    assert(ordre < MAX_FCT);

    Affine3Vector der(nbvar);

    edtfr->edfr->inc_global_step();

    for (int j = 0; j < nbvar; j++)
      der[j] = edtfr->get_derivatives(ordre, yi, j);

    return der;
  };

  IntervalMatrix eval_jacobian_init() { return ydot->jacobian(*yinit); };

  void compute_derivatives(int ordre) {
    assert(ordre < MAX_FCT);

    for (int i = 0; i <= ordre; i++) {
      std::cout << "eval " << i << "^th diff" << std::endl;
      for (int j = 0; j < nbvar; j++)
        std::cout << edtfr->get_derivatives(i, Affine3Vector(*yinit), j).itv()
                  << std::endl;
    }
  };

  void compute_derivatives2(int ordre) {
    assert(ordre < MAX_FCT);

    for (int i = 0; i <= ordre; i++) {
      std::cout << "eval " << i << "^th diff" << std::endl;
      for (int j = 0; j < nbvar; j++)
        std::cout << edtfr->get_derivatives(i, Affine3Vector(*yinit), j)
                  << std::endl;
    }
  };

  ~ivp_ode() {
    delete yinit;
    delete yinit_aff;
    delete ydot;
    delete edtfr;
    delete embedded_ctc;
  };

public:
  int nbvar;
  Function *ydot;
  Affine3Vector *yinit_aff;
  double t0;
  CtcFixPoint *embedded_ctc;

private:
  IntervalVector *yinit;
  edtree_frechet *edtfr;
};

} // namespace ibex

#endif
