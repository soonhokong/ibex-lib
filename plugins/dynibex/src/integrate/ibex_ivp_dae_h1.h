/* ============================================================================
 * D Y N I B E X - Definition of the DAE Hessenberg index 1 structure
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU
 * LGPL. See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Sep, 2015
 * Sponsored   : This research benefited from the support of the "Chair Complex
 * Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT
 * AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation
 * ParisTech and FDO ENSTA"
 * ----------------------------------------------------------------------------
 */
#ifndef IBEX_DAE_H1_H
#define IBEX_DAE_H1_H

#include "ibex_Linear.h"

namespace ibex {

class ivp_dae_h1 {
public:
  ivp_dae_h1(const Function &_ydot, const Function &_g, double _t0,
             const IntervalVector &_yinit, const IntervalVector &_xinit) {

    nbvar = _yinit.size();
    yinit = new IntervalVector(_yinit);
    yinit_aff = new Affine3Vector(_yinit);

    ydot = new Function(_ydot);
    t0 = _t0;
    embedded_ctc = NULL;
    xinit = new IntervalVector(_xinit);
    xinit_aff = new Affine3Vector(_xinit);

    g = new Function(_g);
  };

  ivp_dae_h1(const Function &_ydot, const Function &_g, double _t0,
             const Affine3Vector &_yinit_aff, const Affine3Vector &_xinit_aff) {

    nbvar = _yinit_aff.size();
    yinit = new IntervalVector(_yinit_aff.itv());
    yinit_aff = new Affine3Vector(_yinit_aff);

    ydot = new Function(_ydot);
    t0 = _t0;
    embedded_ctc = NULL;
    xinit = new IntervalVector(_xinit_aff.itv());
    xinit_aff = new Affine3Vector(_xinit_aff);

    g = new Function(_g);
  };

  ivp_dae_h1(const Function &_ydot, const Function &_g, double _t0,
             const IntervalVector &_yinit, IntervalVector &_xinit,
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
    xinit = new IntervalVector(_xinit);
    xinit_aff = new Affine3Vector(_xinit);

    g = new Function(_g);
  };

  ivp_dae_h1(const Function &_ydot, const Function &_g, double _t0,
             const Affine3Vector &_yinit_aff, const Affine3Vector &_xinit_aff,
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
    xinit = new IntervalVector(_xinit_aff.itv());
    xinit_aff = new Affine3Vector(_xinit_aff);

    g = new Function(_g);
  };

  void frechet_precomputing(int order) {
    edtfr = new edtree_frechet(order, ydot, nbvar);
  }

  Affine3Vector computeRADAU3derivative(const Affine3Vector &y) {
    Affine3Vector radau3_deriv(nbvar, 0);

    for (int i = 0; i < nbvar; i++) {
      radau3_deriv[i] = edtfr->lteImplicitRadau3(i, y);
    }
    return radau3_deriv;
  }

  Affine3Vector compute_derivatives_aff(int ordre, const Affine3Vector &yi) {

    Affine3Vector der(nbvar);

    edtfr->edfr->inc_global_step();

    for (int j = 0; j < nbvar; j++)
      der[j] = edtfr->get_derivatives(ordre, yi, j); // frechet generation auto

    return der;
  };

  Affine3Vector compute_algebraic_function(const Affine3Vector &si) {
    Affine3Eval eval_af2(*g);
    return eval_af2.eval(si).v();
  };

  IntervalVector compute_algebraic_krawczyk(IntervalVector si) {

    int n = si.size() - nbvar; // size of algebraic part

    IntervalVector si_old(si);
    IntervalVector kr(n);
    do {
      si_old = si;
      IntervalVector xi = si.subvector(nbvar, si.size() - 1);
      IntervalVector yi = si.subvector(0, nbvar - 1);

      Matrix Id(n, n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          Id[i][j] = i == j;
        }

      kr = xi.mid();

      // krawczyk parametric :
      IntervalMatrix Xinv = g->jacobian(IntervalVector(si.mid()));

      IntervalMatrix Xinv_x = Xinv.submatrix(0, n - 1, nbvar, n + nbvar - 1);

      Matrix Xinv_x_m = Xinv_x.mid();
      Matrix Xxm(n, n);
      try {
        real_inverse(Xinv_x_m, Xxm);
      } catch (SingularMatrixException &e) {
        Xxm = Id;
      }

      IntervalVector si_m(si);
      si_m.put(nbvar, kr);

      // krawczyk parametric :
      IntervalMatrix Y = g->jacobian(IntervalVector(si.mid()));

      IntervalMatrix Y_y = Y.submatrix(0, n - 1, 0, nbvar - 1);

      if (n == 1) {
        kr += -Xxm * IntervalVector(1, g->eval(si.mid()));
      } else {
        Affine3Eval eval_af2(*g);
        kr += -Xxm * eval_af2.eval(Affine3Vector(si.mid())).v().itv();
      }
      kr += (Id - Xxm * Xinv_x) * (xi - xi.mid()) - Xxm * Y_y * (yi - yi.mid());
      xi &= kr;

      si.put(nbvar, xi);

    } while (si.rel_distance(si_old) > 1e-20);

    return kr;
  };

  IntervalVector compute_algebraic_hansen(IntervalVector si) {

    int n = si.size() - nbvar; // size of algebraic part

    IntervalVector si_old(si);
    IntervalVector kr(n);
    do {
      si_old = si;
      IntervalVector xi = si.subvector(nbvar, si.size() - 1);
      IntervalVector yi = si.subvector(0, nbvar - 1);

      Matrix Id(n, n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          Id[i][j] = i == j;
        }

      kr = xi.mid();

      IntervalMatrix Xinv = g->jacobian(si);

      IntervalMatrix Xinv_x = Xinv.submatrix(0, n - 1, nbvar, n + nbvar - 1);

      Matrix Xinv_x_m = Xinv_x.mid();
      Matrix Xxm(n, n);
      try {
        real_inverse(Xinv_x_m, Xxm);
      } catch (SingularMatrixException &e) {
        Xxm = Id;
      }

      IntervalVector si_m(si);
      si_m.put(nbvar, kr);

      kr -= xi;
      kr *= -1.0;

      std::cout << "xi : " << xi << std::endl;
      gauss_seidel(Xxm * Xinv, -Xxm * g->eval_vector(si_m), kr);

      xi &= kr + xi.mid();
      std::cout << "xi ap gauss : " << xi << std::endl;
      si.put(nbvar, xi);

      kr = xi;

    } while (si.rel_distance(si_old) > 1e-20);

    return kr;
  };

  ~ivp_dae_h1() {
    delete yinit;
    delete yinit_aff;
    delete ydot;
    delete embedded_ctc;
    delete xinit;
    delete xinit_aff;
    delete g;
    delete edtfr;
  };

public:
  int nbvar;
  Function *ydot;
  Affine3Vector *yinit_aff;
  double t0;
  CtcFixPoint *embedded_ctc;
  Function *g;
  Affine3Vector *xinit_aff;

private:
  IntervalVector *yinit;
  edtree_frechet *edtfr;
  IntervalVector *xinit;
};

} // namespace ibex

#endif
