/* ============================================================================
 * D Y N I B E X - Definition of the Integrate main include file
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
#ifndef IBEX_INTEG_H
#define IBEX_INTEG_H

#include "ibex_edfrechet.h"
#include "ibex_edtree_frechet.h"
#include "ibex_ivp_ode.h"
#include "ibex_logger.h"

#include "ibex_ivp_dae_h1.h"
#include "ibex_solution_generic.h"

#include "ibex_solution_dae.h"

#include "ibex_simulation.h"
#include "ibex_solution.h"
#include "ibex_solution_HEUN.h"
#include "ibex_solution_LA3.h"
#include "ibex_solution_LC3.h"
#include "ibex_solution_RADAU3.h"
#include "ibex_solution_RK4.h"
#include "ibex_solution_TAYLOR4.h"
#include "ibex_solution_iEULER.h"
#include "ibex_solution_iMIDPOINT.h"

namespace ibex {
const double hmin = 1e-8;
const double hmax = 0.9;
}

#endif
