
#include "CMSSM_semi_two_scale_initial_guesser.hpp"
#include "CMSSM_semi_two_scale_model.hpp"
#include "CMSSM_susy_two_scale_convergence_tester.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"

#include <Eigen/Core>
#include <cassert>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define MODEL model

CMSSM_semianalytic_initial_guesser<Two_scale>::CMSSM_semianalytic_initial_guesser(
   CMSSM_semianalytic<Two_scale>* model_,
   const QedQcd& oneset_,
   const CMSSM_semianalytic_low_scale_constraint<Two_scale>& low_constraint_,
   const CMSSM_semianalytic_susy_scale_constraint<Two_scale>& susy_constraint_,
   const CMSSM_semianalytic_high_scale_constraint<Two_scale>& high_constraint_,
   CMSSM_semianalytic_high_scale_constraint<Two_scale>* main_iteration_high_constraint_
)
   : Initial_guesser<Two_scale>()
   , model(model_)
   , oneset(oneset_)
   , mu_guess(0.)
   , mc_guess(0.)
   , mt_guess(0.)
   , md_guess(0.)
   , ms_guess(0.)
   , mb_guess(0.)
   , me_guess(0.)
   , mm_guess(0.)
   , mtau_guess(0.)
   , running_precision(1.0e-3)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
   , main_iteration_high_constraint(main_iteration_high_constraint_)
{
   assert(model && "CMSSM_semianalytic_initial_guesser: Error: pointer to model"
          " CMSSM_semianalytic<Two_scale> must not be zero");
}

CMSSM_semianalytic_initial_guesser<Two_scale>::~CMSSM_semianalytic_initial_guesser()
{
}

/**
 * Guesses the DR-bar model parameters by calling
 * guess_susy_parameters() and guess_soft_parameters() .
 */
void CMSSM_semianalytic_initial_guesser<Two_scale>::guess()
{
   guess_susy_parameters();
   guess_soft_parameters();
}

/**
 * Guesses the SUSY parameters at \f$m_\text{top}^\text{pole}\f$
 * from the Standard Model gauge couplings and fermion masses,
 * and the given inputs.  Threshold corrections are ignored.
 * The user-defined initial guesses at the low- and high-scales
 * are applied in this routine.  The SUSY parameters are determined
 * by doing a two-scale iteration between the two scales.
 */
void CMSSM_semianalytic_initial_guesser<Two_scale>::guess_susy_parameters()
{
   // initial guess for SUSY iteration
   QedQcd leAtMt(oneset);
   const double MZ = oneset.displayPoleMZ();
   const double MW = oneset.displayPoleMW();
   // use approximate correction (see PDG review) to convert on-shell
   // weak mixing angle to MSbar one at the Z pole mass
   const double sinThetaW2 = 1.0355 * (1.0 - Sqr(MW / MZ));
   const double mtpole = leAtMt.displayPoleMt();

   mu_guess = leAtMt.displayMass(mUp);
   mc_guess = leAtMt.displayMass(mCharm);
   mt_guess = leAtMt.displayMass(mTop) - 30.0;
   md_guess = leAtMt.displayMass(mDown);
   ms_guess = leAtMt.displayMass(mStrange);
   mb_guess = leAtMt.displayMass(mBottom);
   me_guess = leAtMt.displayMass(mElectron);
   mm_guess = leAtMt.displayMass(mMuon);
   mtau_guess = leAtMt.displayMass(mTau);

   // guess gauge couplings at mt
   const DoubleVector alpha_sm(leAtMt.getGaugeMu(mtpole, sinThetaW2));

   const double g1MSbar1lp_at_mt = sqrt(4.0 * M_PI * alpha_sm(1));
   const double g2MSbar1lp_at_mt = sqrt(4.0 * M_PI * alpha_sm(2));
   const double g3MSbar1lp_at_mt = sqrt(4.0 * M_PI * alpha_sm(3));

   model->set_g1(g1MSbar1lp_at_mt);
   model->set_g2(g2MSbar1lp_at_mt);
   model->set_g3(g3MSbar1lp_at_mt);
   model->set_scale(mtpole);

   // apply user-defined initial guess at the matching scale
   // to the SM, i.e. at Q = m_t^pole
   const auto TanBeta = INPUTPARAMETER(TanBeta);

   MODEL->set_vd(Re(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_vu(Re((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta))));
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();

   // initial guess for Mu is at high-scale
   const auto MuInput_at_MS = INPUTPARAMETER(MuInput_at_MS);
   INPUTPARAMETER(MuInput_at_MS) = false;

   guess_high_scale_parameters();

   const double low_scale_guess = low_constraint.get_initial_scale_guess();
   model->run_to(low_scale_guess, running_precision);

   // iterate until converged
   Initial_guess_low_scale_constraint low_scale_constraint
      = { model, leAtMt, g1MSbar1lp_at_mt, g2MSbar1lp_at_mt,
          g3MSbar1lp_at_mt, low_scale_guess };

   std::vector<Constraint<Two_scale>*> constraints(2);
   constraints[0] = &low_scale_constraint;
   constraints[1] = &high_constraint;

   CMSSM_susy_convergence_tester<Two_scale> convergence_tester(model, running_precision);
   Two_scale_increasing_precision running_precision_calculator(10.0, running_precision);

   RGFlow<Two_scale> initial_solver;
   initial_solver.set_convergence_tester(&convergence_tester);
   initial_solver.set_running_precision(&running_precision_calculator);
   initial_solver.add_model(model, constraints);

   initial_solver.solve();

   main_iteration_high_constraint->set_initial_scale_guess(high_constraint.get_scale());

   INPUTPARAMETER(MuInput_at_MS) = MuInput_at_MS;
}

void CMSSM_semianalytic_initial_guesser<Two_scale>::guess_low_scale_parameters(CMSSM_semianalytic<Two_scale>* model, const QedQcd& oneset, double g1, double g2, double g3)
{
   const auto TanBeta = INPUTPARAMETER(TanBeta);

   MODEL->set_vd(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta)));
   MODEL->set_vu((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta)));

   Eigen::Matrix<std::complex<double>,3,3> topDRbar(ZEROMATRIXCOMPLEX(3,3));
   topDRbar(0,0)      = oneset.displayMass(mUp);
   topDRbar(1,1)      = oneset.displayMass(mCharm);
   // note the subtraction is still in the initial iteration
   // without thresholds; afterwards it is unnecessary
   topDRbar(2,2)      = oneset.displayMass(mTop) - 30.;

const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu((Diag((1.4142135623730951*topDRbar)/vu)).real());

   Eigen::Matrix<std::complex<double>,3,3> bottomDRbar(ZEROMATRIXCOMPLEX(3,3));
   bottomDRbar(0,0)   = oneset.displayMass(mDown);
   bottomDRbar(1,1)   = oneset.displayMass(mStrange);
   bottomDRbar(2,2)   = oneset.displayMass(mBottom);

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd((Diag((1.4142135623730951*bottomDRbar)/vd)).real());

   Eigen::Matrix<std::complex<double>,3,3> electronDRbar(ZEROMATRIXCOMPLEX(3,3));
   electronDRbar(0,0) = oneset.displayMass(mElectron);
   electronDRbar(1,1) = oneset.displayMass(mMuon);
   electronDRbar(2,2) = oneset.displayMass(mTau);

   MODEL->set_Ye((Diag((1.4142135623730951*electronDRbar)/vd)).real());

   model->set_g1(g1);
   model->set_g2(g2);
   model->set_g3(g3);
}

void CMSSM_semianalytic_initial_guesser<Two_scale>::guess_high_scale_parameters()
{
   const double high_scale_guess = high_constraint.get_initial_scale_guess();

   model->run_to(high_scale_guess, running_precision);

   // apply high-scale constraint
   high_constraint.set_model(model);
   high_constraint.apply();
}

void CMSSM_semianalytic_initial_guesser<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

/**
 * Calculates the Yukawa couplings Yu of the up-type quarks
 * from the Standard Model up-type quark masses (ignoring threshold
 * corrections).
 */
void CMSSM_semianalytic_initial_guesser<Two_scale>::calculate_Yu_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> topDRbar(ZEROMATRIXCOMPLEX(3,3));
   topDRbar(0,0) = mu_guess;
   topDRbar(1,1) = mc_guess;
   // note subtraction is performed in guess_susy_parameters
   topDRbar(2,2) = mt_guess;

   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu((Diag((1.4142135623730951*topDRbar)/vu)).real());

}

/**
 * Calculates the Yukawa couplings Yd of the down-type
 * quarks from the Standard Model down-type quark masses (ignoring
 * threshold corrections).
 */
void CMSSM_semianalytic_initial_guesser<Two_scale>::calculate_Yd_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> bottomDRbar(ZEROMATRIXCOMPLEX(3,3));
   bottomDRbar(0,0) = md_guess;
   bottomDRbar(1,1) = ms_guess;
   bottomDRbar(2,2) = mb_guess;

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd((Diag((1.4142135623730951*bottomDRbar)/vd)).real());

}

/**
 * Calculates the Yukawa couplings Ye of the leptons
 * from the Standard Model down-type lepton masses (ignoring threshold
 * corrections).
 */
void CMSSM_semianalytic_initial_guesser<Two_scale>::calculate_Ye_DRbar()
{
   Eigen::Matrix<std::complex<double>,3,3> electronDRbar(ZEROMATRIXCOMPLEX(3,3));
   electronDRbar(0,0) = me_guess;
   electronDRbar(1,1) = mm_guess;
   electronDRbar(2,2) = mtau_guess;

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye((Diag((1.4142135623730951*electronDRbar)/vd)).real());

}

/**
 * Guesses the soft-breaking parameters by solving the tree-level EWSB
 * conditions at the low-scale.  The DR-bar mass spectrum is then
 * calculated.
 */
void CMSSM_semianalytic_initial_guesser<Two_scale>::guess_soft_parameters()
{
   const double low_scale_guess = low_constraint.get_initial_scale_guess();
   const double high_scale_guess = high_constraint.get_scale();

   model->run_to(low_scale_guess, running_precision);

   // apply EWSB constraint
   model->calculate_coefficients(high_scale_guess);
   model->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model->calculate_DRbar_masses();
}

void CMSSM_semianalytic_initial_guesser<Two_scale>::Initial_guess_low_scale_constraint::set_model(Two_scale_model* model_)
{
   model = cast_model<CMSSM_semianalytic<Two_scale>*>(model_);
}

} // namespace flexiblesusy
