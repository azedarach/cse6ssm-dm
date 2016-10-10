// ====================================================================
// Test implementation of low scale constraint to be used in
// semianalytic version of the two-scale algorithm
// ====================================================================

#include "CSE6SSM_semi_two_scale_low_scale_constraint.hpp"
#include "CSE6SSM_semi_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "weinberg_angle.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CSE6SSM_semianalytic<Two_scale>
#define CKM ckm
#define PMNS pmns
#define THETAW theta_w
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::CSE6SSM_semianalytic_low_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , oneset()
   , ckm()
   , pmns()
   , MWDRbar(0.)
   , MZDRbar(0.)
   , AlphaS(0.)
   , EDRbar(0.)
   , ThetaWDRbar(0.)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , self_energy_w_at_mw(0.)
   , threshold_corrections_loop_order(1)
{
   ckm << 1., 0., 0.,
          0., 1., 0.,
          0., 0., 1.;

   pmns << 1., 0., 0.,
           0., 1., 0.,
           0., 0., 1.;
}

CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::CSE6SSM_semianalytic_low_scale_constraint(CSE6SSM_semianalytic<Two_scale>* model_, const QedQcd& oneset_)
   : Constraint<Two_scale>()
   , model(model_)
   , oneset(oneset_)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , self_energy_w_at_mw(0.)
{
   initialize();
}

CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::~CSE6SSM_semianalytic_low_scale_constraint()
{
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: CSE6SSM_semianalytic_low_scale_constraint::apply():"
          " model pointer must not be zero");

   // DH: note using saved values to freeze D-terms
   double current_vs = model->get_vs();
   double current_vsb = model->get_vsb();
   double current_vphi = model->get_vphi();

   model->set_vs(model->get_saved_vs());
   model->set_vsb(model->get_saved_vsb());
   model->set_vphi(model->get_saved_vphi());

   model->calculate_DRbar_masses();

   model->set_vs(current_vs);
   model->set_vsb(current_vsb);
   model->set_vphi(current_vphi);

   update_scale();
   calculate_DRbar_gauge_couplings();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto QSInput = INPUTPARAMETER(QSInput);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_vd(Re((2*MZDRbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta)))));
   MODEL->set_vu(Re((2*MZDRbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 +
      Sqr(TanBeta)))));
   MODEL->set_QS(QSInput);

   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);

   recalculate_mw_pole();
}

const Eigen::Matrix<std::complex<double>,3,3>& CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::get_ckm() const
{
   return ckm;
}

const Eigen::Matrix<std::complex<double>,3,3>& CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::get_pmns() const
{
   return pmns;
}

double CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<CSE6SSM_semianalytic<Two_scale>*>(model_);
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::set_sm_parameters(const QedQcd& oneset_)
{
   oneset = oneset_;
}

const QedQcd& CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return oneset;
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   oneset = QedQcd();
   MWDRbar = 0.;
   MZDRbar = 0.;
   AlphaS = 0.;
   EDRbar = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   self_energy_w_at_mw = 0.;
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::initialize()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = LowEnergyConstant(MZ);

   scale = initial_scale_guess;

   MWDRbar = 0.;
   MZDRbar = 0.;
   AlphaS = 0.;
   EDRbar = 0.;
   ThetaWDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
   ckm = oneset.get_complex_ckm();
   pmns = oneset.get_complex_pmns();
   self_energy_w_at_mw = 0.;
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   scale = LowEnergyConstant(MZ);


}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_threshold_corrections()
{
   assert(oneset.displayMu() == get_scale() && "Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_threshold_corrections(): model pointer is zero");

   const double alpha_em = oneset.displayAlpha(ALPHA);
   const double alpha_s  = oneset.displayAlpha(ALPHAS);
   const double mw_pole  = oneset.displayPoleMW();
   const double mz_pole  = oneset.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds()) {
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);
      delta_alpha_s  = calculate_delta_alpha_s(alpha_s);
   }

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   MZDRbar = mz_pole;
   MWDRbar = mw_pole;

   if (model->get_thresholds()) {
      MZDRbar = model->calculate_MVZ_DRbar(mz_pole);
      MWDRbar = model->calculate_MVWm_DRbar(mw_pole);
   }

   AlphaS = alpha_s_drbar;
   EDRbar = e_drbar;
   ThetaWDRbar = calculate_theta_w(alpha_em_drbar);
}

double CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_theta_w(double alpha_em_drbar)
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_theta_w(): model pointer is zero");

   double theta_w = 0.;

   THETAW = ArcSin(Sqrt(1 - Sqr(MWDRbar)/Sqr(MZDRbar)));


   return theta_w;
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_DRbar_gauge_couplings()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_DRbar_gauge_couplings(): model pointer is zero");

   calculate_threshold_corrections();

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaWDRbar);
   new_g2 = EDRbar*Csc(ThetaWDRbar);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

}

double CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_delta_alpha_em(double alphaEm) const
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_delta_alpha_em(): model pointer is zero");

   const double currentScale = model->get_scale();
   const auto MCha = MODELPARAMETER(MCha);
   const auto MChaI = MODELPARAMETER(MChaI);
   const auto MFDX = MODELPARAMETER(MFDX);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSDX = MODELPARAMETER(MSDX);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MSHIPM = MODELPARAMETER(MSHIPM);
   const auto MSHpp = MODELPARAMETER(MSHpp);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MChaP = MODELPARAMETER(MChaP);

   const double delta_alpha_em_SM = 0.15915494309189535*alphaEm*(
      0.3333333333333333 - 1.7777777777777777*FiniteLog(Abs(MFu(2)/currentScale)))
      ;

   const double delta_alpha_em = 0.15915494309189535*alphaEm*(
      -1.3333333333333333*FiniteLog(Abs(MChaP/currentScale)) - 1.3333333333333333*
      FiniteLog(Abs(MCha(0)/currentScale)) - 1.3333333333333333*FiniteLog(Abs(MCha
      (1)/currentScale)) - 1.3333333333333333*FiniteLog(Abs(MChaI(0)/currentScale)
      ) - 1.3333333333333333*FiniteLog(Abs(MChaI(1)/currentScale)) -
      0.14814814814814814*FiniteLog(Abs(MFDX(0)/currentScale)) -
      0.14814814814814814*FiniteLog(Abs(MFDX(1)/currentScale)) -
      0.14814814814814814*FiniteLog(Abs(MFDX(2)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MHpm(1)/currentScale)) - 0.1111111111111111
      *FiniteLog(Abs(MSd(0)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(
      1)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(2)/currentScale)) -
      0.1111111111111111*FiniteLog(Abs(MSd(3)/currentScale)) - 0.1111111111111111
      *FiniteLog(Abs(MSd(4)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(
      5)/currentScale)) - 0.037037037037037035*FiniteLog(Abs(MSDX(0)/currentScale)
      ) - 0.037037037037037035*FiniteLog(Abs(MSDX(1)/currentScale)) -
      0.037037037037037035*FiniteLog(Abs(MSDX(2)/currentScale)) -
      0.037037037037037035*FiniteLog(Abs(MSDX(3)/currentScale)) -
      0.037037037037037035*FiniteLog(Abs(MSDX(4)/currentScale)) -
      0.037037037037037035*FiniteLog(Abs(MSDX(5)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(0)/currentScale)) - 0.3333333333333333*
      FiniteLog(Abs(MSe(1)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(2
      )/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(3)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(4)/currentScale)) - 0.3333333333333333*
      FiniteLog(Abs(MSe(5)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(
      MSHIPM(0)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSHIPM(1)
      /currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSHIPM(2)/currentScale))
      - 0.3333333333333333*FiniteLog(Abs(MSHIPM(3)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSHpp(0)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSHpp(1)/currentScale)) -
      0.4444444444444444*FiniteLog(Abs(MSu(0)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSu(1)/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(2
      )/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(3)/currentScale)) -
      0.4444444444444444*FiniteLog(Abs(MSu(4)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSu(5)/currentScale)));

   return delta_alpha_em + delta_alpha_em_SM;

}

double CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_delta_alpha_s(double alphaS) const
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_delta_alpha_s(): model pointer is zero");

   const double currentScale = model->get_scale();
   const auto MFDX = MODELPARAMETER(MFDX);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSDX = MODELPARAMETER(MSDX);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MGlu = MODELPARAMETER(MGlu);

   const double delta_alpha_s_SM = -0.1061032953945969*alphaS*FiniteLog(Abs(MFu
      (2)/currentScale));

   const double delta_alpha_s = 0.15915494309189535*alphaS*(0.5 - 2*FiniteLog(
      Abs(MGlu/currentScale)) - 0.6666666666666666*FiniteLog(Abs(MFDX(0)
      /currentScale)) - 0.6666666666666666*FiniteLog(Abs(MFDX(1)/currentScale)) -
      0.6666666666666666*FiniteLog(Abs(MFDX(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(5)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSDX(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSDX(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSDX(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSDX(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSDX(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSDX(5)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(5)/currentScale)));

   return delta_alpha_s + delta_alpha_s_SM;

}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_Yu_DRbar()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_Yu_DRbar(): model pointer is zero");

   Eigen::Matrix<std::complex<double>,3,3> topDRbar(ZEROMATRIXCOMPLEX(3,3));
   topDRbar(0,0)      = oneset.displayMass(mUp);
   topDRbar(1,1)      = oneset.displayMass(mCharm);
   topDRbar(2,2)      = oneset.displayMass(mTop);

   if (model->get_thresholds())
      topDRbar(2,2) = model->calculate_MFu_DRbar(oneset.displayPoleMt(), 2);

   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu((Diag((1.4142135623730951*topDRbar)/vu)).real());

}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_Yd_DRbar()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_Yd_DRbar(): model pointer is zero");

   Eigen::Matrix<std::complex<double>,3,3> bottomDRbar(ZEROMATRIXCOMPLEX(3,3));
   bottomDRbar(0,0)   = oneset.displayMass(mDown);
   bottomDRbar(1,1)   = oneset.displayMass(mStrange);
   bottomDRbar(2,2)   = oneset.displayMass(mBottom);

   if (model->get_thresholds())
      bottomDRbar(2,2) = model->calculate_MFd_DRbar(oneset.displayMass(mBottom), 2);

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd((Diag((1.4142135623730951*bottomDRbar)/vd)).real());

}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_Ye_DRbar()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_Ye_DRbar(): model pointer is zero");

   Eigen::Matrix<std::complex<double>,3,3> electronDRbar(ZEROMATRIXCOMPLEX(3,3));
   electronDRbar(0,0) = oneset.displayMass(mElectron);
   electronDRbar(1,1) = oneset.displayMass(mMuon);
   electronDRbar(2,2) = oneset.displayMass(mTau);

   if (model->get_thresholds()) {
      electronDRbar(0,0) = model->calculate_MFe_DRbar(oneset.displayMass(mElectron), 0);
      electronDRbar(1,1) = model->calculate_MFe_DRbar(oneset.displayMass(mMuon), 1);
      electronDRbar(2,2) = model->calculate_MFe_DRbar(oneset.displayMass(mTau), 2);
   }

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye((Diag((1.4142135623730951*electronDRbar)/vd)).real());

}

void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::calculate_MNeutrino_DRbar()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "calculate_MNeutrino_DRbar(): model pointer is zero");

   neutrinoDRbar.setZero();
   neutrinoDRbar(0,0) = oneset.displayNeutrinoPoleMass(1);
   neutrinoDRbar(1,1) = oneset.displayNeutrinoPoleMass(2);
   neutrinoDRbar(2,2) = oneset.displayNeutrinoPoleMass(3);
}

/**
 * Recalculates the W boson pole mass using the new gauge couplings.
 */
void CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::recalculate_mw_pole()
{
   assert(model && "CSE6SSM_semianalytic_low_scale_constraint<Two_scale>::"
          "recalculate_mw_pole(): model pointer is zero");

   if (!model->get_thresholds())
      return;


}

} // namespace flexiblesusy
