
#include "CSE6SSM_higgs_upper_bound.hpp"

#include "pv.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

CSE6SSM_higgs_upper_bound::CSE6SSM_higgs_upper_bound(const CSE6SSM_mass_eigenstates& model_)
   : model(model_)
   , upper_bound(0.)
   , include_all_gens(false)
   , include_ups(true)
   , include_downs(true)
   , include_exotics(true)
   , include_inert_singlets(true)
   , include_inert_neutral_higgs(true)
   , include_inert_charged_higgs(true)
{
   // ensure tree-level masses are used in calculations
   model.calculate_DRbar_masses();
}

CSE6SSM_higgs_upper_bound::~CSE6SSM_higgs_upper_bound()
{
}

double CSE6SSM_higgs_upper_bound::calculate_tree_level_upper_bound()
{
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double TanBeta = vu / vd;
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;
   const double Cos2Beta = Sqr(CosBeta) - Sqr(SinBeta);

   upper_bound = 0.5 * Sqr(Lambdax) * Sqr(v) * Sqr(Sin2Beta)
      + 0.25 * Sqr(Sqr(g2) + 0.6 * Sqr(g1)) * Sqr(v)
      * Sqr(Cos2Beta) + 0.025 * Sqr(g1p) * Sqr(v)
      * Sqr(-3.0 * Sqr(CosBeta) - 2.0 * Sqr(SinBeta));

   return Sqrt(upper_bound);
}

double CSE6SSM_higgs_upper_bound::calculate_one_loop_upper_bound()
{
   upper_bound = Sqr(calculate_tree_level_upper_bound());

   // additional contributions at 1-loop order
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double CosBeta = Cos(ArcTan(TanBeta));

   upper_bound += Sqr(CosBeta) * get_tadpole_vd() / vd;
   upper_bound += Sqr(SinBeta) * get_tadpole_vu() / vu;
   if (include_ups) {
      if (include_all_gens) {
         upper_bound += get_up_contribution(0);
         upper_bound += get_up_contribution(1);
      }
      upper_bound += get_up_contribution(2);
   }
   if (include_downs) {
      if (include_all_gens) {
         upper_bound += get_down_contribution(0);
         upper_bound += get_down_contribution(1);
      }
      upper_bound += get_down_contribution(2);
   }
   if (include_exotics) {
      for (unsigned gen = 0; gen < 3; ++gen)
         upper_bound += get_exotic_contribution(gen);
   }
   if (include_inert_singlets) {
      for (unsigned gen = 0; gen < 3; ++gen)
         upper_bound += get_inert_singlet_contribution(gen);
   }
   if (include_inert_neutral_higgs) {
      for (unsigned gen = 0; gen < 2; ++gen)
         upper_bound += get_inert_neutral_higgs_contribution(gen);
   }
   if (include_inert_charged_higgs) {
      for (unsigned gen = 0; gen < 2; ++gen)
         upper_bound += get_inert_charged_higgs_contribution(gen);
   }

   upper_bound = Sqrt(upper_bound);

   return upper_bound;
}

double CSE6SSM_higgs_upper_bound::get_tadpole_vd() const
{
   double result = 0.;

   if (include_ups) {
      if (include_all_gens) {
         result -= get_dV1lp_up_dvd(0);
         result -= get_dV1lp_up_dvd(1);
      }
      result -= get_dV1lp_up_dvd(2);
   }

   if (include_downs) {
      if (include_all_gens) {
         result -= get_dV1lp_down_dvd(0);
         result -= get_dV1lp_down_dvd(1);
      }
      result -= get_dV1lp_down_dvd(2);
   }

   if (include_exotics) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_exotic_dvd(gen);
   }

   if (include_inert_singlets) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_inert_singlet_dvd(gen);
   }
   if (include_inert_neutral_higgs) {
      for (unsigned gen = 0; gen < 2; ++gen)
         result -= get_dV1lp_inert_neutral_higgs_dvd(gen);
   }
   if (include_inert_charged_higgs) {
      for (unsigned gen = 0; gen < 2; ++ gen)
         result -= get_dV1lp_inert_charged_higgs_dvd(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_tadpole_vu() const
{
   double result = 0.;

   if (include_ups) {
      if (include_all_gens) {
         result -= get_dV1lp_up_dvu(0);
         result -= get_dV1lp_up_dvu(1);
      }
      result -= get_dV1lp_up_dvu(2);
   }

   if (include_downs) {
      if (include_all_gens) {
         result -= get_dV1lp_down_dvu(0);
         result -= get_dV1lp_down_dvu(1);
      }
      result -= get_dV1lp_down_dvu(2);
   }

   if (include_exotics) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_exotic_dvu(gen);
   }

   if (include_inert_singlets) {
      for (unsigned gen = 0; gen < 3; ++gen)
         result -= get_dV1lp_inert_singlet_dvu(gen);
   }
   if (include_inert_neutral_higgs) {
      for (unsigned gen = 0; gen < 2; ++ gen)
         result -= get_dV1lp_inert_neutral_higgs_dvu(gen);
   }
   if (include_inert_charged_higgs) {
      for (unsigned gen = 0; gen < 2; ++ gen)
         result -= get_dV1lp_inert_charged_higgs_dvu(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_up_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_up_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_up_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_up_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_down_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_down_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_down_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_down_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_exotic_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_exotic_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_exotic_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_exotic_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_inert_singlet_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_inert_singlet_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_inert_singlet_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_inert_singlet_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_inert_neutral_higgs_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_inert_neutral_higgs_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_inert_neutral_higgs_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_inert_neutral_higgs_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_inert_charged_higgs_contribution(unsigned gen) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double Sin2Beta = 2.0 * SinBeta * CosBeta;

   const double Delta00Pr = get_unrotated_inert_charged_higgs_contribution(gen, 0, 0);
   const double Delta01Pr = get_unrotated_inert_charged_higgs_contribution(gen, 0, 1);
   const double Delta11Pr = get_unrotated_inert_charged_higgs_contribution(gen, 1, 1);

   const double result = Sqr(CosBeta) * Delta00Pr + Sqr(SinBeta)
      * Delta11Pr + Sin2Beta * Delta01Pr;

   return result;
}


double CSE6SSM_higgs_upper_bound::get_unrotated_up_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_up_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_up_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_up_dvu_dvu(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_unrotated_down_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_down_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_down_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_down_dvu_dvu(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_unrotated_exotic_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_exotic_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_exotic_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_exotic_dvu_dvu(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_unrotated_inert_singlet_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_inert_singlet_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_inert_singlet_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_inert_singlet_dvu_dvu(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_unrotated_inert_neutral_higgs_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_inert_neutral_higgs_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_inert_neutral_higgs_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_inert_neutral_higgs_dvu_dvu(gen);
   }

   return result;
}

double CSE6SSM_higgs_upper_bound::get_unrotated_inert_charged_higgs_contribution(unsigned gen, unsigned i, unsigned j) const
{
   double result = 0.;

   if (i == 0 && j == 0) {
      result = get_d2V1lp_inert_charged_higgs_dvd_dvd(gen);
   } else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
      result = get_d2V1lp_inert_charged_higgs_dvd_dvu(gen);
   } else if (i == 1 && j == 1) {
      result = get_d2V1lp_inert_charged_higgs_dvu_dvu(gen);
   }

   return result;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_mass_matrix_Su(unsigned gen) const
{
   const double mq2 = model.get_mq2(gen, gen);
   const double mu2 = model.get_mu2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double yf = model.get_Yu(gen, gen);
   const double Tyf = model.get_TYu(gen, gen);
   const double QS = model.get_QS();

   const double mf2 = calculate_MFu2(gen);

   const double DeltaQ = 0.0125 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double Deltau = 0.0125 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2 + 0.125 * (Sqr(g2) - 0.2 *
      Sqr(g1)) * (Sqr(vd) - Sqr(vu)) + DeltaQ + mf2;
   mass_matrix_Su(0,1) = Tyf * vu / Sqrt(2.0) - 0.5 *
      Lambdax * yf * vd * vs;
   mass_matrix_Su(1,0) = mass_matrix_Su(0,1);
   mass_matrix_Su(1,1) = mu2 + 0.1 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + Deltau + mf2;

   return mass_matrix_Su;
}

Eigen::Array<double,2,1> CSE6SSM_higgs_upper_bound::calculate_MSu2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Su(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MSu2;

   MSu2(0) = 0.5 * (mass_sum - mass_diff);
   MSu2(1) = 0.5 * (mass_sum + mass_diff);

   return MSu2;
}

double CSE6SSM_higgs_upper_bound::calculate_Sin2ThetaSu(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Su(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_Cos2ThetaSu(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Su(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_MFu2(unsigned gen) const
{
   const double yf = model.get_Yu(gen, gen);
   const double vev = model.get_vu();

   return 0.5 * Sqr(yf * vev);
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_mass_matrix_Sd(unsigned gen) const
{
   const double mq2 = model.get_mq2(gen, gen);
   const double md2 = model.get_md2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double Lambdax = model.get_Lambdax();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double yf = model.get_Yd(gen, gen);
   const double Tyf = model.get_TYd(gen, gen);
   const double QS = model.get_QS();

   const double mf2 = calculate_MFd2(gen);

   const double DeltaQ = 0.0125 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double Deltad = 0.025 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2 - 0.125 * (Sqr(g2) + 0.2 *
      Sqr(g1)) * (Sqr(vd) - Sqr(vu)) + DeltaQ + mf2;
   mass_matrix_Sd(0,1) = Tyf * vd / Sqrt(2.0) - 0.5 *
      Lambdax * yf * vu * vs;
   mass_matrix_Sd(1,0) = mass_matrix_Sd(0,1);
   mass_matrix_Sd(1,1) = md2 - 0.05 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + Deltad + mf2;

   return mass_matrix_Sd;
}

Eigen::Array<double,2,1> CSE6SSM_higgs_upper_bound::calculate_MSd2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Sd(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MSd2;

   MSd2(0) = 0.5 * (mass_sum - mass_diff);
   MSd2(1) = 0.5 * (mass_sum + mass_diff);

   return MSd2;
}

double CSE6SSM_higgs_upper_bound::calculate_Sin2ThetaSd(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Sd(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_Cos2ThetaSd(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_Sd(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_MFd2(unsigned gen) const
{
   const double yf = model.get_Yd(gen, gen);
   const double vev = model.get_vd();

   return 0.5 * Sqr(yf * vev);
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_mass_matrix_SDX(unsigned gen) const
{
   const double mDx2 = model.get_mDx2(gen, gen);
   const double mDxbar2 = model.get_mDxbar2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double Lambdax = model.get_Lambdax();
   const double Sigmax = model.get_Sigmax();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();
   const double Kappa = model.get_Kappa(gen, gen);
   const double TKappa = model.get_TKappa(gen, gen);
   const double QS = model.get_QS();

   const double mf2 = calculate_MFDX2(gen);

   const double DeltaDx = -0.025 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaDxbar = -0.0375 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_SDX;

   mass_matrix_SDX(0,0) = mDx2 + 0.05 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + DeltaDx + mf2;
   mass_matrix_SDX(0,1) = TKappa * vs / Sqrt(2.0) - 0.5 *
      Kappa * (Lambdax * vd * vu + Sigmax * vphi * vsb);
   mass_matrix_SDX(1,0) = mass_matrix_SDX(0,1);
   mass_matrix_SDX(1,1) = mDxbar2 - 0.05 * Sqr(g1) * (Sqr(vd)
      - Sqr(vu)) + DeltaDxbar + mf2;

   return mass_matrix_SDX;
}

Eigen::Array<double,2,1> CSE6SSM_higgs_upper_bound::calculate_MSDX2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_SDX(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MSDX2;

   MSDX2(0) = 0.5 * (mass_sum - mass_diff);
   MSDX2(1) = 0.5 * (mass_sum + mass_diff);

   return MSDX2;
}

double CSE6SSM_higgs_upper_bound::calculate_Sin2ThetaSDX(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_SDX(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_Cos2ThetaSDX(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_SDX(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_MFDX2(unsigned gen) const
{
   const double yf = model.get_Kappa(gen, gen);
   const double vev = model.get_vs();

   return 0.5 * Sqr(yf * vev);
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_mass_matrix_HI0(unsigned gen) const
{
   const double mH1I2 = model.get_mH1I2(gen, gen);
   const double mH2I2 = model.get_mH2I2(gen, gen);
   const double Lambda = model.get_Lambda12(gen, gen);
   const double TLambda = model.get_TLambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double Sigmax = model.get_Sigmax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double mf2 = calculate_MFHI02(gen);

   const double DeltaHd = -0.0375 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaHu = -0.025 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_HI0;

   mass_matrix_HI0(0,0) = mH1I2 + 0.125 * (Sqr(g2) + 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHd + mf2;
   mass_matrix_HI0(0,1) = -TLambda * vs / Sqrt(2.0) + 0.5 * Lambda
      * (Lambdax * vd * vu + Sigmax * vphi * vsb);
   mass_matrix_HI0(1,0) = mass_matrix_HI0(0,1);
   mass_matrix_HI0(1,1) = mH2I2 - 0.125 * (Sqr(g2) + 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHu + mf2;

   return mass_matrix_HI0;
}

Eigen::Array<double,2,1> CSE6SSM_higgs_upper_bound::calculate_MHI02(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HI0(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MHI02;

   MHI02(0) = 0.5 * (mass_sum - mass_diff);
   MHI02(1) = 0.5 * (mass_sum + mass_diff);

   return MHI02;
}

double CSE6SSM_higgs_upper_bound::calculate_Sin2ThetaHI0(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HI0(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_Cos2ThetaHI0(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HI0(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_MFHI02(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double vev = model.get_vs();

   return 0.5 * Sqr(Lambda * vev);
}

double CSE6SSM_higgs_upper_bound::calculate_MSI02(unsigned gen) const
{
   const double mSI2 = model.get_mSI2(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();

   const double DeltaS = 0.0625 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   const double MSI02 = mSI2 + DeltaS;

   return MSI02;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_mass_matrix_HIPM(unsigned gen) const
{
   const double mH1I2 = model.get_mH1I2(gen, gen);
   const double mH2I2 = model.get_mH2I2(gen, gen);
   const double Lambda = model.get_Lambda12(gen, gen);
   const double TLambda = model.get_TLambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double Sigmax = model.get_Sigmax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double mf2 = calculate_MFHI02(gen);

   const double DeltaHd = -0.0375 * Sqr(g1p) * (-3.0 * Sqr(vd)
      - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));
   const double DeltaHu = -0.025 * Sqr(g1p) * (-3.0 *
      Sqr(vd) - 2.0 * Sqr(vu) + QS * Sqr(vs) - QS * Sqr(vsb));

   Eigen::Matrix<double,2,2> mass_matrix_HIPM;

   mass_matrix_HIPM(0,0) = mH1I2 - 0.125 * (Sqr(g2) - 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHd + mf2;
   mass_matrix_HIPM(0,1) = TLambda * vs / Sqrt(2.0) - 0.5 * Lambda
      * (Lambdax * vd * vu + Sigmax * vphi * vsb);
   mass_matrix_HIPM(1,0) = mass_matrix_HIPM(0,1);
   mass_matrix_HIPM(1,1) = mH2I2 + 0.125 * (Sqr(g2) - 0.6 * Sqr(g1))
      * (Sqr(vd) - Sqr(vu)) + DeltaHu + mf2;

   return mass_matrix_HIPM;
}

Eigen::Array<double,2,1> CSE6SSM_higgs_upper_bound::calculate_MHIPM2(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HIPM(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   const double mass_sum = mass_matrix(0,0) + mass_matrix(1,1);

   Eigen::Array<double,2,1> MHIPM2;

   MHIPM2(0) = 0.5 * (mass_sum - mass_diff);
   MHIPM2(1) = 0.5 * (mass_sum + mass_diff);

   return MHIPM2;
}

double CSE6SSM_higgs_upper_bound::calculate_Sin2ThetaHIPM(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HIPM(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return 2.0 * mass_matrix(0,1) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_Cos2ThetaHIPM(unsigned gen) const
{
   Eigen::Matrix<double,2,2> mass_matrix = get_mass_matrix_HIPM(gen);

   const double mass_diff = Sqrt(Sqr(mass_matrix(0,0) -
      mass_matrix(1,1)) + 4.0 * Sqr(mass_matrix(0,1)));

   return (mass_matrix(1,1) - mass_matrix(0,0)) / mass_diff;
}

double CSE6SSM_higgs_upper_bound::calculate_MFHIPM2(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double vev = model.get_vs();

   return 0.5 * Sqr(Lambda * vev);
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_Su_dvd(unsigned gen) const
{
   const double yf = model.get_Yu(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vs = model.get_vs();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) * vd - 0.075 * Sqr(g1p) * vd;
   derivatives(0,1) = -0.5 * Lambdax * yf * vs;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.2 * Sqr(g1) * vd - 0.075 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_Su_dvu(unsigned gen) const
{
   const double yf = model.get_Yu(gen, gen);
   const double Tyf = model.get_TYu(gen, gen);
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) * vu - 0.05 * Sqr(g1p) * vu
      + Sqr(yf) * vu;
   derivatives(0,1) = Tyf / Sqrt(2.0);
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.2 * Sqr(g1) * vu - 0.05 * Sqr(g1p) * vu + Sqr(yf) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_Sd_dvd(unsigned gen) const
{
   const double yf = model.get_Yd(gen, gen);
   const double Tyf = model.get_TYd(gen, gen);
   const double vd = model.get_vd();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) + 0.2 * Sqr(g1)) * vd - 0.075 * Sqr(g1p) * vd
      + Sqr(yf) * vd;
   derivatives(0,1) = Tyf / Sqrt(2.0);
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.1 * Sqr(g1) * vd - 0.15 * Sqr(g1p) * vd + Sqr(yf) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_Sd_dvu(unsigned gen) const
{
   const double yf = model.get_Yd(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) + 0.2 * Sqr(g1)) * vu - 0.05 * Sqr(g1p) * vu;
   derivatives(0,1) = -0.5 * Lambdax * yf * vs;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.1 * Sqr(g1) * vu - 0.1 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_SDX_dvd(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Kappa = model.get_Kappa(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.1 * Sqr(g1) * vd + 0.15 * Sqr(g1p) * vd;
   derivatives(0,1) = -0.5 * Kappa * Lambdax * vu;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.1 * Sqr(g1) * vd + 0.225 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_SDX_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Kappa = model.get_Kappa(gen, gen);
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.1 * Sqr(g1) * vu + 0.1 * Sqr(g1p) * vu;
   derivatives(0,1) = -0.5 * Kappa * Lambdax * vd;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.1 * Sqr(g1) * vu + 0.15 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_HI0_dvd(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vd + 0.225 * Sqr(g1p) * vd;
   derivatives(0,1) = 0.5 * Lambda * Lambdax * vu;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vd + 0.15 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_HI0_dvu(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vu + 0.15 * Sqr(g1p) * vu;
   derivatives(0,1) = 0.5 * Lambda * Lambdax * vd;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * vu + 0.1 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_HIPM_dvd(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vd + 0.225 * Sqr(g1p) * vd;
   derivatives(0,1) = -0.5 * Lambda * Lambdax * vu;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vd + 0.15 * Sqr(g1p) * vd;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_dmass_matrix_HIPM_dvu(unsigned gen) const
{
   const double Lambda = model.get_Lambda12(gen, gen);
   const double Lambdax = model.get_Lambdax();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vu + 0.15 * Sqr(g1p) * vu;
   derivatives(0,1) = -0.5 * Lambda * Lambdax * vd;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) * vu + 0.1 * Sqr(g1p) * vu;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_Su_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) - 0.075 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.2 * Sqr(g1) - 0.075 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_Su_dvd_dvu(unsigned) const
{
   return Eigen::Matrix<double,2,2>::Zero();
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_Su_dvu_dvu(unsigned gen) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double yf = model.get_Yu(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.2 * Sqr(g1)) - 0.05 * Sqr(g1p) + Sqr(yf);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.2 * Sqr(g1) - 0.05 * Sqr(g1p) + Sqr(yf);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_Sd_dvd_dvd(unsigned gen) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();
   const double yf = model.get_Yd(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) + 0.2 * Sqr(g1)) - 0.075 * Sqr(g1p) + Sqr(yf);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.1 * Sqr(g1) - 0.15 * Sqr(g1p) + Sqr(yf);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_Sd_dvd_dvu(unsigned) const
{
   return Eigen::Matrix<double,2,2>::Zero();
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_Sd_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) + 0.2 * Sqr(g1)) - 0.05 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.1 * Sqr(g1) - 0.1 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_SDX_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.1 * Sqr(g1) + 0.15 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.1 * Sqr(g1) + 0.225 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_SDX_dvd_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Kappa = model.get_Kappa(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.;
   derivatives(0,1) = -0.5 * Kappa * Lambdax;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_SDX_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.1 * Sqr(g1) + 0.1 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.1 * Sqr(g1) + 0.15 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_HI0_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.225 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_HI0_dvd_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Lambda = model.get_Lambda12(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.;
   derivatives(0,1) = 0.5 * Lambda * Lambdax;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_HI0_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) + 0.1 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_HIPM_dvd_dvd(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.225 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_HIPM_dvd_dvu(unsigned gen) const
{
   const double Lambdax = model.get_Lambdax();
   const double Lambda = model.get_Lambda12(gen, gen);

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.;
   derivatives(0,1) = -0.5 * Lambda * Lambdax;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = 0.;

   return derivatives;
}

Eigen::Matrix<double,2,2> CSE6SSM_higgs_upper_bound::get_d2mass_matrix_HIPM_dvu_dvu(unsigned) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double g1p = model.get_g1p();

   Eigen::Matrix<double,2,2> derivatives;

   derivatives(0,0) = 0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.15 * Sqr(g1p);
   derivatives(0,1) = 0.;
   derivatives(1,0) = derivatives(0,1);
   derivatives(1,1) = -0.25 * (Sqr(g2) - 0.6 * Sqr(g1)) + 0.1 * Sqr(g1p);

   return derivatives;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_up_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_Su_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_up_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_Su_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));
   const double yf = model.get_Yu(gen, gen);
   const double vu = model.get_vu();
   const double dMF2_dvu = Sqr(yf) * vu;

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double MF2 = calculate_MFu2(gen);
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));
   const double A0MF = passarino_veltman::ReA0(MF2, Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1 - 2.0 * dMF2_dvu * A0MF);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_down_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_Sd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSd(gen);
   const double Cos2Theta = calculate_Cos2ThetaSd(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));
   const double yf = model.get_Yd(gen, gen);
   const double vd = model.get_vd();
   const double dMF2_dvd = Sqr(yf) * vd;

   const Eigen::Array<double,2,1> MS2(calculate_MSd2(gen));
   const double MF2 = calculate_MFd2(gen);
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));
   const double A0MF = passarino_veltman::ReA0(MF2, Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1 - 2.0 * dMF2_dvd * A0MF);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_down_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_Sd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSd(gen);
   const double Cos2Theta = calculate_Cos2ThetaSd(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSd2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_exotic_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_SDX_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_exotic_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_SDX_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * 3.0 * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_inert_singlet_dvd(unsigned gen) const
{
   const double vd = model.get_vd();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double dMS2_dvd = -0.075 * QS * Sqr(g1p) * vd;

   const double MS2 = calculate_MSI02(gen);
   const double scale = model.get_scale();
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));

   const double result = -oneOver16PiSqr * dMS2_dvd * A0MS;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_inert_singlet_dvu(unsigned gen) const
{
   const double vu = model.get_vu();
   const double g1p = model.get_g1p();
   const double QS = model.get_QS();
   const double dMS2_dvu = -0.05 * QS * Sqr(g1p) * vu;

   const double MS2 = calculate_MSI02(gen);
   const double scale = model.get_scale();
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));

   const double result = -oneOver16PiSqr * dMS2_dvu * A0MS;

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_inert_neutral_higgs_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HI0_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_inert_neutral_higgs_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HI0_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_inert_charged_higgs_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HIPM_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);

   const double dMS20_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvd = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvd * A0MS0
      + dMS21_dvd * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_dV1lp_inert_charged_higgs_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> derivatives(get_dmass_matrix_HIPM_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);

   const double dMS20_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) + Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) - 2.0 * Sin2Theta * derivatives(0,1));
   const double dMS21_dvu = 0.5 * (derivatives(0,0) +
      derivatives(1,1) - Cos2Theta * (derivatives(0,0)
      - derivatives(1,1)) + 2.0 * Sin2Theta * derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double scale = model.get_scale();
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = -oneOver16PiSqr * (dMS20_dvu * A0MS0
      + dMS21_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_up_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_Su_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Su_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_up_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_Su_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_Su_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Su_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_up_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_Su_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Su_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSu(gen);
   const double Cos2Theta = calculate_Cos2ThetaSu(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double yf = model.get_Yu(gen, gen);
   const double vu = model.get_vu();

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMF2_dvu = Sqr(yf) * vu;

   const Eigen::Array<double,2,1> MS2(calculate_MSu2(gen));
   const double MF2 = calculate_MFu2(gen);
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MF2_dvu_dvu = Sqr(yf);

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double logMF2Q2 = Log(MF2 / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));
   const double A0MF = passarino_veltman::ReA0(MF2, Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1 - 2.0 * Sqr(dMF2_dvu) * logMF2Q2 + 2.0 * d2MF2_dvu_dvu * A0MF);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_down_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_Sd_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Sd_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSd(gen);
   const double Cos2Theta = calculate_Cos2ThetaSd(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double yf = model.get_Yd(gen, gen);
   const double vd = model.get_vd();

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMF2_dvd = Sqr(yf) * vd;

   const Eigen::Array<double,2,1> MS2(calculate_MSd2(gen));
   const double MF2 = calculate_MFd2(gen);
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MF2_dvd_dvd = Sqr(yf);

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double logMF2Q2 = Log(MF2 / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));
   const double A0MF = passarino_veltman::ReA0(MF2, Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1 - 2.0 * Sqr(dMF2_dvd) * logMF2Q2 + 2.0 * d2MF2_dvd_dvd * A0MF);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_down_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_Sd_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_Sd_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Sd_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSd(gen);
   const double Cos2Theta = calculate_Cos2ThetaSd(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSd2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_down_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_Sd_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_Sd_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSd(gen);
   const double Cos2Theta = calculate_Cos2ThetaSd(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSd2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_exotic_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_SDX_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_SDX_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_exotic_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_SDX_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_SDX_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_SDX_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_exotic_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_SDX_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_SDX_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaSDX(gen);
   const double Cos2Theta = calculate_Cos2ThetaSDX(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MSDX2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * 3.0 * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_singlet_dvd_dvd(unsigned gen) const
{
   const double g1p = model.get_g1p();
   const double vd = model.get_vd();
   const double QS = model.get_QS();

   const double dMS2_dvd = -0.075 * QS * Sqr(g1p) * vd;
   const double d2MS2_dvd_dvd = -0.075 * QS * Sqr(g1p);

   const double MS2 = calculate_MSI02(gen);

   const double scale = model.get_scale();
   const double logMS2Q2 = Log(MS2 / Sqr(scale));
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS2_dvd) * logMS2Q2
      - d2MS2_dvd_dvd * A0MS);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_singlet_dvd_dvu(unsigned gen) const
{
   const double g1p = model.get_g1p();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double QS = model.get_QS();

   const double dMS2_dvd = - 0.075 * QS * Sqr(g1p) * vd;
   const double dMS2_dvu = -0.05 * QS * Sqr(g1p) * vu;

   const double MS2 = calculate_MSI02(gen);

   const double scale = model.get_scale();
   const double logMS2Q2 = Log(MS2 / Sqr(scale));

   const double result = oneOver16PiSqr * (dMS2_dvd * dMS2_dvu *
      logMS2Q2);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_singlet_dvu_dvu(unsigned gen) const
{
   const double g1p = model.get_g1p();
   const double vu = model.get_vu();
   const double QS = model.get_QS();

   const double dMS2_dvu = -0.05 * QS * Sqr(g1p) * vu;
   const double d2MS2_dvu_dvu = -0.05 * QS * Sqr(g1p);

   const double MS2 = calculate_MSI02(gen);

   const double scale = model.get_scale();
   const double logMS2Q2 = Log(MS2 / Sqr(scale));
   const double A0MS = passarino_veltman::ReA0(MS2, Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS2_dvu) * logMS2Q2
      - d2MS2_dvu_dvu * A0MS);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_neutral_higgs_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HI0_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HI0_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_neutral_higgs_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_HI0_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_HI0_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HI0_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_neutral_higgs_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HI0_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HI0_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHI0(gen);
   const double Cos2Theta = calculate_Cos2ThetaHI0(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHI02(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_charged_higgs_dvd_dvd(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HIPM_dvd(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HIPM_dvd_dvd(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvd_dvd = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvd) * logMS20Q2
      - d2MS20_dvd_dvd * A0MS0 + Sqr(dMS21_dvd) * logMS21Q2 - d2MS21_dvd_dvd
      * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_charged_higgs_dvd_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> vd_first_derivatives(get_dmass_matrix_HIPM_dvd(gen));
   const Eigen::Matrix<double,2,2> vu_first_derivatives(get_dmass_matrix_HIPM_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HIPM_dvd_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) + Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) - 2.0 * Sin2Theta * vd_first_derivatives(0,1));
   const double dMS21_dvd = 0.5 * (vd_first_derivatives(0,0) +
      vd_first_derivatives(1,1) - Cos2Theta * (vd_first_derivatives(0,0)
      - vd_first_derivatives(1,1)) + 2.0 * Sin2Theta * vd_first_derivatives(0,1));

   const double dMS20_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) + Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) - 2.0 * Sin2Theta * vu_first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (vu_first_derivatives(0,0) +
      vu_first_derivatives(1,1) - Cos2Theta * (vu_first_derivatives(0,0)
      - vu_first_derivatives(1,1)) + 2.0 * Sin2Theta * vu_first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) + Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) - 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double d2MS21_dvd_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      (vd_first_derivatives(0,0) - vd_first_derivatives(1,1)) * (
      vu_first_derivatives(0,0) - vu_first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * vd_first_derivatives(0,1) * vu_first_derivatives(0,1)
      + Sin4Theta * (vd_first_derivatives(0,1) * (vu_first_derivatives(0,0) -
      vu_first_derivatives(1,1)) + vu_first_derivatives(0,1) * (
      vd_first_derivatives(0,0) - vd_first_derivatives(1,1)))) - Cos2Theta *
      (second_derivatives(0,0) - second_derivatives(1,1)) + 2.0 * Sin2Theta *
      second_derivatives(0,1));

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (dMS20_dvd * dMS20_dvu *
      logMS20Q2 - d2MS20_dvd_dvu * A0MS0 + dMS21_dvd * dMS21_dvu * logMS21Q2
      - d2MS21_dvd_dvu * A0MS1);

   return result;
}

double CSE6SSM_higgs_upper_bound::get_d2V1lp_inert_charged_higgs_dvu_dvu(unsigned gen) const
{
   const Eigen::Matrix<double,2,2> first_derivatives(get_dmass_matrix_HIPM_dvu(gen));
   const Eigen::Matrix<double,2,2> second_derivatives(get_d2mass_matrix_HIPM_dvu_dvu(gen));

   const double Sin2Theta = calculate_Sin2ThetaHIPM(gen);
   const double Cos2Theta = calculate_Cos2ThetaHIPM(gen);
   const double Sin4Theta = 2.0 * Sin2Theta * Cos2Theta;

   const double dMS20_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) + Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) - 2.0 * Sin2Theta * first_derivatives(0,1));
   const double dMS21_dvu = 0.5 * (first_derivatives(0,0) +
      first_derivatives(1,1) - Cos2Theta * (first_derivatives(0,0)
      - first_derivatives(1,1)) + 2.0 * Sin2Theta * first_derivatives(0,1));

   const Eigen::Array<double,2,1> MS2(calculate_MHIPM2(gen));
   const double inverse_mass_diff = 1.0 / (MS2(1) - MS2(0));

   const double d2MS20_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) - inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) + Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) - 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double d2MS21_dvu_dvu = 0.5 * (second_derivatives(0,0) +
      second_derivatives(1,1) + inverse_mass_diff * (Sqr(Sin2Theta) *
      Sqr(first_derivatives(0,0) - first_derivatives(1,1)) + 4.0 *
      Sqr(Cos2Theta) * Sqr(first_derivatives(0,1)) + 2.0 * Sin4Theta *
      first_derivatives(0,1) * (first_derivatives(0,0) -
      first_derivatives(1,1))) - Cos2Theta * (second_derivatives(0,0) -
      second_derivatives(1,1)) + 2.0 * Sin2Theta * second_derivatives(0,1)
      );

   const double scale = model.get_scale();
   const double logMS20Q2 = Log(MS2(0) / Sqr(scale));
   const double logMS21Q2 = Log(MS2(1) / Sqr(scale));
   const double A0MS0 = passarino_veltman::ReA0(MS2(0), Sqr(scale));
   const double A0MS1 = passarino_veltman::ReA0(MS2(1), Sqr(scale));

   const double result = oneOver16PiSqr * (Sqr(dMS20_dvu) * logMS20Q2
      - d2MS20_dvu_dvu * A0MS0 + Sqr(dMS21_dvu) * logMS21Q2 - d2MS21_dvu_dvu
      * A0MS1);

   return result;
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_full_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = model.self_energy_hh(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_full_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   return model.self_energy_hh(p, gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = get_unrotated_self_energy(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const std::complex<double> up = get_unrotated_up_self_energy(p, gO1, gO2);
   const std::complex<double> down = get_unrotated_down_self_energy(p, gO1, gO2);
   const std::complex<double> exotic = get_unrotated_exotic_self_energy(p, gO1, gO2);
   const std::complex<double> inert_higgs = get_unrotated_inert_higgs_self_energy(p, gO1, gO2);
   const std::complex<double> inert_charged_higgs
      = get_unrotated_inert_charged_higgs_self_energy(p, gO1, gO2);

   return up + down + exotic + inert_higgs + inert_charged_higgs;
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_up_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = get_unrotated_up_self_energy(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_down_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = get_unrotated_down_self_energy(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_exotic_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = get_unrotated_exotic_self_energy(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_inert_higgs_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = get_unrotated_inert_higgs_self_energy(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_inert_charged_higgs_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();

   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));

   Eigen::Matrix<double,5,5> rot_matrix(Eigen::Matrix<double,5,5>::Zero());

   rot_matrix(0,0) = CosBeta;
   rot_matrix(0,1) = -SinBeta;
   rot_matrix(1,0) = SinBeta;
   rot_matrix(1,1) = CosBeta;
   rot_matrix(2,2) = CosTheta;
   rot_matrix(2,3) = SinTheta;
   rot_matrix(3,2) = -SinTheta;
   rot_matrix(3,3) = CosTheta;
   rot_matrix(4,4) = 1.0;

   Eigen::Matrix<std::complex<double>,5,5> unrot_self_energy;

   for (unsigned i = 0; i < 5; ++i) {
      for (unsigned j = 0; j < 5; ++j) {
         unrot_self_energy(i,j) = get_unrotated_inert_charged_higgs_self_energy(p, i, j);
      }
   }

   Eigen::Matrix<std::complex<double>,5,5> self_energy
      = rot_matrix.transpose() * unrot_self_energy * rot_matrix;

   return self_energy(gO1, gO2);
}

double CSE6SSM_higgs_upper_bound::A0(double m) const
{
   return passarino_veltman::ReA0(m*m, Sqr(model.get_scale()));
}

double CSE6SSM_higgs_upper_bound::B0(double p, double m1, double m2) const
{
   return passarino_veltman::ReB0(p*p, m1*m1, m2*m2, Sqr(model.get_scale()));
}

double CSE6SSM_higgs_upper_bound::G0(double p, double m1, double m2) const
{
   return passarino_veltman::ReG0(p*p, m1*m1, m2*m2, Sqr(model.get_scale()));
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_up_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const Eigen::Array<double,6,1> MSu(model.get_MSu());
   const Eigen::Array<double,3,1> MFu(model.get_MFu());

   std::complex<double> result;

   std::complex<double> tmp_7663;
   std::complex<double> tmp_7664;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_7664 += A0(MSu(gI1))*model.CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_7663 += tmp_7664;
   result += (-3) * tmp_7663;
   std::complex<double> tmp_7673;
   std::complex<double> tmp_7674;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_7675;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_7675 += B0(p,MSu(gI1),MSu(gI2))*Conj(model.CpUhhconjSuSu(gO2,gI1,
            gI2))*model.CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_7674 += tmp_7675;
   }
   tmp_7673 += tmp_7674;
   result += (3) * tmp_7673;
   std::complex<double> tmp_7622;
   std::complex<double> tmp_7623;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_7624;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_7624 += (Conj(model.CpUhhbarFuFuPL(gO2,gI1,gI2))*model.CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(model.CpUhhbarFuFuPR(gO2,gI1,gI2))*model.CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_7623 += tmp_7624;
   }
   tmp_7622 += tmp_7623;
   result += (3) * tmp_7622;
   std::complex<double> tmp_7637;
   std::complex<double> tmp_7638;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_7639;
      std::complex<double> tmp_7640;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_7640 += B0(p,MFu(gI1),MFu(gI2))*(Conj(model.CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*model.CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(model.CpUhhbarFuFuPL(gO2,gI1,gI2))
            *model.CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_7639 += tmp_7640;
      tmp_7638 += (MFu(gI1)) * tmp_7639;
   }
   tmp_7637 += tmp_7638;
   result += (-6) * tmp_7637;

   return result * oneOver16PiSqr;
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_down_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const Eigen::Array<double,6,1> MSd(model.get_MSd());
   const Eigen::Array<double,3,1> MFd(model.get_MFd());

   std::complex<double> result;

   std::complex<double> tmp_7657;
   std::complex<double> tmp_7658;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_7658 += A0(MSd(gI1))*model.CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_7657 += tmp_7658;
   result += (-3) * tmp_7657;
   std::complex<double> tmp_7665;
   std::complex<double> tmp_7666;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_7667;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_7667 += B0(p,MSd(gI1),MSd(gI2))*Conj(model.CpUhhconjSdSd(gO2,gI1,
            gI2))*model.CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_7666 += tmp_7667;
   }
   tmp_7665 += tmp_7666;
   result += (3) * tmp_7665;
   std::complex<double> tmp_7614;
   std::complex<double> tmp_7615;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_7616;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_7616 += (Conj(model.CpUhhbarFdFdPL(gO2,gI1,gI2))*model.CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(model.CpUhhbarFdFdPR(gO2,gI1,gI2))*model.CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_7615 += tmp_7616;
   }
   tmp_7614 += tmp_7615;
   result += (3) * tmp_7614;
   std::complex<double> tmp_7625;
   std::complex<double> tmp_7626;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_7627;
      std::complex<double> tmp_7628;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_7628 += B0(p,MFd(gI1),MFd(gI2))*(Conj(model.CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*model.CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(model.CpUhhbarFdFdPL(gO2,gI1,gI2))
            *model.CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_7627 += tmp_7628;
      tmp_7626 += (MFd(gI1)) * tmp_7627;
   }
   tmp_7625 += tmp_7626;
   result += (-6) * tmp_7625;

   return result * oneOver16PiSqr;
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_exotic_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const Eigen::Array<double,6,1> MSDX(model.get_MSDX());
   const Eigen::Array<double,3,1> MFDX(model.get_MFDX());

   std::complex<double> result;

   std::complex<double> tmp_7659;
   std::complex<double> tmp_7660;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_7660 += A0(MSDX(gI1))*model.CpUhhUhhconjSDXSDX(gO1,gO2,gI1,gI1);
   }
   tmp_7659 += tmp_7660;
   result += (-3) * tmp_7659;
   std::complex<double> tmp_7668;
   std::complex<double> tmp_7669;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_7670;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_7670 += B0(p,MSDX(gI1),MSDX(gI2))*Conj(model.CpUhhconjSDXSDX(gO2,
            gI1,gI2))*model.CpUhhconjSDXSDX(gO1,gI1,gI2);
      }
      tmp_7669 += tmp_7670;
   }
   tmp_7668 += tmp_7669;
   result += (3) * tmp_7668;
   std::complex<double> tmp_7617;
   std::complex<double> tmp_7618;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_7619;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_7619 += (Conj(model.CpUhhbarFDXFDXPL(gO2,gI1,gI2))*
            model.CpUhhbarFDXFDXPL(gO1,gI1,gI2) + Conj(model.CpUhhbarFDXFDXPR(gO2,gI1,gI2))*
            model.CpUhhbarFDXFDXPR(gO1,gI1,gI2))*G0(p,MFDX(gI1),MFDX(gI2));
      }
      tmp_7618 += tmp_7619;
   }
   tmp_7617 += tmp_7618;
   result += (3) * tmp_7617;
   std::complex<double> tmp_7629;
   std::complex<double> tmp_7630;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_7631;
      std::complex<double> tmp_7632;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_7632 += B0(p,MFDX(gI1),MFDX(gI2))*(Conj(model.CpUhhbarFDXFDXPR(gO2
            ,gI1,gI2))*model.CpUhhbarFDXFDXPL(gO1,gI1,gI2) + Conj(model.CpUhhbarFDXFDXPL(gO2,
            gI1,gI2))*model.CpUhhbarFDXFDXPR(gO1,gI1,gI2))*MFDX(gI2);
      }
      tmp_7631 += tmp_7632;
      tmp_7630 += (MFDX(gI1)) * tmp_7631;
   }
   tmp_7629 += tmp_7630;
   result += (-6) * tmp_7629;

   return result * oneOver16PiSqr;
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_inert_higgs_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const Eigen::Array<double,7,1> MSHI0(model.get_MSHI0());
   const Eigen::Array<double,7,1> MChiI(model.get_MChiI());

   std::complex<double> result;

   std::complex<double> tmp_7676;
   std::complex<double> tmp_7677;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      tmp_7677 += A0(MSHI0(gI1))*model.CpUhhUhhconjSHI0SHI0(gO1,gO2,gI1,gI1);
   }
   tmp_7676 += tmp_7677;
   result += (-1) * tmp_7676;
   std::complex<double> tmp_7678;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      std::complex<double> tmp_7679;
      for (unsigned gI2 = 0; gI2 < 7; ++gI2) {
         tmp_7679 += B0(p,MSHI0(gI1),MSHI0(gI2))*Conj(model.CpUhhconjSHI0SHI0(
            gO2,gI1,gI2))*model.CpUhhconjSHI0SHI0(gO1,gI1,gI2);
      }
      tmp_7678 += tmp_7679;
   }
   result += tmp_7678;
   std::complex<double> tmp_7680;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      std::complex<double> tmp_7681;
      for (unsigned gI2 = 0; gI2 < 7; ++gI2) {
         tmp_7681 += B0(p,MSHI0(gI1),MSHI0(gI2))*Conj(model.CpUhhSHI0SHI0(gO2,
            gI1,gI2))*model.CpUhhSHI0SHI0(gO1,gI1,gI2);
      }
      tmp_7680 += tmp_7681;
   }
   result += tmp_7680;
   std::complex<double> tmp_7682;
   std::complex<double> tmp_7683;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      std::complex<double> tmp_7684;
      for (unsigned gI2 = 0; gI2 < 7; ++gI2) {
         tmp_7684 += (Conj(model.CpUhhChiIChiIPL(gO2,gI1,gI2))*model.CpUhhChiIChiIPL(
            gO1,gI1,gI2) + Conj(model.CpUhhChiIChiIPR(gO2,gI1,gI2))*model.CpUhhChiIChiIPR(gO1,
            gI1,gI2))*G0(p,MChiI(gI1),MChiI(gI2));
      }
      tmp_7683 += tmp_7684;
   }
   tmp_7682 += tmp_7683;
   result += (0.5) * tmp_7682;
   std::complex<double> tmp_7685;
   std::complex<double> tmp_7686;
   for (unsigned gI1 = 0; gI1 < 7; ++gI1) {
      std::complex<double> tmp_7687;
      std::complex<double> tmp_7688;
      for (unsigned gI2 = 0; gI2 < 7; ++gI2) {
         tmp_7688 += B0(p,MChiI(gI1),MChiI(gI2))*(Conj(model.CpUhhChiIChiIPR(
            gO2,gI1,gI2))*model.CpUhhChiIChiIPL(gO1,gI1,gI2) + Conj(model.CpUhhChiIChiIPL(gO2,
            gI1,gI2))*model.CpUhhChiIChiIPR(gO1,gI1,gI2))*MChiI(gI2);
      }
      tmp_7687 += tmp_7688;
      tmp_7686 += (MChiI(gI1)) * tmp_7687;
   }
   tmp_7685 += tmp_7686;
   result += (-1) * tmp_7685;

   return result * oneOver16PiSqr;
}

std::complex<double> CSE6SSM_higgs_upper_bound::get_unrotated_inert_charged_higgs_self_energy(double p, unsigned gO1, unsigned gO2) const
{
   const Eigen::Array<double,4,1> MSHIPM(model.get_MSHIPM());
   const Eigen::Array<double,2,1> MChaI(model.get_MChaI());

   std::complex<double> result;

   std::complex<double> tmp_7641;
   std::complex<double> tmp_7642;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_7642 += A0(MSHIPM(gI1))*model.CpUhhUhhconjSHIPMSHIPM(gO1,gO2,gI1,gI1);
   }
   tmp_7641 += tmp_7642;
   result += (-1) * tmp_7641;
   std::complex<double> tmp_7643;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_7644;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_7644 += B0(p,MSHIPM(gI1),MSHIPM(gI2))*Conj(
            model.CpUhhconjSHIPMSHIPM(gO2,gI1,gI2))*model.CpUhhconjSHIPMSHIPM(gO1,gI1,gI2);
      }
      tmp_7643 += tmp_7644;
   }
   result += tmp_7643;
   std::complex<double> tmp_7593;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_7594;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_7594 += (Conj(model.CpUhhbarChaIChaIPL(gO2,gI1,gI2))*
            model.CpUhhbarChaIChaIPL(gO1,gI1,gI2) + Conj(model.CpUhhbarChaIChaIPR(gO2,gI1,gI2)
            )*model.CpUhhbarChaIChaIPR(gO1,gI1,gI2))*G0(p,MChaI(gI1),MChaI(gI2));
      }
      tmp_7593 += tmp_7594;
   }
   result += tmp_7593;
   std::complex<double> tmp_7602;
   std::complex<double> tmp_7603;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_7604;
      std::complex<double> tmp_7605;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_7605 += B0(p,MChaI(gI1),MChaI(gI2))*(Conj(model.CpUhhbarChaIChaIPR
            (gO2,gI1,gI2))*model.CpUhhbarChaIChaIPL(gO1,gI1,gI2) + Conj(
            model.CpUhhbarChaIChaIPL(gO2,gI1,gI2))*model.CpUhhbarChaIChaIPR(gO1,gI1,gI2))*
            MChaI(gI2);
      }
      tmp_7604 += tmp_7605;
      tmp_7603 += (MChaI(gI1)) * tmp_7604;
   }
   tmp_7602 += tmp_7603;
   result += (-2) * tmp_7602;

   return result * oneOver16PiSqr;
}

} // namespace flexiblesusy
