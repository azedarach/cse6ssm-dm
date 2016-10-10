// ====================================================================
// Executable for doing simple scans using the semianalytic algorithm
// ====================================================================

#include "CSE6SSM_semi_two_scale_input_parameters.hpp"
#include "CSE6SSM_semianalytic_spectrum_generator.hpp"
#include "CSE6SSM_semi_two_scale_model_slha.hpp"

#include "command_line_options.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <Eigen/Core>

#include <iostream>
#include <iomanip>
#include <cstring>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: scan_semianalytic_CSE6SSM.x [options]\n"
      "Options:\n"
      "  --m12=<value>\n"
      "  --Azero=<value>\n"
      "  --TanBeta=<value>\n"
      "  --sInput=<value>\n"
      "  --QSInput=<value>\n"
      "  --SigmaLInput=<value>\n"
      "  --KappaPrInput=<value>\n"
      "  --SigmaxInput=<value>\n"
      "  --LambdaxInput=<value>\n"
      "  --MuPrInput=<value>\n"
      "  --MuPhiInput=<value>\n"
      "  --BMuPrInput=<value>\n"
      "  --BMuPhiInput=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if(Command_line_options::get_parameter_value(option, "--m12=", input.m12))
         continue;

      if(Command_line_options::get_parameter_value(option, "--Azero=", input.Azero))
         continue;

      if(Command_line_options::get_parameter_value(option, "--TanBeta=", input.TanBeta))
         continue;

      if(Command_line_options::get_parameter_value(option, "--sInput=", input.sInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--QSInput=", input.QSInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--SigmaLInput=", input.SigmaLInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--KappaPrInput=", input.KappaPrInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--SigmaxInput=", input.SigmaxInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--LambdaxInput=", input.LambdaxInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuPrInput=", input.MuPrInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--MuPhiInput=", input.MuPhiInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BMuPrInput=", input.BMuPrInput))
         continue;

      if(Command_line_options::get_parameter_value(option, "--BMuPhiInput=", input.BMuPhiInput))
         continue;


      if (strcmp(option,"--help") == 0 || strcmp(option,"-h") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

void set_default_parameters(CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   input.m12 = 3500.;
   input.TanBeta = 10.0;
   input.Azero = 1100.;
   
   input.sInput = 650000.0; // GeV
   input.QSInput = 5.;
   
   input.hEInput(0,0) = 0.;
   input.hEInput(0,1) = 0.;
   input.hEInput(1,0) = 0.;
   input.hEInput(1,1) = 0.;
   input.hEInput(2,0) = 0.;
   input.hEInput(2,1) = 0.;
   
   input.SigmaLInput = 3.0e-2;
   input.KappaPrInput = 1.0e-2;
   input.SigmaxInput = 1.0e-2;
   
   input.gDInput(0,0) = 0.;
   input.gDInput(0,1) = 0.;
   input.gDInput(0,2) = 0.;
   input.gDInput(1,0) = 0.;
   input.gDInput(1,1) = 0.;
   input.gDInput(1,2) = 0.;
   input.gDInput(2,0) = 0.;
   input.gDInput(2,1) = 0.;
   input.gDInput(2,2) = 0.;
   
   input.KappaInput(0,0) = 1.0e-3;
   input.KappaInput(0,1) = 0.;
   input.KappaInput(0,2) = 0.;
   input.KappaInput(1,0) = 0.;
   input.KappaInput(1,1) = 1.0e-3;
   input.KappaInput(1,2) = 0.;
   input.KappaInput(2,0) = 0.;
   input.KappaInput(2,1) = 0.;
   input.KappaInput(2,2) = 1.0e-3;
   
   input.Lambda12Input(0,0) = 6.0e-3;
   input.Lambda12Input(0,1) = 0.;
   input.Lambda12Input(1,0) = 0.;
   input.Lambda12Input(1,1) = 6.0e-3;
   
   input.fuInput(0,0) = 1.0e-7;
   input.fuInput(0,1) = 0.;
   input.fuInput(1,0) = 0.;
   input.fuInput(1,1) = 1.0e-7;
   input.fuInput(2,0) = 1.0e-7;
   input.fuInput(2,1) = 0.;
   
   input.fdInput(0,0) = 1.0e-7;
   input.fdInput(0,1) = 0.;
   input.fdInput(1,0) = 0.;
   input.fdInput(1,1) = 1.0e-7;
   input.fdInput(2,0) = 0.;
   input.fdInput(2,1) = 1.0e-7;
   
   input.MuPrInput = 1.0e4;
   input.MuPhiInput = 0.;
   input.BMuPrInput = 1.0e4;
   input.BMuPhiInput = 0.;
}

Eigen::Matrix<double,5,5> get_unrotated_mass_matrix_hh(CSE6SSM_mass_eigenstates model)
{
   Eigen::Array<double,5,1> soft_masses = model.get_ewsb_tree_level_soft_masses();

   model.set_mHd2(soft_masses(0));
   model.set_mHu2(soft_masses(1));
   model.set_ms2(soft_masses(2));
   model.set_msbar2(soft_masses(3));
   model.set_mphi2(soft_masses(4));

   return model.get_mass_matrix_hh();
}

double get_singlet_mixing_contribution_1(const CSE6SSM_mass_eigenstates& model)
{
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double Sigmax = model.get_Sigmax();

   const double s = Sqrt(Sqr(vs) + Sqr(vsb));
   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));
   const double CosTwoTheta = Sqr(CosTheta) - Sqr(SinTheta);
   const double SinTwoTheta = 2.0 * SinTheta * CosTheta;

   const double result = 0.5 * Sqr(Sigmax) * Sqr(s) * SinTwoTheta *
      CosTwoTheta;

   return result;
}

double get_singlet_mixing_contribution_2(const CSE6SSM_mass_eigenstates& model)
{
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double TSigmax = model.get_TSigmax();

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));
   const double CosTwoTheta = Sqr(CosTheta) - Sqr(SinTheta);

   const double result = -Sqrt(2.0) * TSigmax * vphi * CosTwoTheta;

   return result;
}

double get_singlet_mixing_contribution_3(const CSE6SSM_mass_eigenstates& model)
{
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double Sigmax = model.get_Sigmax();
   const double KappaPr = model.get_KappaPr();

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));
   const double CosTwoTheta = Sqr(CosTheta) - Sqr(SinTheta);

   const double result = -CosTwoTheta * KappaPr * Sigmax * Sqr(vphi);

   return result;
}

double get_singlet_mixing_contribution_4(const CSE6SSM_mass_eigenstates& model)
{
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double Sigmax = model.get_Sigmax();
   const double MuPhi = model.get_MuPhi();

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));
   const double CosTwoTheta = Sqr(CosTheta) - Sqr(SinTheta);

   const double result = -CosTwoTheta * Sqrt(2.0) * Sigmax * MuPhi * vphi;

   return result;
}

double get_singlet_mixing_contribution_5(const CSE6SSM_mass_eigenstates& model)
{
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double Sigmax = model.get_Sigmax();
   const double XiF = model.get_XiF();

   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));
   const double SinTheta = Sin(ArcTan(TanTheta));
   const double CosTwoTheta = Sqr(CosTheta) - Sqr(SinTheta);

   const double result = -CosTwoTheta * 2.0 * Sigmax * XiF;

   return result;
}

double get_singlet_mixing_contribution_6(const CSE6SSM_mass_eigenstates& model)
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double TLambdax = model.get_TLambdax();

   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double s = Sqrt(Sqr(vs) + Sqr(vsb));
   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double SinTwoBeta = 2.0 * SinBeta * CosBeta;
   const double TanTheta = vsb / vs;
   const double SinTheta = Sin(ArcTan(TanTheta));

   const double result = 0.5 * TLambdax * Sqr(v) * SinTheta * SinTwoBeta
      / (Sqrt(2.0) * s);

   return result;
}

double get_singlet_mixing_contribution_7(const CSE6SSM_mass_eigenstates& model)
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double vs = model.get_vs();
   const double vsb = model.get_vsb();
   const double vphi = model.get_vphi();
   const double Sigmax = model.get_Sigmax();
   const double Lambdax = model.get_Lambdax();

   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double s = Sqrt(Sqr(vs) + Sqr(vsb));
   const double TanBeta = vu / vd;
   const double CosBeta = Cos(ArcTan(TanBeta));
   const double SinBeta = Sin(ArcTan(TanBeta));
   const double SinTwoBeta = 2.0 * SinBeta * CosBeta;
   const double TanTheta = vsb / vs;
   const double CosTheta = Cos(ArcTan(TanTheta));

   const double result = 0.25 * Lambdax * Sigmax * vphi * Sqr(v) *
      CosTheta * SinTwoBeta / s;

   return result;
}

} // namespace flexiblesusy

int main(int argc, char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   CSE6SSM_semianalytic_input_parameters<algorithm_type> input;
   set_default_parameters(input);
   set_command_line_parameters(argc, argv, input);

   input.LambdaxInput = 7.2e-4;

   QedQcd oneset;

   // set SM inputs to match default SLHA inputs
   oneset.setAlpha(ALPHA, 1.0 / 127.934);
   oneset.setFermiConstant(1.16637e-5);
   oneset.setAlpha(ALPHAS, 0.1176);
   oneset.setPoleMZ(91.1876);
   oneset.setMbMb(4.2);
   oneset.setPoleMt(173.3);
   oneset.setMass(mTau, 1.777);
   oneset.setPoleMtau(1.777);
   oneset.setNeutrinoPoleMass(3, 0.);
   oneset.setMass(mElectron, 5.10998902e-4);
   oneset.setNeutrinoPoleMass(1, 0.);
   oneset.setMass(mMuon, 1.05658357e-1);
   oneset.setNeutrinoPoleMass(2, 0.);
   oneset.setMass(mDown, 4.75e-3);
   oneset.setMass(mUp, 2.4e-3);
   oneset.setMass(mStrange, 1.04e-1);
   oneset.setMass(mCharm, 1.27);

   oneset.toMz();

   const double range_lower = -10000.0;
   const double range_upper = 10000.0;
   const std::size_t range_length = 800;

   const std::vector<double> range(float_range(range_lower, range_upper, range_length));

   std::cout << "# "
             << std::setw(18) << std::left << "m12/GeV" << ' '
             << std::setw(18) << std::left << "Azero/GeV" << ' '
             << std::setw(18) << std::left << "KappaPrInput" << ' '
             << std::setw(18) << std::left << "LambdaxInput" << ' '
             << std::setw(18) << std::left << "SigmaxInput" << ' '
             << std::setw(18) << std::left << "MuPhiInput/GeV" << ' '
             << std::setw(18) << std::left << "vs/GeV" << ' '
             << std::setw(18) << std::left << "vsb/GeV" << ' '
             << std::setw(18) << std::left << "vphi/GeV"
             << std::setw(18) << std::left << "Mhh(0)/GeV" << ' '
             << std::setw(18) << std::left << "Mhh(1)/GeV" << ' '
             << std::setw(18) << std::left << "MhhDRbar(0)/GeV" << ' '
             << std::setw(18) << std::left << "MhhDRbar(1)/GeV" << ' '
             << std::setw(18) << std::left << "MhhDRbar(2)/GeV" << ' '
             << std::setw(18) << std::left << "MhhDRbar(3)/GeV" << ' '
             << std::setw(18) << std::left << "MhhDRbar(4)/GeV" << ' '
             << std::setw(18) << std::left << "Mh(0,0)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(0,1)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(0,2)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(0,3)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(0,4)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(1,1)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(1,2)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(1,3)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(1,4)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(2,2)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(2,3)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(2,4)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(3,3)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(3,4)/GeV^2" << ' '
             << std::setw(18) << std::left << "Mh(4,4)/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T1/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T2/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T3/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T4/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T5/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T6/GeV^2" << ' '
             // << std::setw(18) << std::left << "Mh(0,1)T7/GeV^2" << ' '
             << std::setw(18) << std::left << "ZH(0,0)" << ' '
             << std::setw(18) << std::left << "ZH(0,1)" << ' '
             << std::setw(18) << std::left << "ZH(0,2)" << ' '
             << std::setw(18) << std::left << "ZH(0,3)" << ' '
             << std::setw(18) << std::left << "ZH(0,4)" << ' '
             << std::setw(18) << std::left << "ZH(1,0)" << ' '
             << std::setw(18) << std::left << "ZH(1,1)" << ' '
             << std::setw(18) << std::left << "ZH(1,2)" << ' '
             << std::setw(18) << std::left << "ZH(1,3)" << ' '
             << std::setw(18) << std::left << "ZH(1,4)" << ' '
             << std::setw(18) << std::left << "ZH(2,0)" << ' '
             << std::setw(18) << std::left << "ZH(2,1)" << ' '
             << std::setw(18) << std::left << "ZH(2,2)" << ' '
             << std::setw(18) << std::left << "ZH(2,3)" << ' '
             << std::setw(18) << std::left << "ZH(2,4)" << ' '
             << std::setw(18) << std::left << "ZH(3,0)" << ' '
             << std::setw(18) << std::left << "ZH(3,1)" << ' '
             << std::setw(18) << std::left << "ZH(3,2)" << ' '
             << std::setw(18) << std::left << "ZH(3,3)" << ' '
             << std::setw(18) << std::left << "ZH(3,4)" << ' '
             << std::setw(18) << std::left << "ZH(4,0)" << ' '
             << std::setw(18) << std::left << "ZH(4,1)" << ' '
             << std::setw(18) << std::left << "ZH(4,2)" << ' '
             << std::setw(18) << std::left << "ZH(4,3)" << ' '
             << std::setw(18) << std::left << "ZH(4,4)" << ' '
             << std::setw(18) << std::left << "error"
             << '\n';

   for (std::vector<double>::const_iterator it = range.begin(),
           end = range.end(); it != end; ++it) {
      input.Azero = *it;

      CSE6SSM_semianalytic_spectrum_generator<algorithm_type>
         spectrum_generator;
      spectrum_generator.set_precision_goal(1.0e-4);
      spectrum_generator.set_max_iterations(0);         // 0 == automatic
      spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
      spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale
      spectrum_generator.set_ewsb_loop_order(2);
      spectrum_generator.set_pole_mass_loop_order(2);
      spectrum_generator.set_beta_loop_order(2);
      spectrum_generator.set_threshold_corrections_loop_order(1);

      spectrum_generator.run(oneset, input);

      const CSE6SSM_semianalytic_slha<algorithm_type> model(spectrum_generator.get_model());
      const CSE6SSM_physical& pole_masses = model.get_physical_slha();
      const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();
      const double higgs = pole_masses.Mhh(0);
      const double second_higgs = pole_masses.Mhh(1);
      const Eigen::Matrix<double,5,5> mass_matrix_hh = get_unrotated_mass_matrix_hh(model);
      // const double term_1 = get_singlet_mixing_contribution_1(model);
      // const double term_2 = get_singlet_mixing_contribution_2(model);
      // const double term_3 = get_singlet_mixing_contribution_3(model);
      // const double term_4 = get_singlet_mixing_contribution_4(model);
      // const double term_5 = get_singlet_mixing_contribution_5(model);
      // const double term_6 = get_singlet_mixing_contribution_6(model);
      // const double term_7 = get_singlet_mixing_contribution_7(model);
      Eigen::Array<double,5,1> Mhh_DRbar;
      Eigen::Matrix<double,5,5> ZH;
      fs_diagonalize_hermitian(mass_matrix_hh, Mhh_DRbar, ZH); 
      for (std::size_t i = 0; i < 5; ++i)
         Mhh_DRbar(i) = SignedAbsSqrt(Mhh_DRbar(i));
      const bool error = problems.have_problem();

      std::cout << std::scientific << std::setprecision(12) << " "
                << std::setw(18) << std::left << input.m12 << ' '
                << std::setw(18) << std::left << input.Azero << ' '
                << std::setw(18) << std::left << input.KappaPrInput << ' '
                << std::setw(18) << std::left << input.LambdaxInput << ' '
                << std::setw(18) << std::left << input.SigmaxInput << ' '
                << std::setw(18) << std::left << input.MuPhiInput << ' '
                << std::setw(18) << std::left << model.get_vs() << ' '
                << std::setw(18) << std::left << model.get_vsb() << ' '
                << std::setw(18) << std::left << model.get_vphi() << ' '
                << std::setw(18) << std::left << higgs << ' '
                << std::setw(18) << std::left << second_higgs << ' '
                << std::setw(18) << std::left << Mhh_DRbar(0) << ' '
                << std::setw(18) << std::left << Mhh_DRbar(1) << ' '
                << std::setw(18) << std::left << Mhh_DRbar(2) << ' '
                << std::setw(18) << std::left << Mhh_DRbar(3) << ' '
                << std::setw(18) << std::left << Mhh_DRbar(4) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(0,0) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(0,1) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(0,2) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(0,3) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(0,4) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(1,1) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(1,2) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(1,3) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(1,4) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(2,2) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(2,3) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(2,4) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(3,3) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(3,4) << ' '
                << std::setw(18) << std::left << mass_matrix_hh(4,4) << ' '
                // << std::setw(18) << std::left << term_1 << ' '
                // << std::setw(18) << std::left << term_2 << ' '
                // << std::setw(18) << std::left << term_3 << ' '
                // << std::setw(18) << std::left << term_4 << ' '
                // << std::setw(18) << std::left << term_5 << ' '
                // << std::setw(18) << std::left << term_6 << ' '
                // << std::setw(18) << std::left << term_7 << ' '
                << std::setw(18) << std::left << ZH(0,0) << ' '
                << std::setw(18) << std::left << ZH(0,1) << ' '
                << std::setw(18) << std::left << ZH(0,2) << ' '
                << std::setw(18) << std::left << ZH(0,3) << ' '
                << std::setw(18) << std::left << ZH(0,4) << ' '
                << std::setw(18) << std::left << ZH(1,0) << ' '
                << std::setw(18) << std::left << ZH(1,1) << ' '
                << std::setw(18) << std::left << ZH(1,2) << ' '
                << std::setw(18) << std::left << ZH(1,3) << ' '
                << std::setw(18) << std::left << ZH(1,4) << ' '
                << std::setw(18) << std::left << ZH(2,0) << ' '
                << std::setw(18) << std::left << ZH(2,1) << ' '
                << std::setw(18) << std::left << ZH(2,2) << ' '
                << std::setw(18) << std::left << ZH(2,3) << ' '
                << std::setw(18) << std::left << ZH(2,4) << ' '
                << std::setw(18) << std::left << ZH(3,0) << ' '
                << std::setw(18) << std::left << ZH(3,1) << ' '
                << std::setw(18) << std::left << ZH(3,2) << ' '
                << std::setw(18) << std::left << ZH(3,3) << ' '
                << std::setw(18) << std::left << ZH(3,4) << ' '
                << std::setw(18) << std::left << ZH(4,0) << ' '
                << std::setw(18) << std::left << ZH(4,1) << ' '
                << std::setw(18) << std::left << ZH(4,2) << ' '
                << std::setw(18) << std::left << ZH(4,3) << ' '
                << std::setw(18) << std::left << ZH(4,4) << ' '
                << std::setw(18) << std::left << error;
      if (error) {
         std::cout << "\t# " << problems;
      }
      std::cout << '\n';
   }

   return 0;
}
