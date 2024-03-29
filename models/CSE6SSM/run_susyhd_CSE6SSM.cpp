// ====================================================================
// Program to calculate pole masses and mixings of a CSE6SSM parameter
// space point using a semianalytic solver
// ====================================================================

#include "CSE6SSM_semi_two_scale_input_parameters.hpp"
#include "CSE6SSM_slha_io.hpp"
#include "CSE6SSM_semianalytic_spectrum_generator.hpp"

#include "susyhd_call.hpp"

#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

SUSYHD::SUSYHD_input_parameters match_to_MSSM(flexiblesusy::CSE6SSM_semianalytic<flexiblesusy::Two_scale> model)
{
   SUSYHD::SUSYHD_input_parameters input;

   input.TanBeta = model.get_vu() / model.get_vd();
   input.M1 = model.get_MassB();
   input.M2 = model.get_MassWB();
   input.M3 = model.get_MassG();
   input.Mu = model.get_Lambdax() * model.get_vs() / flexiblesusy::Sqrt(2.0);
   input.At = model.get_TYu(2,2) / model.get_Yu(2,2);
   input.mQ3 = flexiblesusy::Sqrt(model.get_mq2(2,2));
   input.mU3 = flexiblesusy::Sqrt(model.get_mu2(2,2));
   input.mD3 = flexiblesusy::Sqrt(model.get_md2(2,2));
   input.mQ2 = flexiblesusy::Sqrt(model.get_mq2(1,1));
   input.mU2 = flexiblesusy::Sqrt(model.get_mu2(1,1));
   input.mD2 = flexiblesusy::Sqrt(model.get_md2(1,1));
   input.mQ1 = flexiblesusy::Sqrt(model.get_mq2(0,0));
   input.mU1 = flexiblesusy::Sqrt(model.get_mu2(0,0));
   input.mD1 = flexiblesusy::Sqrt(model.get_md2(0,0));
   input.mL3 = flexiblesusy::Sqrt(model.get_ml2(2,2));
   input.mE3 = flexiblesusy::Sqrt(model.get_me2(2,2));
   input.mL2 = flexiblesusy::Sqrt(model.get_ml2(1,1));
   input.mE2 = flexiblesusy::Sqrt(model.get_me2(1,1));
   input.mL1 = flexiblesusy::Sqrt(model.get_ml2(0,0));
   input.mE1 = flexiblesusy::Sqrt(model.get_me2(0,0));
   input.mA = model.get_MAh(2);

   return input;
}

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      CSE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   CSE6SSM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;
   QedQcd oneset;
   CSE6SSM_semianalytic_input_parameters<algorithm_type> input;

   if (slha_input_file.empty()) {
      ERROR("No SLHA input file given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_file(slha_input_file);
      slha_io.fill(oneset);
      slha_io.fill(input);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   oneset.toMz(); // run SM fermion masses to MZ

   CSE6SSM_semianalytic_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   spectrum_generator.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_force_output(
      spectrum_generator_settings.get(Spectrum_generator_settings::force_output) >= 1.0);
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order));
   spectrum_generator.set_beta_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::beta_loop_order));
   spectrum_generator.set_threshold_corrections_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::threshold_corrections_loop_order));
   spectrum_generator.set_two_loop_corrections(
      spectrum_generator_settings.get_two_loop_corrections());
   // note
   spectrum_generator.set_ewsb_iteration_precision(0.5);

   spectrum_generator.run(oneset, input);

   const CSE6SSM_semianalytic_slha<algorithm_type> model(spectrum_generator.get_model());
   const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   CSE6SSM_scales scales;
   scales.HighScale = spectrum_generator.get_high_scale();
   scales.SUSYScale = spectrum_generator.get_susy_scale();
   scales.LowScale = spectrum_generator.get_low_scale();

   // SUSYHD calculation of Higgs mass
   double eft_mhh = 0.;
   try {
      mathematica::MathematicaLink link("-linkname \"math -mathlink\"");
      SUSYHD::SUSYHDLink susyhd(&link);

      eft_mhh = susyhd.calculate_MHiggs(match_to_MSSM(model));
   } catch (const mathematica::Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   // output
   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(oneset);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!problems.have_problem() ||
       spectrum_generator_settings.get(Spectrum_generator_settings::force_output)) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model, scales, eft_mhh);
   }

   if (slha_output_file.empty()) {
      slha_io.write_to_stream(std::cout);
   } else {
      slha_io.write_to_file(slha_output_file);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;
}
