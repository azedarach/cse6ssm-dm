// ====================================================================
// Does a grid scan of the CSE6SSM parameter space using the
// semianalytic algorithm
// ====================================================================

#include "CSE6SSM_semi_two_scale_input_parameters.hpp"
#include "CSE6SSM_semianalytic_spectrum_generator.hpp"
#include "CSE6SSM_scan_utilities.hpp"

#include "error.hpp"
#include "grid_scanner.hpp"
#include "lowe.h"
#include "scan_command_line_options.hpp"
#include "scan_parser.hpp"

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <sys/time.h>

std::size_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

static std::minstd_rand generator(seed);

namespace flexiblesusy {

void set_default_parameter_values(CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   input.m12 = 3500.;
   input.Azero = 8000.;
   input.TanBeta = 10.;
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
   input.LambdaxInput = 7.2e-4;
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

void set_input_parameter_by_name(const std::string& name, double value,
                                 CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   if (name == "m12") {
      input.m12 = value;
   } else if (name == "Azero") {
      input.Azero = value;
   } else if (name == "TanBeta") {
      input.TanBeta = value;
   } else if (name == "sInput") {
      input.sInput = value;
   } else if (name == "QSInput") {
      input.QSInput = value;
   } else if (name == "hEInput(0,0)") {
      input.hEInput(0,0) = value;
   } else if (name == "hEInput(1,0)") {
      input.hEInput(1,0) = value;
   } else if (name == "hEInput(2,0)") {
      input.hEInput(2,0) = value;
   } else if (name == "hEInput(0,1)") {
      input.hEInput(0,1) = value;
   } else if (name == "hEInput(1,1)") {
      input.hEInput(1,1) = value;
   } else if (name == "hEInput(2,1)") {
      input.hEInput(2,1) = value;
   } else if (name == "SigmaLInput") {
      input.SigmaLInput = value;
   } else if (name == "KappaPrInput") {
      input.KappaPrInput = value;
   } else if (name == "SigmaxInput") {
      input.SigmaxInput = value;
   } else if (name == "gDInput(0,0)") {
      input.gDInput(0,0) = value;
   } else if (name == "gDInput(1,0)") {
      input.gDInput(1,0) = value;
   } else if (name == "gDInput(2,0)") {
      input.gDInput(2,0) = value;
   } else if (name == "gDInput(0,1)") {
      input.gDInput(0,1) = value;
   } else if (name == "gDInput(1,1)") {
      input.gDInput(1,1) = value;
   } else if (name == "gDInput(2,1)") {
      input.gDInput(2,1) = value;
   } else if (name == "gDInput(0,2)") {
      input.gDInput(0,2) = value;
   } else if (name == "gDInput(1,2)") {
      input.gDInput(1,2) = value;
   } else if (name == "gDInput(2,2)") {
      input.gDInput(2,2) = value;
   } else if (name == "KappaInput(0,0)") {
      input.KappaInput(0,0) = value;
   } else if (name == "KappaInput(1,0)") {
      input.KappaInput(1,0) = value;
   } else if (name == "KappaInput(2,0)") {
      input.KappaInput(2,0) = value;
   } else if (name == "KappaInput(0,1)") {
      input.KappaInput(0,1) = value;
   } else if (name == "KappaInput(1,1)") {
      input.KappaInput(1,1) = value;
   } else if (name == "KappaInput(2,1)") {
      input.KappaInput(2,1) = value;
   } else if (name == "KappaInput(0,2)") {
      input.KappaInput(0,2) = value;
   } else if (name == "KappaInput(1,2)") {
      input.KappaInput(1,2) = value;
   } else if (name == "KappaInput(2,2)") {
      // NOTE: temporary modification
      input.KappaInput(0,0) = value;
      input.KappaInput(1,1) = value;
      input.KappaInput(2,2) = value;
   }  else if (name == "Lambda12Input(0,0)") {
      input.Lambda12Input(0,0) = value;
   } else if (name == "Lambda12Input(1,0)") {
      input.Lambda12Input(1,0) = value;
   } else if (name == "Lambda12Input(0,1)") {
      input.Lambda12Input(0,1) = value;
   } else if (name == "Lambda12Input(1,1)") {
      // NOTE: temporary modification
      input.Lambda12Input(0,0) = value;
      input.Lambda12Input(1,1) = value;
   } else if (name == "LambdaxInput") {
      input.LambdaxInput = value;
   } else if (name == "fuInput(0,0)") {
      input.fuInput(0,0) = value;
   } else if (name == "fuInput(1,0)") {
      input.fuInput(1,0) = value;
   } else if (name == "fuInput(2,0)") {
      input.fuInput(2,0) = value;
   } else if (name == "fuInput(0,1)") {
      input.fuInput(0,1) = value;
   } else if (name == "fuInput(1,1)") {
      input.fuInput(1,1) = value;
   } else if (name == "fuInput(2,1)") {
      input.fuInput(2,1) = value;
   } else if (name == "fdInput(0,0)") {
      input.fdInput(0,0) = value;
   } else if (name == "fdInput(1,0)") {
      input.fdInput(1,0) = value;
   } else if (name == "fdInput(2,0)") {
      input.fdInput(2,0) = value;
   } else if (name == "fdInput(0,1)") {
      input.fdInput(0,1) = value;
   } else if (name == "fdInput(1,1)") {
      input.fdInput(1,1) = value;
   } else if (name == "fdInput(2,1)") {
      input.fdInput(2,1) = value;
   } else if (name == "MuPrInput") {
      input.MuPrInput = value;
   } else if (name == "MuPhiInput") {
      input.MuPhiInput = value;
   } else if (name == "BMuPrInput") {
      input.BMuPrInput = value;
   } else if (name == "BMuPhiInput") {
      input.BMuPhiInput = value;
   } else {
      WARNING("Unrecognised parameter '" + name + "'!");
   }
}

void set_fixed_parameter_values(const std::vector<Fixed_parameter>& params,
                                CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   for (std::vector<Fixed_parameter>::const_iterator it = params.begin(),
           end = params.end(); it != end; ++it) {
      set_input_parameter_by_name(it->name, it->value, input);
   }
}

void set_gridscan_parameter_values(const std::vector<Scanned_parameter>& params,
                                   const std::vector<std::size_t>& posn,
                                   CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   const std::size_t num_params = params.size();
   for (std::size_t i = 0; i < num_params; ++i) {
      double value = params[i].lower_bound;
      if (params[i].num_points > 1) {
         value += posn[i] * (params[i].upper_bound - params[i].lower_bound) / (params[i].num_points - 1.0);
      }
      set_input_parameter_by_name(params[i].name, value, input);
   }
}

void set_random_parameter_values(const std::vector<Scanned_parameter>& params,
                                 const std::vector<std::size_t>& posn,
                                 CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   std::uniform_real_distribution<double> distribution(0., 1.);
   const std::size_t num_params = params.size();
   for (std::size_t i = 0; i < num_params; ++i) {
      double value = params[i].lower_bound + (params[i].upper_bound - params[i].lower_bound) * distribution(generator);
      set_input_parameter_by_name(params[i].name, value, input);
   }
}

} // namespace flexiblesusy

double get_wall_time()
{
   struct timeval time;
   if (gettimeofday(&time,NULL)) {
      return 0;
   }
   return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}

double get_cpu_time()
{
   return (double)clock() / CLOCKS_PER_SEC;
}

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;
   typedef std::chrono::duration<int,std::micro> microseconds_t;

   std::chrono::high_resolution_clock::time_point start_point = std::chrono::high_resolution_clock::now();
   double wall_start = get_wall_time();
   double cpu_start = get_cpu_time();

   Scan_command_line_options options(argc, argv);
   if (options.must_print_model_info())
      CSE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string scan_input_file(options.get_scan_input_file());
   const std::string pole_mass_output_file(options.get_pole_mass_output_file());
   const std::string drbar_mass_output_file(options.get_drbar_mass_output_file());
   const std::string drbar_susy_pars_output_file(options.get_drbar_susy_pars_output_file());
   const std::string drbar_soft_pars_output_file(options.get_drbar_soft_pars_output_file());
   const std::string drbar_mixings_output_file(options.get_drbar_mixings_output_file());
   const std::string slha_pole_mass_output_file(options.get_slha_pole_mass_output_file());
   const std::string slha_running_mass_output_file(options.get_slha_running_mass_output_file());
   const std::string slha_susy_pars_output_file(options.get_slha_susy_pars_output_file());
   const std::string slha_soft_pars_output_file(options.get_slha_soft_pars_output_file());
   const std::string slha_pole_mixings_output_file(options.get_slha_pole_mixings_output_file());
   const std::string slha_running_mixings_output_file(options.get_slha_running_mixings_output_file());
   const std::string coefficients_output_file(options.get_coefficients_output_file());

   Scan_parser parser;
   if (scan_input_file.empty()) {
      WARNING("No scan input file given!\n"
              "   Default scan parameters will be used.\n"
              "   You can provide them via the option --scan-input-file="); 
   } else {
      try {
         parser.parse_scan_inputs_file(scan_input_file);
      } catch (const ReadError& error) {
         ERROR(error.what());
         return EXIT_FAILURE;
      }
   }

   // output streams
   std::ofstream pole_mass_out_stream;
   if (!pole_mass_output_file.empty()) {
      pole_mass_out_stream.open(pole_mass_output_file, std::ofstream::out);
   }

   std::ostream & pole_mass_out = pole_mass_output_file.empty() ? std::cout : pole_mass_out_stream;

   bool must_write_drbar_masses = false;
   bool must_write_drbar_susy_pars = false;
   bool must_write_drbar_soft_pars = false;
   bool must_write_drbar_mixings = false;
   bool must_write_slha_pole_masses = false;
   bool must_write_slha_running_masses = false;
   bool must_write_slha_susy_pars = false;
   bool must_write_slha_soft_pars = false;
   bool must_write_slha_pole_mixings = false;
   bool must_write_slha_running_mixings = false;
   bool must_write_coefficients = false;

   std::ofstream drbar_mass_out_stream;
   std::ofstream drbar_susy_pars_out_stream;
   std::ofstream drbar_soft_pars_out_stream;
   std::ofstream drbar_mixings_out_stream;
   std::ofstream slha_pole_mass_out_stream;
   std::ofstream slha_running_mass_out_stream;
   std::ofstream slha_susy_pars_out_stream;
   std::ofstream slha_soft_pars_out_stream;
   std::ofstream slha_pole_mixings_out_stream;
   std::ofstream slha_running_mixings_out_stream;
   std::ofstream coefficients_out_stream;
   if (!drbar_mass_output_file.empty()) {
      must_write_drbar_masses = true;
      drbar_mass_out_stream.open(drbar_mass_output_file, std::ofstream::out);
   }
   if (!drbar_susy_pars_output_file.empty()) {
      must_write_drbar_susy_pars = true;
      drbar_susy_pars_out_stream.open(drbar_susy_pars_output_file, std::ofstream::out);
   }
   if (!drbar_soft_pars_output_file.empty()) {
      must_write_drbar_soft_pars = true;
      drbar_soft_pars_out_stream.open(drbar_soft_pars_output_file, std::ofstream::out);
   }
   if (!drbar_mixings_output_file.empty()) {
      must_write_drbar_mixings = true;
      drbar_mixings_out_stream.open(drbar_mixings_output_file, std::ofstream::out);
   }
   if (!slha_pole_mass_output_file.empty()) {
      must_write_slha_pole_masses = true;
      slha_pole_mass_out_stream.open(slha_pole_mass_output_file, std::ofstream::out);
   }
   if (!slha_running_mass_output_file.empty()) {
      must_write_slha_running_masses = true;
      slha_running_mass_out_stream.open(slha_running_mass_output_file, std::ofstream::out);
   }
   if (!slha_susy_pars_output_file.empty()) {
      must_write_slha_susy_pars = true;
      slha_susy_pars_out_stream.open(slha_susy_pars_output_file, std::ofstream::out);
   }
   if (!slha_soft_pars_output_file.empty()) {
      must_write_slha_soft_pars = true;
      slha_soft_pars_out_stream.open(slha_soft_pars_output_file, std::ofstream::out);
   }
   if (!slha_pole_mixings_output_file.empty()) {
      must_write_slha_pole_mixings = true;
      slha_pole_mixings_out_stream.open(slha_pole_mixings_output_file, std::ofstream::out);
   }
   if (!slha_running_mixings_output_file.empty()) {
      must_write_slha_running_mixings = true;
      slha_running_mixings_out_stream.open(slha_running_mixings_output_file, std::ofstream::out);
   }
   if (!coefficients_output_file.empty()) {
      must_write_coefficients = true;
      coefficients_out_stream.open(coefficients_output_file, std::ofstream::out);
   }

   CSE6SSM_semianalytic_input_parameters<algorithm_type> input;
   set_default_parameter_values(input);
   set_fixed_parameter_values(parser.get_fixed_parameters(), input);

   QedQcd oneset;
   oneset.toMz();

   std::vector<std::size_t> scan_dimensions;

   std::vector<Scanned_parameter> scanned_parameters
      = parser.get_scanned_parameters();

   const bool is_grid_scan = parser.is_grid_scan();
   if (is_grid_scan) {
      for (std::size_t i = 0; i < scanned_parameters.size(); ++i) {
         scan_dimensions.push_back(scanned_parameters[i].num_points);
      }
   } else {
      for (std::size_t i = 0; i < scanned_parameters.size(); ++i) {
         scan_dimensions.push_back(1);
      }
      scan_dimensions.back() = parser.get_number_of_points();
   }

   Grid_scanner scan(scan_dimensions);

   CSE6SSM_semianalytic_pole_mass_writer pole_mass_writer;
   CSE6SSM_semianalytic_drbar_values_writer drbar_values_writer;
   CSE6SSM_semianalytic_slha_values_writer slha_values_writer;
   CSE6SSM_semianalytic_coefficients_writer coefficients_writer;
   bool must_write_comment_line = true;

   while (!scan.has_finished()) {
      if (is_grid_scan) {
         set_gridscan_parameter_values(scanned_parameters, scan.get_position(), input);
      } else {
         set_random_parameter_values(scanned_parameters, scan.get_position(), input);
      }
      CSE6SSM_semianalytic_spectrum_generator<algorithm_type>
         spectrum_generator;
      spectrum_generator.set_precision_goal(1.0e-3);
      spectrum_generator.set_max_iterations(0);   // 0 == automatic
      spectrum_generator.set_calculate_sm_masses(1); // 0 == no
      spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale
      spectrum_generator.set_ewsb_loop_order(2);
      spectrum_generator.set_pole_mass_loop_order(2);
      spectrum_generator.set_beta_loop_order(2);
      spectrum_generator.set_threshold_corrections_loop_order(1);
      spectrum_generator.set_force_output(1);

      // note
      spectrum_generator.set_ewsb_iteration_precision(100.0);

      spectrum_generator.run(oneset, input);

      const CSE6SSM_semianalytic<algorithm_type>& model
         = spectrum_generator.get_model();
      const CSE6SSM_semianalytic_slha<algorithm_type> model_slha(model);

      slha_values_writer.set_high_scale(spectrum_generator.get_high_scale());
      slha_values_writer.set_susy_scale(spectrum_generator.get_susy_scale());
      slha_values_writer.set_low_scale(spectrum_generator.get_low_scale());

      coefficients_writer.set_high_scale(spectrum_generator.get_high_scale());
      coefficients_writer.set_susy_scale(spectrum_generator.get_susy_scale());
      coefficients_writer.set_low_scale(spectrum_generator.get_low_scale());

      pole_mass_writer.extract_pole_masses(model);

      if (must_write_drbar_masses)
         drbar_values_writer.extract_drbar_masses(model);
      if (must_write_drbar_susy_pars)
         drbar_values_writer.extract_drbar_susy_pars(model);
      if (must_write_drbar_soft_pars)
         drbar_values_writer.extract_drbar_soft_pars(model);
      if (must_write_drbar_mixings)
         drbar_values_writer.extract_drbar_mixings(model);
      if (must_write_slha_pole_masses)
         slha_values_writer.extract_slha_pole_masses(model_slha);
      if (must_write_slha_running_masses)
         slha_values_writer.extract_slha_running_masses(model_slha);
      if (must_write_slha_susy_pars)
         slha_values_writer.extract_slha_susy_pars(model_slha);
      if (must_write_slha_soft_pars)
         slha_values_writer.extract_slha_soft_pars(model_slha);
      if (must_write_slha_pole_mixings)
         slha_values_writer.extract_slha_pole_mixings(model_slha);
      if (must_write_slha_running_mixings)
         slha_values_writer.extract_slha_running_mixings(model_slha);
      if (must_write_coefficients)
         coefficients_writer.extract_coefficients(model);

      if (must_write_comment_line) {
         pole_mass_writer.write_pole_masses_comment_line(pole_mass_out);

         if (must_write_drbar_masses)
            drbar_values_writer.write_drbar_masses_comment_line(drbar_mass_out_stream);
         if (must_write_drbar_susy_pars)
            drbar_values_writer.write_drbar_susy_pars_comment_line(drbar_susy_pars_out_stream);
         if (must_write_drbar_soft_pars)
            drbar_values_writer.write_drbar_soft_pars_comment_line(drbar_soft_pars_out_stream);
         if (must_write_drbar_mixings)
            drbar_values_writer.write_drbar_mixings_comment_line(drbar_mixings_out_stream);
         if (must_write_slha_pole_masses)
            slha_values_writer.write_slha_pole_masses_comment_line(slha_pole_mass_out_stream);
         if (must_write_slha_running_masses)
            slha_values_writer.write_slha_running_masses_comment_line(slha_running_mass_out_stream);
         if (must_write_slha_susy_pars)
            slha_values_writer.write_slha_susy_pars_comment_line(slha_susy_pars_out_stream);
         if (must_write_slha_soft_pars)
            slha_values_writer.write_slha_soft_pars_comment_line(slha_soft_pars_out_stream);
         if (must_write_slha_pole_mixings)
            slha_values_writer.write_slha_pole_mixings_comment_line(slha_pole_mixings_out_stream);
         if (must_write_slha_running_mixings)
            slha_values_writer.write_slha_running_mixings_comment_line(slha_running_mixings_out_stream);
         if (must_write_coefficients)
            coefficients_writer.write_coefficients_comment_line(coefficients_out_stream);

         must_write_comment_line = false;
      }
      pole_mass_writer.write_pole_masses_line(pole_mass_out);

      if (must_write_drbar_masses)
         drbar_values_writer.write_drbar_masses_line(drbar_mass_out_stream);
      if (must_write_drbar_susy_pars)
         drbar_values_writer.write_drbar_susy_pars_line(drbar_susy_pars_out_stream);
      if (must_write_drbar_soft_pars)
         drbar_values_writer.write_drbar_soft_pars_line(drbar_soft_pars_out_stream);
      if (must_write_drbar_mixings)
         drbar_values_writer.write_drbar_mixings_line(drbar_mixings_out_stream);
      if (must_write_slha_pole_masses)
         slha_values_writer.write_slha_pole_masses_line(slha_pole_mass_out_stream);
      if (must_write_slha_running_masses)
         slha_values_writer.write_slha_running_masses_line(slha_running_mass_out_stream);
      if (must_write_slha_susy_pars)
         slha_values_writer.write_slha_susy_pars_line(slha_susy_pars_out_stream);
      if (must_write_slha_soft_pars)
         slha_values_writer.write_slha_soft_pars_line(slha_soft_pars_out_stream);
      if (must_write_slha_pole_mixings)
         slha_values_writer.write_slha_pole_mixings_line(slha_pole_mixings_out_stream);
      if (must_write_slha_running_mixings)
         slha_values_writer.write_slha_running_mixings_line(slha_running_mixings_out_stream);
      if (must_write_coefficients)
         coefficients_writer.write_coefficients_line(coefficients_out_stream);

      scan.step_forward();
   }

   if (pole_mass_out_stream.is_open())
      pole_mass_out_stream.close();

   if (drbar_mass_out_stream.is_open())
      drbar_mass_out_stream.close();

   if (drbar_susy_pars_out_stream.is_open())
      drbar_susy_pars_out_stream.close();

   if (drbar_soft_pars_out_stream.is_open())
      drbar_soft_pars_out_stream.close();

   if (drbar_mixings_out_stream.is_open())
      drbar_mixings_out_stream.close();

   if (slha_pole_mass_out_stream.is_open())
      slha_pole_mass_out_stream.close();

   if (slha_running_mass_out_stream.is_open())
      slha_running_mass_out_stream.close();

   if (slha_susy_pars_out_stream.is_open())
      slha_susy_pars_out_stream.close();

   if (slha_soft_pars_out_stream.is_open())
      slha_soft_pars_out_stream.close();

   if (slha_pole_mixings_out_stream.is_open())
      slha_pole_mixings_out_stream.close();

   if (slha_running_mixings_out_stream.is_open())
      slha_running_mixings_out_stream.close();

   if (coefficients_out_stream.is_open())
      coefficients_out_stream.close();

   std::chrono::high_resolution_clock::time_point end_point = std::chrono::high_resolution_clock::now();
   microseconds_t duration(std::chrono::duration_cast<microseconds_t>(end_point - start_point));
   double time_in_seconds = duration.count() * 0.000001;
   double wall_end = get_wall_time();
   double cpu_end = get_cpu_time();

   cout << "# Scan completed in " << time_in_seconds << " seconds\n";
   cout << "# Wall time = " << wall_end - wall_start << " seconds\n";
   cout << "# CPU time  = " << cpu_end - cpu_start << " seconds\n";
   cout << "# Random seed = " << seed << "\n";

   return 0;
}
