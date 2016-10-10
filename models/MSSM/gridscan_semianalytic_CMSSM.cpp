
#include "CMSSM_semi_two_scale_input_parameters.hpp"
#include "CMSSM_semianalytic_spectrum_generator.hpp"
#include "CMSSM_scan_utilities.hpp"

#include "error.hpp"
#include "grid_scanner.hpp"
#include "lowe.h"
#include "scan_command_line_options.hpp"
#include "scan_parser.hpp"

#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/time.h>

std::size_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

static std::minstd_rand generator(seed);

namespace flexiblesusy {

void set_default_parameter_values(CMSSM_semianalytic_input_parameters<Two_scale>& input)
{
   input.m12 = 3500.;
   input.Azero = 8000.;
   input.TanBeta = 10.;
   input.MuInput = 1000.;
   input.MuInput_at_MS = false;
}

void set_input_parameter_by_name(const std::string& name, double value,
                                 CMSSM_semianalytic_input_parameters<Two_scale>& input)
{
   if (name == "m12") {
      input.m12 = value;
   } else if (name == "Azero") {
      input.Azero = value;
   } else if (name == "TanBeta") {
      input.TanBeta = value;
   } else if (name == "MuInput") {
      input.MuInput = value;
   } else if (name == "MuInput_at_MS") {
      input.MuInput_at_MS = value;
   } else {
      WARNING("Unrecognized parameter '" + name + "'!");
   }
}

void set_fixed_parameter_values(const std::vector<Fixed_parameter>& params,
                                CMSSM_semianalytic_input_parameters<Two_scale>& input)
{
   for (std::vector<Fixed_parameter>::const_iterator it = params.begin(),
           end = params.end(); it != end; ++it) {
      set_input_parameter_by_name(it->name, it->value, input);
   }
}

void set_gridscan_parameter_values(const std::vector<Scanned_parameter>& params,
                                   const std::vector<std::size_t>& posn,
                                   CMSSM_semianalytic_input_parameters<Two_scale>& input)
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
                                 CMSSM_semianalytic_input_parameters<Two_scale>& input)
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
      MSSM_info::CMSSM_info::print(std::cout);
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

   CMSSM_semianalytic_input_parameters<algorithm_type> input;
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

   CMSSM_semianalytic_pole_mass_writer pole_mass_writer;
   CMSSM_semianalytic_drbar_values_writer drbar_values_writer;
   CMSSM_semianalytic_slha_values_writer slha_values_writer;
   bool must_write_comment_line = true;

   while (!scan.has_finished()) {
      if (is_grid_scan) {
         set_gridscan_parameter_values(scanned_parameters, scan.get_position(), input);
      } else {
         set_random_parameter_values(scanned_parameters, scan.get_position(), input);
      }
      CMSSM_semianalytic_spectrum_generator<algorithm_type>
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

      spectrum_generator.run(oneset, input);

      const CMSSM_semianalytic<algorithm_type>& model
         = spectrum_generator.get_model();
      const CMSSM_semianalytic_slha<algorithm_type> model_slha(model);

      slha_values_writer.set_high_scale(spectrum_generator.get_high_scale());
      slha_values_writer.set_susy_scale(spectrum_generator.get_susy_scale());
      slha_values_writer.set_low_scale(spectrum_generator.get_low_scale());

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
