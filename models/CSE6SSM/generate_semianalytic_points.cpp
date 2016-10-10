// ====================================================================
// Generates points from an input file, but does not run them in the
// spectrum generator
// ====================================================================

#include "CSE6SSM_semi_two_scale_input_parameters.hpp"
#include "CSE6SSM_scan_utilities.hpp"

#include "command_line_options.hpp"
#include "error.hpp"
#include "grid_scanner.hpp"
#include "scan_parser.hpp"

#include <iostream>
#include <fstream>
#include <chrono>

std::size_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

static std::minstd_rand generator(seed);

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: ./generate_semianalytic_points.x [options]\n"
      "Options:\n"
      "  --scan-input-file=<filename>              input file defining scan\n"
      "  --points-output-file=<filename>           file to output points to\n"
      "  --help,-h                                 print this help message"
             << std::endl;
}

void parse_command_line_options(int argc, const char* argv[], 
                                std::string& scan_input_file,
                                std::string& points_output_file)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if (Command_line_options::starts_with(option,"--scan-input-file")) {
         scan_input_file
            = Command_line_options::get_option_value(option,"=");
         continue;
      }

      if (Command_line_options::starts_with(option,"--points-output-file")) {
         points_output_file
            = Command_line_options::get_option_value(option,"=");
         continue;
      }

      if (strcmp(option,"--help") == 0 || strcmp(option,"-h") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}


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
      input.KappaInput(2,2) = value;
   }  else if (name == "Lambda12Input(0,0)") {
      input.Lambda12Input(0,0) = value;
   } else if (name == "Lambda12Input(1,0)") {
      input.Lambda12Input(1,0) = value;
   } else if (name == "Lambda12Input(0,1)") {
      input.Lambda12Input(0,1) = value;
   } else if (name == "Lambda12Input(1,1)") {
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

void write_comment_line(std::ostream & out)
{
   const std::size_t width = 18;

   out << "# ";
   write_CSE6SSM_semianalytic_inputs_list(out, width);
   out << '\n';
}

void write_point(const CSE6SSM_semianalytic_input_parameters<Two_scale>& input,
                 std::ostream & out)
{
   const std::size_t width = 18;

   write_CSE6SSM_inputs(input, out, width);
   out << '\n';
}

} // namespace flexiblesusy

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   typedef Two_scale algorithm_type;

   std::string scan_input_file;
   std::string points_output_file;

   parse_command_line_options(argc, argv, scan_input_file, points_output_file);

   if (scan_input_file.empty()) {
      ERROR("No scan input file given!\n"
            "   You can provide one via the option --scan-input-file=");
      return EXIT_FAILURE;
   }

   Scan_parser parser;
   try {
      parser.parse_scan_inputs_file(scan_input_file);
   } catch (const ReadError& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   std::ofstream points_out_stream;
   if (!points_output_file.empty()) {
      points_out_stream.open(points_output_file, std::ofstream::out);
   }

   std::ostream & points_out = points_output_file.empty() ? std::cout
      : points_out_stream;

   CSE6SSM_semianalytic_input_parameters<algorithm_type> input;
   set_default_parameter_values(input);
   set_fixed_parameter_values(parser.get_fixed_parameters(), input);

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

   bool must_write_comment_line = true;

   while (!scan.has_finished()) {
      if (is_grid_scan) {
         set_gridscan_parameter_values(scanned_parameters, scan.get_position(), input);
      } else {
         set_random_parameter_values(scanned_parameters, scan.get_position(), input);
      }

      if (must_write_comment_line) {
         write_comment_line(points_out);
         must_write_comment_line = false;
      }

      write_point(input, points_out);

      scan.step_forward();
   }

   if (points_out_stream.is_open())
      points_out_stream.close();

   return 0;
}
