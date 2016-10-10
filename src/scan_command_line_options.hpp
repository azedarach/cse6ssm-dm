// ====================================================================
// Class providing command line options for grid scan.
// - better way to do this would be to inherit from command_line_options
// ====================================================================

#ifndef SCAN_COMMAND_LINE_OPTIONS_H
#define SCAN_COMMAND_LINE_OPTIONS_H

#include <iosfwd>
#include <string>
#include "error.hpp"

namespace flexiblesusy {

/**
 * @class Scan_command_line_options
 * @brief parses the command line options
 *
 * Usage:
 * @code
 * int main(int argc, const char* argv[])
 * {
 *    using namespace flexiblesusy;
 *    Scan_command_line_options options(argc, argv);
 *    if (options.must_exit())
 *       return options.status();
 *    cout << "Scan input file: " << options.get_scan_input_file() << endl;
 * }
 * @endcode
 */
class Scan_command_line_options {
public:
   Scan_command_line_options();
   Scan_command_line_options(int, const char*[]);
   ~Scan_command_line_options();

   bool must_exit() const { return do_exit; }
   bool must_print_model_info() const { return do_print_model_info; }
   int status() const { return exit_status; }
   void parse(int, const char*[]);
   void print_build_info(std::ostream&) const;
   void print_usage(std::ostream&) const;
   void print_version(std::ostream&) const;
   void reset();

   const std::string& get_scan_input_file() const { return scan_input_file; }
   const std::string& get_pole_mass_output_file() const { return pole_mass_output_file; }
   const std::string& get_drbar_mass_output_file() const { return drbar_mass_output_file; }
   const std::string& get_drbar_susy_pars_output_file() const { return drbar_susy_pars_output_file; }
   const std::string& get_drbar_soft_pars_output_file() const { return drbar_soft_pars_output_file; }
   const std::string& get_drbar_mixings_output_file() const { return drbar_mixings_output_file; }
   const std::string& get_slha_pole_mass_output_file() const { return slha_pole_mass_output_file; }
   const std::string& get_slha_running_mass_output_file() const { return slha_running_mass_output_file; }
   const std::string& get_slha_susy_pars_output_file() const { return slha_susy_pars_output_file; }
   const std::string& get_slha_soft_pars_output_file() const { return slha_soft_pars_output_file; }
   const std::string& get_slha_pole_mixings_output_file() const { return slha_pole_mixings_output_file; }
   const std::string& get_slha_running_mixings_output_file() const { return slha_running_mixings_output_file; }
   const std::string& get_coefficients_output_file() const { return coefficients_output_file; }
   const std::string& get_program_name() const { return program; }

   static std::string get_option_value(const std::string&, const std::string&);
   static bool get_parameter_value(const std::string&, const std::string&, double&);
   static bool get_parameter_value(const std::string&, const std::string&, int&);
   static bool starts_with(const std::string&, const std::string&);

private:
   bool do_exit;
   bool do_print_model_info;
   int exit_status;
   std::string program;
   std::string scan_input_file;
   std::string pole_mass_output_file;
   std::string drbar_mass_output_file;
   std::string drbar_susy_pars_output_file;
   std::string drbar_soft_pars_output_file;
   std::string drbar_mixings_output_file;
   std::string slha_pole_mass_output_file;
   std::string slha_running_mass_output_file;
   std::string slha_susy_pars_output_file;
   std::string slha_soft_pars_output_file;
   std::string slha_pole_mixings_output_file;
   std::string slha_running_mixings_output_file;
   std::string coefficients_output_file;
};

} // namespace flexiblesusy

#endif
