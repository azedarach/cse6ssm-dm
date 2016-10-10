#include "scan_command_line_options.hpp"
#include "build_info.hpp"
#include "logger.hpp"

#include <cstdio>
#include <cstdlib>

namespace flexiblesusy {

Scan_command_line_options::Scan_command_line_options()
   : do_exit(false)
   , do_print_model_info(false)
   , exit_status(EXIT_SUCCESS)
   , scan_input_file()
   , pole_mass_output_file()
   , drbar_mass_output_file()
   , drbar_susy_pars_output_file()
   , drbar_soft_pars_output_file()
   , drbar_mixings_output_file()
   , slha_pole_mass_output_file()
   , slha_running_mass_output_file()
   , slha_susy_pars_output_file()
   , slha_soft_pars_output_file()
   , slha_pole_mixings_output_file()
   , slha_running_mixings_output_file()
   , coefficients_output_file()
{
}

Scan_command_line_options::Scan_command_line_options(int argc, const char* argv[])
   : do_exit(false)
   , do_print_model_info(false)
   , exit_status(EXIT_SUCCESS)
   , program()
   , scan_input_file()
{
   parse(argc, argv);
}

Scan_command_line_options::~Scan_command_line_options()
{
}

/**
 * Parse string of program options and store the given option values
 * in the member variables of this class.
 *
 * @param argc number of program arguments
 * @param argv program arguments
 */
void Scan_command_line_options::parse(int argc, const char* argv[])
{
   assert(argc > 0);
   reset();
   program = argv[0];

   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);
      if (starts_with(option,"--scan-input-file=")) {
         scan_input_file = get_option_value(option, "=");
         if (scan_input_file.empty())
            WARNING("no scan input file name given");
      } else if (starts_with(option,"--pole-mass-output-file=")) {
         pole_mass_output_file = get_option_value(option, "=");
         if (pole_mass_output_file.empty())
            WARNING("no pole mass output file name given");
      } else if (starts_with(option,"--drbar-mass-output-file=")) {
         drbar_mass_output_file = get_option_value(option, "=");
         if (drbar_mass_output_file.empty())
            WARNING("no DRbar mass output file name given");
      } else if (starts_with(option,"--drbar-susy-pars-output-file=")) {
         drbar_susy_pars_output_file = get_option_value(option, "=");
         if (drbar_susy_pars_output_file.empty())
            WARNING("no DRbar SUSY parameters output file name given");
      } else if (starts_with(option,"--drbar-soft-pars-output-file=")) {
         drbar_soft_pars_output_file = get_option_value(option, "=");
         if (drbar_soft_pars_output_file.empty())
            WARNING("no DRbar soft parameters output file name given");
      } else if (starts_with(option,"--drbar-mixings-output-file=")) {
         drbar_mixings_output_file = get_option_value(option, "=");
         if (drbar_mixings_output_file.empty())
            WARNING("no DRbar mixings output file name given");
      } else if (starts_with(option,"--slha-pole-mass-output-file=")) {
         slha_pole_mass_output_file = get_option_value(option, "=");
         if (slha_pole_mass_output_file.empty())
            WARNING("no SLHA pole mass output file name given");
      } else if (starts_with(option,"--slha-running-mass-output-file=")) {
         slha_running_mass_output_file = get_option_value(option, "=");
         if (slha_running_mass_output_file.empty())
            WARNING("no SLHA running mass output file name given");
      } else if (starts_with(option,"--slha-susy-pars-output-file=")) {
         slha_susy_pars_output_file = get_option_value(option, "=");
         if (slha_susy_pars_output_file.empty())
            WARNING("no SLHA SUSY parameters output file name given");
      } else if (starts_with(option,"--slha-soft-pars-output-file=")) {
         slha_soft_pars_output_file = get_option_value(option, "=");
         if (slha_soft_pars_output_file.empty())
            WARNING("no SLHA soft parameters output file name given");
      } else if (starts_with(option,"--slha-pole-mixings-output-file=")) {
         slha_pole_mixings_output_file = get_option_value(option, "=");
         if (slha_pole_mixings_output_file.empty())
            WARNING("no SLHA pole mixings output file name given");
      } else if (starts_with(option,"--slha-running-mixings-output-file=")) {
         slha_running_mixings_output_file = get_option_value(option, "=");
         if (slha_running_mixings_output_file.empty())
            WARNING("no SLHA running mixings output file name given");
      } else if (starts_with(option,"--coefficients-output-file=")) {
         coefficients_output_file = get_option_value(option, "=");
         if (coefficients_output_file.empty())
            WARNING("no coefficients output file name given");
      } else if (option == "--help" || option == "-h") {
         print_usage(std::cout);
         do_exit = true;
      } else if (option == "--build-info") {
         print_build_info(std::cout);
         do_exit = true;
      } else if (option == "--model-info") {
         do_print_model_info = true;
         do_exit = true;
      } else if (option == "--version" || option == "-v") {
         print_version(std::cout);
         do_exit = true;
      } else {
         ERROR("Unrecognized command line option: " << option);
         do_exit = true;
         exit_status = EXIT_FAILURE;
      }
   }
}

void Scan_command_line_options::print_build_info(std::ostream& ostr) const
{
   flexiblesusy::print_all_info(ostr);
}

void Scan_command_line_options::print_version(std::ostream& ostr) const
{
   ostr << "FlexibleSUSY ";
   print_flexiblesusy_version(ostr);
   ostr << std::endl;
}

void Scan_command_line_options::print_usage(std::ostream& ostr) const
{
   ostr << "Usage: " << program << " [options]\n"
           "Options:\n"
           "  --scan-input-file=<filename>                 scan input file\n"
           "                                               If not given, default values\n"
           "                                               are used.\n"
           "  --pole-mass-output-file=<filename>           SLHA output file\n"
           "                                               If not given, the output is\n"
           "                                               printed to stdout.\n"
           "  --drbar-mass-output-file=<filename>          DRbar masses output file\n"
           "                                               If not given, DR bar masses\n"
           "                                               are not printed.\n"
           "  --drbar-susy-pars-output-file=<filename>     DRbar SUSY parameters output file\n"
           "                                               If not given, SUSY parameters\n"
           "                                               are not printed.\n"
           "  --drbar-soft-pars-output-file=<filename>     DRbar soft parameters output file\n"
           "                                               If not given, soft parameters\n"
           "                                               are not printed.\n"
           "  --drbar-mixings-output-file=<filename>       DRbar mixings output file\n"
           "                                               If not given, DR mixings\n"
           "                                               are not printed.\n"
           " --slha-pole-mass-output-file=<filename>       SLHA pole masses output file\n"
           "                                               If not given, SLHA pole masses\n"
           "                                               are not printed.\n"
           " --slha-running-mass-output-file=<filename>    SLHA running masses output file\n"
           "                                               If not given, SLHA running masses\n"
           "                                               are not printed.\n"
           " --slha-susy-pars-output-file=<filename>       SLHA SUSY parameters output file\n"
           "                                               If not given, SLHA SUSY parameters\n"
           "                                               are not printed.\n"
           " --slha-soft-pars-output-file=<filename>       SLHA soft parameters output file\n"
           "                                               If not given, SLHA soft parameters\n"
           "                                               are not printed.\n"
           " --slha-pole-mixings-output-file=<filename>    SLHA pole mixings output file\n"
           "                                               If not given, SLHA pole mixings\n"
           "                                               are not printed.\n"
           " --slha-running-mixings-output-file=<filename> SLHA running mixings output file\n"
           "                                               If not given, SLHA running mixings\n"
           "                                               are not printed.\n"
           " --coefficients-output-file=<filename>         semi-analytic coefficients output\n"
           "                                               file. If not given, coefficients\n"
           "                                               are not printed.\n"
           "  --build-info                                 print build information\n"
           "  --model-info                                 print model information\n"
           "  --help,-h                                    print this help message\n"
           "  --version,-v                                 print program version"
        << std::endl;
}

/**
 * Resets all command line options to their initial values.
 */
void Scan_command_line_options::reset()
{
   do_exit = false;
   do_print_model_info = false;
   exit_status = EXIT_SUCCESS;
   program.clear();
   scan_input_file.clear();
   pole_mass_output_file.clear();
   drbar_mass_output_file.clear();
   drbar_susy_pars_output_file.clear();
   drbar_soft_pars_output_file.clear();
   drbar_mixings_output_file.clear();
   slha_pole_mass_output_file.clear();
   slha_running_mass_output_file.clear();
   slha_susy_pars_output_file.clear();
   slha_soft_pars_output_file.clear();
   slha_pole_mixings_output_file.clear();
   slha_running_mixings_output_file.clear();
   coefficients_output_file.clear();
}

/**
 * Returns true if the string str starts with prefix, false otherwise.
 *
 * @param str string to search in
 * @param prefix string to search for
 *
 * @return true if the string str starts with prefix, false otherwise
 */
bool Scan_command_line_options::starts_with(const std::string& str,
                                       const std::string& prefix)
{
   return !str.compare(0, prefix.size(), prefix);
}

/**
 * Returns a string containing the option value, for options
 * of the form '--option<delim><value>
 *
 * @param str string to search in
 * @param delim character to use a delimiter
 *
 * @return string containing substring starting from delimiter (empty if not found)
 */
   std::string Scan_command_line_options::get_option_value(const std::string& str, const std::string& delim)
{
   std::string value_str("");
   std::size_t posn = str.find(delim);
   if (posn != std::string::npos) {
      value_str = str.substr(++posn);
   } 
   return value_str;
}

/**
 * Extracts the parameter value from a command line option string of
 * the form --m0=125 .
 *
 * @param str full option string, including the parameter value (--m0=125)
 * @param prefix option string, without the parameter value (--m0=)
 * @param[out] parameter output parameter value (of type double)
 *
 * @return true, if str starts with prefix, false otherwise
 */
bool Scan_command_line_options::get_parameter_value(const std::string& str,
                                               const std::string& prefix,
                                               double& parameter)
{
   if (starts_with(str, prefix)) {
      parameter = atof(str.substr(prefix.length()).c_str());
      return true;
   }
   return false;
}

/**
 * Extracts the parameter value from a command line option string of
 * the form --m0=125 .
 *
 * @param str full option string, including the parameter value (--m0=125)
 * @param prefix option string, without the parameter value (--m0=)
 * @param[out] parameter output parameter value (of type int)
 *
 * @return true, if str starts with prefix, false otherwise
 */
bool Scan_command_line_options::get_parameter_value(const std::string& str,
                                               const std::string& prefix,
                                               int& parameter)
{
   if (starts_with(str, prefix)) {
      parameter = atoi(str.substr(prefix.length()).c_str());
      return true;
   }
   return false;
}

} // namespace flexiblesusy
