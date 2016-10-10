// ====================================================================
// Does a grid scan of the CSE6SSM parameter space
//   - trial version before generalising it to arbitrary input
//     parameters (by which I mean arbitrary scans over the 
//     parameters defined in CSE6SSM_input_parameters)
// ====================================================================

#include "CSE6SSM_two_scale_input_parameters.hpp"
#include "CSE6SSM_scan_parameters.hpp"
#include "CSE6SSM_scan_utilities.hpp"
#include "CSE6SSM_spectrum_generator.hpp"

#include "scan_command_line_options.hpp"
#include "error.hpp"
#include "grid_scanner.hpp"
#include "lowe.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/time.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

std::size_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

static std::minstd_rand generator(seed);

namespace flexiblesusy {
   
   void set_default_parameter_values(CSE6SSM_input_parameters<Two_scale>& input)
   {
      if (is_zero(input.TanBeta))
         input.TanBeta = 10.0;

      if (is_zero(input.SignLambdax))
         input.SignLambdax = 1;

      input.sInput = 40000.0; // GeV
      input.QSInput = 5.;
      
      input.hEInput(0,0) = 0.;
      input.hEInput(0,1) = 0.;
      input.hEInput(1,0) = 0.;
      input.hEInput(1,1) = 0.;
      input.hEInput(2,0) = 0.;
      input.hEInput(2,1) = 0.;
      
      input.SigmaLInput = 3.0e-1;
      input.KappaPrInput = 2.0e-2;
      input.SigmaxInput = 1.0e-1;
      
      input.gDInput(0,0) = 0.;
      input.gDInput(0,1) = 0.;
      input.gDInput(0,2) = 0.;
      input.gDInput(1,0) = 0.;
      input.gDInput(1,1) = 0.;
      input.gDInput(1,2) = 0.;
      input.gDInput(2,0) = 0.;
      input.gDInput(2,1) = 0.;
      input.gDInput(2,2) = 0.;
      
      input.KappaInput(0,0) = 2.0e-1;
      input.KappaInput(0,1) = 0.;
      input.KappaInput(0,2) = 0.;
      input.KappaInput(1,0) = 0.;
      input.KappaInput(1,1) = 2.0e-1;
      input.KappaInput(1,2) = 0.;
      input.KappaInput(2,0) = 0.;
      input.KappaInput(2,1) = 0.;
      input.KappaInput(2,2) = 2.0e-1;
      
      input.Lambda12Input(0,0) = 5.0e-1;
      input.Lambda12Input(0,1) = 0.;
      input.Lambda12Input(1,0) = 0.;
      input.Lambda12Input(1,1) = 5.0e-1;
      
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

   void set_minpar_values(CSE6SSM_scan_parameters params, const std::vector<std::size_t>& posn, CSE6SSM_input_parameters<Two_scale>& input)
   {
      if (params.get_is_grid_scan()) {

         input.m0 = params.get_m0_lower() + params.get_m0_incr() * posn.at(0);
         input.m12 = params.get_m12_lower() + params.get_m12_incr() * posn.at(1);
         input.TanBeta = params.get_TanBeta_lower() + params.get_TanBeta_incr() * posn.at(2);
         input.SignLambdax = params.get_SignLambdax_lower() + params.get_SignLambdax_incr() * posn.at(3);
         input.Azero = params.get_Azero_lower() + params.get_Azero_incr() * posn.at(4);
         input.KappaInput(0,0) = params.get_Kappa_lower() + params.get_Kappa_incr() * posn.at(5);
         input.KappaInput(1,1) = params.get_Kappa_lower() + params.get_Kappa_incr() * posn.at(5);
         input.KappaInput(2,2) = params.get_Kappa_lower() + params.get_Kappa_incr() * posn.at(5);
         input.Lambda12Input(0,0) = params.get_Lambda12_lower() + params.get_Lambda12_incr() * posn.at(6);
         input.Lambda12Input(1,1) = params.get_Lambda12_lower() + params.get_Lambda12_incr() * posn.at(6);
      } else {
         input.m0 = params.get_random_m0(generator);
         input.m12 = params.get_random_m12(generator);
         input.TanBeta = params.get_random_TanBeta(generator);
         input.SignLambdax = params.get_random_SignLambdax(generator);
         input.Azero = params.get_random_Azero(generator);
         const double random_Kappa = params.get_random_Kappa(generator);
         input.KappaInput(0,0) = random_Kappa;
         input.KappaInput(1,1) = random_Kappa;
         input.KappaInput(2,2) = random_Kappa;
         const double random_Lambda12 = params.get_random_Lambda12(generator);
         input.Lambda12Input(0,0) = random_Lambda12;
         input.Lambda12Input(1,1) = random_Lambda12;
      }
   }

   inline void trim(std::string& str)
   {
      std::size_t startpos = str.find_first_not_of(" \t\n\v\f\r");
      if (startpos != std::string::npos) str.erase(0, startpos);
      std::size_t endpos = str.find_last_not_of(" \t\n\v\f\r");
      if (endpos != std::string::npos) str.erase(endpos+1);
   }

   CSE6SSM_scan_parameters parse_scan_inputs_file(const std::string& scan_input_file)
   {
      std::ifstream ifs(scan_input_file, std::ifstream::in);
      if (ifs.fail()) {
         throw ReadError("unable to open file " + scan_input_file);
      }
      
      double m0_lower = 0.;
      double m0_upper = 0.;
      int m0_npts = 1;
      double m12_lower = 0.;
      double m12_upper = 0.;
      int m12_npts = 1;
      double TanBeta_lower = 1.;
      double TanBeta_upper = 1.;
      int TanBeta_npts = 1;
      int SignLambdax_lower = 1;
      int SignLambdax_upper = 1;
      int SignLambdax_npts = 1;
      double Azero_lower = 0.;
      double Azero_upper = 0.;
      int Azero_npts = 1;

      // extension to allow scanning over universal kappa
      // and lambda couplings (not very neat, but works for now)
      double Kappa_lower = 0.2;
      double Kappa_upper = 0.2;
      int Kappa_npts = 1;

      double Lambda12_lower = 0.5;
      double Lambda12_upper = 0.5;
      int Lambda12_npts = 1;

      int total_npts = 1;
      bool is_grid_scan = true;

      // read from file
      // # starts a (single line) comment
      // \n, ;, and , are delimiters
      std::string line;
      while (std::getline(ifs, line)) {
         if (line.find("#") != std::string::npos) {
            line = line.substr(0, line.find("#"));
         }
         if (line.empty())
            continue;

         // break up into individual fields
         std::vector<std::string> fields;
         boost::split(fields, line, boost::is_any_of(",;"));

         for (std::size_t i = 0; i < fields.size(); ++i) {
            // remove whitespace
            trim(fields[i]);

            // get field name and value
            if (!fields[i].empty()) {
               std::vector<std::string> field;
               boost::split(field, fields[i], boost::is_any_of("="));
               if (field.size() < 2) {
                  WARNING("Ignoring invalid input '" + fields[i] + "'");
               } else {
                  trim(field[0]);
                  trim(field[1]);

                  // compare against valid inputs
                  if (field[0] == "is_grid_scan") {
                     boost::to_lower(field[1]);
                     if (field[1] == "false" || field[1] == "f") {
                        is_grid_scan = false;
                     } else if (field[1] == "true" || field[1] == "t") {
                        is_grid_scan = true;
                     } else {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                     }
                  } else if (field[0] == "m0_lower") {
                     try {
                        m0_lower = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        m0_lower = 0.;
                     }
                  } else if (field[0] == "m0_upper") {
                     try {
                        m0_upper = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        m0_upper = 0.;
                     }
                  } else if (field[0] == "m0_npts") {
                     try {
                        m0_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (m0_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           m0_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        m0_npts = 1;
                     }
                  } else if (field[0] == "m12_lower") {
                     try {
                        m12_lower = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        m12_lower = 0.;
                     }
                  } else if (field[0] == "m12_upper") {
                     try {
                        m12_upper = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        m12_upper = 0.;
                     }
                  } else if (field[0] == "m12_npts") {
                     try {
                        m12_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (m12_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           m12_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        m12_npts = 1;
                     }
                  } else if (field[0] == "TanBeta_lower") {
                     try {
                        TanBeta_lower = boost::lexical_cast<double>(field[1]);
                        if (TanBeta_lower < 1. || TanBeta_lower > 1000.) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           TanBeta_lower = 1.;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        TanBeta_lower = 1.;
                     }
                  } else if (field[0] == "TanBeta_upper") {
                     try {
                        TanBeta_upper = boost::lexical_cast<double>(field[1]);
                        if (TanBeta_upper < 1. || TanBeta_upper > 1000.) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           TanBeta_upper = 1.;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        TanBeta_upper = 1.;
                     }
                  } else if (field[0] == "TanBeta_npts") {
                     try {
                        TanBeta_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (TanBeta_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           TanBeta_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        TanBeta_npts = 1;
                     }
                  } else if (field[0] == "SignLambdax_lower") {
                     try {
                        SignLambdax_lower = Sign(boost::lexical_cast<double>(field[1]));
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        SignLambdax_lower = -1;
                     }
                  } else if (field[0] == "SignLambdax_upper") {
                     try {
                        SignLambdax_upper = Sign(boost::lexical_cast<double>(field[1]));
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        SignLambdax_upper = 1;
                     }
                  } else if (field[0] == "SignLambdax_npts") {
                     try {
                        SignLambdax_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (SignLambdax_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           SignLambdax_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        SignLambdax_npts = 1;
                     }
                  } else if (field[0] == "Azero_lower") {
                     try {
                        Azero_lower = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Azero_lower = 0.;
                     }
                  } else if (field[0] == "Azero_upper") {
                     try {
                        Azero_upper = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Azero_upper = 0.;
                     }
                  } else if (field[0] == "Azero_npts") {
                     try {
                        Azero_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (Azero_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           Azero_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Azero_npts = 1;
                     }
                  } else if (field[0] == "Kappa_lower") {
                     try {
                        Kappa_lower = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Kappa_lower = 0.2;
                     }
                  } else if (field[0] == "Kappa_upper") {
                     try {
                        Kappa_upper = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Kappa_upper = 0.2;
                     }
                  } else if (field[0] == "Kappa_npts") {
                     try {
                        Kappa_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (Kappa_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           Kappa_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Kappa_npts = 1;
                     }
                  } else if (field[0] == "Lambda12_lower") {
                     try {
                        Lambda12_lower = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Lambda12_lower = 0.5;
                     }
                  } else if (field[0] == "Lambda12_upper") {
                     try {
                        Lambda12_upper = boost::lexical_cast<double>(field[1]);
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Lambda12_upper = 0.5;
                     }
                  } else if (field[0] == "Lambda12_npts") {
                     try {
                        Lambda12_npts = boost::lexical_cast<int>(field[1]);
                        if (!is_grid_scan) {
                           WARNING("Random scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (Lambda12_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           Lambda12_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        Lambda12_npts = 1;
                     }
                  } else if (field[0] == "total_npts") {
                     try {
                        total_npts = boost::lexical_cast<int>(field[1]);
                        if (is_grid_scan) {
                           WARNING("Grid scan requested, input '" + fields[i] + "' will be ignored");
                        } else if (total_npts <= 0) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           total_npts = 1;
                        }
                     } catch (const boost::bad_lexical_cast& error) {
                        WARNING("Ignoring invalid input '" + fields[i] + "'");
                        total_npts = 1;
                     }
                  } else {
                     WARNING("Ignoring invalid input '" + fields[i] + "'");
                  }
               } //< if (field.size() < 2)
            } //< if (!fields[i].empty())
         } //< for(std::size_t i = 0; i < fields.size(); ++i)
      } //< while(std::getline(ifs, line))

      // initialise scan parameters
      if (is_grid_scan) {
         return CSE6SSM_scan_parameters(m0_lower, m0_upper, m0_npts,
                                        m12_lower, m12_upper, m12_npts,
                                        TanBeta_lower, TanBeta_upper, TanBeta_npts,
                                        SignLambdax_lower, SignLambdax_upper, SignLambdax_npts,
                                        Azero_lower, Azero_upper, Azero_npts, Kappa_lower,
                                        Kappa_upper, Kappa_npts, Lambda12_lower, Lambda12_upper,
                                        Lambda12_npts);
      } else {
         return CSE6SSM_scan_parameters(m0_lower, m0_upper,
                                        m12_lower, m12_upper,
                                        TanBeta_lower, TanBeta_upper, 
                                        SignLambdax_lower, SignLambdax_upper,
                                        Azero_lower, Azero_upper, Kappa_lower, 
                                        Kappa_upper, Lambda12_lower, Lambda12_upper, total_npts);
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

//   std::size_t seed = start_point.time_since_epoch().count();

//   static std::minstd_rand generator(seed);
   std::uniform_real_distribution<double> distribution(0., 1.);

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

   CSE6SSM_scan_parameters parameters;

   if (scan_input_file.empty()) {
      WARNING("No scan input file given!\n"
              "   Default scan parameters will be used.\n"
              "   You can provide them via the option --scan-input-file=");
   } else {
      try {
         parameters = parse_scan_inputs_file(scan_input_file);
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

   std::ofstream drbar_mass_out_stream;
   std::ofstream drbar_susy_pars_out_stream;
   std::ofstream drbar_soft_pars_out_stream;
   std::ofstream drbar_mixings_out_stream;
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

   CSE6SSM_input_parameters<Two_scale> input;
   set_default_parameter_values(input);

   // attempt to read scan input file

   QedQcd oneset;
   oneset.toMz();

   CSE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-3);
   spectrum_generator.set_max_iterations(0);   // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0); // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale 

   std::vector<std::size_t> scan_dimensions = {parameters.get_m0_npts(), parameters.get_m12_npts(), 
                                               parameters.get_TanBeta_npts(), 
                                               parameters.get_SignLambdax_npts(), 
                                               parameters.get_Azero_npts(),
                                               parameters.get_Kappa_npts(),
                                               parameters.get_Lambda12_npts()};

   Grid_scanner scan(scan_dimensions);

   CSE6SSM_pole_mass_writer pole_mass_writer;
   CSE6SSM_drbar_values_writer drbar_values_writer;
   bool must_write_comment_line = true;
   while (!scan.has_finished()) {
      set_minpar_values(parameters, scan.get_position(), input);

      spectrum_generator.run(oneset, input);

      const CSE6SSM<algorithm_type>& model = spectrum_generator.get_model();

      pole_mass_writer.extract_pole_masses(model);

      if (must_write_drbar_masses)
         drbar_values_writer.extract_drbar_masses(model);
      if (must_write_drbar_susy_pars)
         drbar_values_writer.extract_drbar_susy_pars(model);
      if (must_write_drbar_soft_pars)
         drbar_values_writer.extract_drbar_soft_pars(model);
      if (must_write_drbar_mixings)
         drbar_values_writer.extract_drbar_mixings(model);

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

