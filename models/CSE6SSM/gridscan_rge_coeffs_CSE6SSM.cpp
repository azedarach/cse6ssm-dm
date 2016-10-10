// ====================================================================
// Does a grid scan of the CSE6SSM parameter space, printing
// out the coefficients in the expansions of the soft masses
// mHd2 and mHu2, and in the tree level expansion of Lambdax
// ====================================================================

#include "CSE6SSM_two_scale_input_parameters.hpp"
#include "CSE6SSM_scan_parameters.hpp"
#include "CSE6SSM_scan_utilities.hpp"
#include "CSE6SSM_spectrum_generator.hpp"
#include "CSE6SSM_slha_io.hpp"

#include "scan_command_line_options.hpp"
#include "error.hpp"
#include "grid_scanner.hpp"
#include "lowe.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <sys/time.h>
#include <map>
#include <vector>

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
      double output_scale = -1.;
      bool output_at_susy_scale = true;
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
                  } else if (field[0] == "output_scale") {
                     if (field[1] == "susy_scale") {
                        output_scale = -1.;
                        output_at_susy_scale = true;
                     } else {
                        try {
                           output_scale = boost::lexical_cast<double>(field[1]);
                           if (output_scale > 0.)
                              output_at_susy_scale = false;
                        } catch (const boost::bad_lexical_cast& error) {
                           WARNING("Ignoring invalid input '" + fields[i] + "'");
                           output_scale = -1.;
                           output_at_susy_scale = true;
                        }
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
         if (output_at_susy_scale) {
            return CSE6SSM_scan_parameters(m0_lower, m0_upper, m0_npts,
                                           m12_lower, m12_upper, m12_npts,
                                           TanBeta_lower, TanBeta_upper, TanBeta_npts,
                                           SignLambdax_lower, SignLambdax_upper, SignLambdax_npts,
                                           Azero_lower, Azero_upper, Azero_npts, Kappa_lower,
                                           Kappa_upper, Kappa_npts, Lambda12_lower, Lambda12_upper,
                                           Lambda12_npts);
         } else {
            return CSE6SSM_scan_parameters(m0_lower, m0_upper, m0_npts,
                                           m12_lower, m12_upper, m12_npts,
                                           TanBeta_lower, TanBeta_upper, TanBeta_npts,
                                           SignLambdax_lower, SignLambdax_upper, SignLambdax_npts,
                                           Azero_lower, Azero_upper, Azero_npts, Kappa_lower, 
                                           Kappa_upper, Kappa_npts, Lambda12_lower, Lambda12_upper, 
                                           Lambda12_npts, output_scale);
         }
      } else {
         if (output_at_susy_scale) {
            return CSE6SSM_scan_parameters(m0_lower, m0_upper,
                                           m12_lower, m12_upper,
                                           TanBeta_lower, TanBeta_upper, 
                                           SignLambdax_lower, SignLambdax_upper,
                                           Azero_lower, Azero_upper, Kappa_lower, 
                                           Kappa_upper, Lambda12_lower, Lambda12_upper,
                                           total_npts);
         } else {
            return CSE6SSM_scan_parameters(m0_lower, m0_upper,
                                           m12_lower, m12_upper,
                                           TanBeta_lower, TanBeta_upper, 
                                           SignLambdax_lower, SignLambdax_upper,
                                           Azero_lower, Azero_upper, Kappa_lower, 
                                           Kappa_upper, Lambda12_lower, Lambda12_upper,
                                           total_npts, output_scale);
         }
      }
   }

   double get_soft_parameter_from_inputs(CSE6SSM<Two_scale> model, CSE6SSM_info::Parameters p, 
                                         double m0, double m12, double Azero, double output_scale)
   {
      const auto Ye = model.get_Ye();
      const auto Yd = model.get_Yd();
      const auto Yu = model.get_Yu();
      const auto KappaPr = model.get_KappaPr();
      const auto Sigmax = model.get_Sigmax();
      const auto hE = model.get_hE();
      const auto SigmaL = model.get_SigmaL();
      const auto gD = model.get_gD();
      const auto fu = model.get_fu();
      const auto fd = model.get_fd();
      const auto Kappa = model.get_Kappa();
      const auto Lambda12 = model.get_Lambda12();
      const auto Lambdax = model.get_Lambdax();
      
      model.set_TYe(Azero * Ye);
      model.set_TYd(Azero * Yd);
      model.set_TYu(Azero * Yu);
      model.set_TKappaPr(Azero * KappaPr);
      model.set_TSigmax(Azero * Sigmax);
      model.set_ThE(Azero * hE);
      model.set_TSigmaL(Azero * SigmaL);
      model.set_TgD(Azero * gD);
      model.set_Tfu(Azero * fu);
      model.set_Tfd(Azero * fd);
      model.set_TKappa(Azero * Kappa);
      model.set_TLambda12(Azero * Lambda12);
      model.set_TLambdax(Azero * Lambdax);
      model.set_mHd2(Sqr(m0));
      model.set_mHu2(Sqr(m0));
      model.set_ms2(Sqr(m0));
      model.set_msbar2(Sqr(m0));
      model.set_mphi2(Sqr(m0));
      model.set_mHp2(Sqr(m0));
      model.set_mHpbar2(Sqr(m0));
      model.set_mH1I2(Sqr(m0)*UNITMATRIX(2));
      model.set_mH2I2(Sqr(m0)*UNITMATRIX(2));
      model.set_mSI2(Sqr(m0)*UNITMATRIX(3));
      model.set_mq2(Sqr(m0)*UNITMATRIX(3));
      model.set_ml2(Sqr(m0)*UNITMATRIX(3));
      model.set_md2(Sqr(m0)*UNITMATRIX(3));
      model.set_mu2(Sqr(m0)*UNITMATRIX(3));
      model.set_me2(Sqr(m0)*UNITMATRIX(3));
      model.set_mDx2(Sqr(m0)*UNITMATRIX(3));
      model.set_mDxbar2(Sqr(m0)*UNITMATRIX(3));
      model.set_MassB(m12);
      model.set_MassWB(m12);
      model.set_MassG(m12);
      model.set_MassBp(m12);
      
      model.run_to(output_scale);
      
      return model.get_parameter(p);
   }

   double get_tree_level_Lambdax_soft_term(CSE6SSM<Two_scale> model, double mHd2_coeff, double mHu2_coeff)
   {
      const auto vd = model.get_vd();
      const auto vu = model.get_vu();
      const auto vs = model.get_vs();
      
      double coeff = 2.0 * (Sqr(vd) * mHd2_coeff - Sqr(vu) * mHu2_coeff) / (Sqr(vs) * (Sqr(vu) - Sqr(vd)));
      
      return coeff;
   }
   
   double get_tree_level_Lambdax_constant_term(CSE6SSM<Two_scale> model)
   {
      const auto g1 = model.get_g1();
      const auto g2 = model.get_g2();
      const auto g1p = model.get_g1p();
      const auto vd = model.get_vd();
      const auto vu = model.get_vu();
      const auto vs = model.get_vs();
      const auto vsb = model.get_vsb();
      const auto QS = model.get_input().QSInput;
      
      double coeff = - 0.25 * (Sqr(g2) + 0.6 * Sqr(g1)) * (Sqr(vd) + Sqr(vu)) / Sqr(vs)
         + 0.025 * Sqr(g1p) * (2.0 * Sqr(vu) - 3.0 * Sqr(vd)) * 
         (-3.0 * Sqr(vd) - 2.0 * Sqr(vu) + QS * (Sqr(vs) - Sqr(vsb))) / (Sqr(vs) * (Sqr(vu) - Sqr(vd)));
      
      return coeff;
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
   bool must_write_pole_masses = false;

   std::ofstream pole_mass_out_stream;
   if (!pole_mass_output_file.empty()) {
      must_write_pole_masses = true;
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

   // soft scalar masses to calculate coefficients for
   std::vector<CSE6SSM_info::Parameters> soft_scalar_masses 
      = {CSE6SSM_info::mHd2, CSE6SSM_info::mHu2, CSE6SSM_info::ms2,
         CSE6SSM_info::msbar2, CSE6SSM_info::mphi2, CSE6SSM_info::mq200,
         CSE6SSM_info::mu200, CSE6SSM_info::md200};

   bool is_calculating_mHd2_coeffs = false;
   bool is_calculating_mHu2_coeffs = false;
   bool can_calculate_Lambdax_coeffs = false;
   for (std::size_t i = 0; i < soft_scalar_masses.size(); ++i) {
      if (soft_scalar_masses[i] == CSE6SSM_info::mHd2) {
         is_calculating_mHd2_coeffs = true;
      } else if (soft_scalar_masses[i] == CSE6SSM_info::mHu2) {
         is_calculating_mHu2_coeffs = true;
      } 

      if (is_calculating_mHd2_coeffs && is_calculating_mHu2_coeffs) {
         can_calculate_Lambdax_coeffs = true;
         break;
      }
   }

   // soft gaugino masses to calculate coefficients for
   std::vector<CSE6SSM_info::Parameters> soft_gaugino_masses
      = {CSE6SSM_info::MassB, CSE6SSM_info::MassWB, CSE6SSM_info::MassG, 
         CSE6SSM_info::MassBp};
   
   // soft trilinears to calculate coefficients for
   std::vector<CSE6SSM_info::Parameters> soft_trilinears
      = {CSE6SSM_info::TYu22, CSE6SSM_info::TYu00, CSE6SSM_info::TYu11, 
         CSE6SSM_info::TYd22, CSE6SSM_info::TYd11, CSE6SSM_info::TYd00, 
         CSE6SSM_info::TSigmax, CSE6SSM_info::TLambdax, CSE6SSM_info::TKappa22,
         CSE6SSM_info::TLambda1211};

   // print comment line to standard output for coefficients
   std::size_t width = 18;
   std::cout << "# "
             << std::left << std::setw(width) << "m0/GeV" << ' '
             << std::left << std::setw(width) << "m12/GeV" << ' '
             << std::left << std::setw(width) << "TanBeta" << ' '
             << std::left << std::setw(width) << "SignLambdax" << ' '
             << std::left << std::setw(width) << "Azero/GeV" << ' '
             << std::left << std::setw(width) << "sInput/GeV" << ' '
             << std::left << std::setw(width) << "QSInput" << ' '
             << std::left << std::setw(width) << "hEInput(0,0)" << ' '
             << std::left << std::setw(width) << "hEInput(1,0)" << ' '
             << std::left << std::setw(width) << "hEInput(2,0)" << ' '
             << std::left << std::setw(width) << "hEInput(0,1)" << ' '
             << std::left << std::setw(width) << "hEInput(1,1)" << ' '
             << std::left << std::setw(width) << "hEInput(2,1)" << ' '
             << std::left << std::setw(width) << "SigmaLInput" << ' '
             << std::left << std::setw(width) << "KappaPrInput" << ' '
             << std::left << std::setw(width) << "SigmaxInput" << ' '
             << std::left << std::setw(width) << "gDInput(0,0)" << ' '
             << std::left << std::setw(width) << "gDInput(1,0)" << ' '
             << std::left << std::setw(width) << "gDInput(2,0)" << ' '
             << std::left << std::setw(width) << "gDInput(0,1)" << ' '
             << std::left << std::setw(width) << "gDInput(1,1)" << ' '
             << std::left << std::setw(width) << "gDInput(2,1)" << ' '
             << std::left << std::setw(width) << "gDInput(0,2)" << ' '
             << std::left << std::setw(width) << "gDInput(1,2)" << ' '
             << std::left << std::setw(width) << "gDInput(2,2)" << ' '
             << std::left << std::setw(width) << "KappaInput(0,0)" << ' '
             << std::left << std::setw(width) << "KappaInput(1,0)" << ' '
             << std::left << std::setw(width) << "KappaInput(2,0)" << ' '
             << std::left << std::setw(width) << "KappaInput(0,1)" << ' '
             << std::left << std::setw(width) << "KappaInput(1,1)" << ' '
             << std::left << std::setw(width) << "KappaInput(2,1)" << ' '
             << std::left << std::setw(width) << "KappaInput(0,2)" << ' '
             << std::left << std::setw(width) << "KappaInput(1,2)" << ' '
             << std::left << std::setw(width) << "KappaInput(2,2)" << ' '
             << std::left << std::setw(width) << "Lambda12Input(0,0)" << ' '
             << std::left << std::setw(width) << "Lambda12Input(1,0)" << ' '
             << std::left << std::setw(width) << "Lambda12Input(0,1)" << ' '
             << std::left << std::setw(width) << "Lambda12Input(1,1)" << ' '
             << std::left << std::setw(width) << "fuInput(0,0)" << ' '
             << std::left << std::setw(width) << "fuInput(1,0)" << ' '
             << std::left << std::setw(width) << "fuInput(2,0)" << ' '
             << std::left << std::setw(width) << "fuInput(0,1)" << ' '
             << std::left << std::setw(width) << "fuInput(1,1)" << ' '
             << std::left << std::setw(width) << "fuInput(2,1)" << ' '
             << std::left << std::setw(width) << "fdInput(0,0)" << ' '
             << std::left << std::setw(width) << "fdInput(1,0)" << ' '
             << std::left << std::setw(width) << "fdInput(2,0)" << ' '
             << std::left << std::setw(width) << "fdInput(0,1)" << ' '
             << std::left << std::setw(width) << "fdInput(1,1)" << ' '
             << std::left << std::setw(width) << "fdInput(2,1)" << ' '
             << std::left << std::setw(width) << "MuPrInput/GeV" << ' '
             << std::left << std::setw(width) << "MuPhiInput/GeV" << ' '
             << std::left << std::setw(width) << "BMuPrInput/GeV^2" << ' '
             << std::left << std::setw(width) << "BMuPhiInput/GeV^2" << ' '
             << std::left << std::setw(width) << "Q/GeV" << ' '
             << std::left << std::setw(width) << "g1(Q)" << ' '
             << std::left << std::setw(width) << "g2(Q)" << ' '
             << std::left << std::setw(width) << "g3(Q)" << ' '
             << std::left << std::setw(width) << "g1p(Q)" << ' '
             << std::left << std::setw(width) << "Lambdax(Q)" << ' '
             << std::left << std::setw(width) << "vd(Q)/GeV" << ' '
             << std::left << std::setw(width) << "vu(Q)/GeV" << ' '
             << std::left << std::setw(width) << "vs(Q)/GeV" << ' '
             << std::left << std::setw(width) << "vsb(Q)/GeV" << ' '
             << std::left << std::setw(width) << "vphi(Q)/GeV" << ' ';
   for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_scalar_masses.begin(),
           end = soft_scalar_masses.end(); it != end; ++it) {
      const std::string scalar_mass_name(CSE6SSM_info::parameter_names[*it]);
      std::cout << std::left << std::setw(width) << "coeff:a(" + scalar_mass_name + ",Q)" << ' '
                << std::left << std::setw(width) << "coeff:b(" + scalar_mass_name + ",Q)" << ' '
                << std::left << std::setw(width) << "coeff:c(" + scalar_mass_name + ",Q)" << ' '
                << std::left << std::setw(width) << "coeff:d(" + scalar_mass_name + ",Q)" << ' '
                << std::left << std::setw(width) << scalar_mass_name + "(Q)" << ' '
                << std::left << std::setw(width) << scalar_mass_name + "Approx(Q)" << ' ';
   }
   for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_gaugino_masses.begin(),
           end = soft_gaugino_masses.end(); it != end; ++it) {
      const std::string gaugino_mass_name(CSE6SSM_info::parameter_names[*it]);
      std::cout << std::left << std::setw(width) << "coeff:p(" + gaugino_mass_name + ",Q)" << ' '
                << std::left << std::setw(width) << "coeff:q(" + gaugino_mass_name + ",Q)" << ' '
                << std::left << std::setw(width) << gaugino_mass_name + "(Q)" << ' '
                << std::left << std::setw(width) << gaugino_mass_name + "Approx(Q)" << ' ';
   }
   for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_trilinears.begin(),
           end = soft_trilinears.end(); it != end; ++it) {
      const std::string trilinear_name(CSE6SSM_info::parameter_names[*it]);
      std::cout << std::left << std::setw(width) << "coeff:e(" + trilinear_name + ",Q)" << ' '
                << std::left << std::setw(width) << "coeff:f(" + trilinear_name + ",Q)" << ' '
                << std::left << std::setw(width) << trilinear_name + "(Q)" << ' '
                << std::left << std::setw(width) << trilinear_name + "Approx(Q)" << ' ';
   }
   if (can_calculate_Lambdax_coeffs) {
      std::cout << std::left << std::setw(width) << "aLambdax(Q)" << ' '
                << std::left << std::setw(width) << "bLambdax(Q)" << ' '
                << std::left << std::setw(width) << "cLambdax(Q)" << ' '
                << std::left << std::setw(width) << "dLambdax(Q)" << ' '
                << std::left << std::setw(width) << "lLambdax(Q)" << ' '
                << std::left << std::setw(width) << "Lambdax0lp(Q)" << ' '
                << std::left << std::setw(width) << "Lambdax0lpApprox(Q)" << ' ';
   }
   std::cout << std::left << std::setw(width) << "error" << '\n';
   
   // note
   std::size_t count = 0;
   while (!scan.has_finished()) {
      ++count;
      set_minpar_values(parameters, scan.get_position(), input);

      QedQcd oneset;
      oneset.toMz();

      CSE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
      spectrum_generator.set_precision_goal(1.0e-3);
      spectrum_generator.set_max_iterations(0);   // 0 == automatic
      spectrum_generator.set_calculate_sm_masses(0); // 0 == no
      // note
      spectrum_generator.set_threshold_corrections_loop_order(1);
      if (parameters.get_output_scale() <= 0.) {
         spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale 
      } else {
         spectrum_generator.set_parameter_output_scale(parameters.get_output_scale());
      }
      
      spectrum_generator.run(oneset, input);

      const CSE6SSM<algorithm_type>& model = spectrum_generator.get_model();

      // calculate coefficients
      const Problems<CSE6SSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();

      const bool error = problems.have_problem();

      // note
      if (!error) {
         std::string slha_file = "/home/dylan/Documents/Postgraduate/CSE6SSM-Spectrum/scans/cne6ssm_small_lambdax_searches/slha_files/cne6ssm_small_lambdax_m0_scan_7_slha_" + boost::lexical_cast<std::string>(count) + ".txt";
         CSE6SSM_slha_io slha_io;

         slha_io.set_spinfo(problems);
         slha_io.set_sminputs(oneset);
         slha_io.set_minpar(input);
         slha_io.set_extpar(input);
         slha_io.set_spectrum(model);
         slha_io.write_to_file(slha_file);
      }

      double high_scale = spectrum_generator.get_high_scale();
      double susy_scale = spectrum_generator.get_susy_scale();

      double output_scale = susy_scale;
      if (parameters.get_output_scale() > 0.) {
         output_scale = parameters.get_output_scale();
      }

      double aLambdax = 0.;
      double bLambdax = 0.;
      double cLambdax = 0.;
      double dLambdax = 0.;
      double lLambdax = 0.;

      const double coeffs_precision = 0.5;
      bool coeffs_error = false;

      double tree_level_Lambdax = 0.;

      std::map<CSE6SSM_info::Parameters, std::vector<double> > soft_scalar_mass_coeffs;
      std::map<CSE6SSM_info::Parameters, double> soft_scalar_mass_values;
      std::map<CSE6SSM_info::Parameters, double> soft_scalar_mass_errors;
      for (std::size_t i = 0; i < soft_scalar_masses.size(); ++i) {
         soft_scalar_mass_coeffs[soft_scalar_masses[i]] = {0., 0., 0., 0.};
         soft_scalar_mass_values[soft_scalar_masses[i]] = 0.;
         soft_scalar_mass_errors[soft_scalar_masses[i]] = 0.;
      }

      std::map<CSE6SSM_info::Parameters, std::vector<double> > soft_gaugino_mass_coeffs;
      std::map<CSE6SSM_info::Parameters, double> soft_gaugino_mass_values;
      std::map<CSE6SSM_info::Parameters, double> soft_gaugino_mass_errors;
      for (std::size_t i = 0; i < soft_gaugino_masses.size(); ++i) {
         soft_gaugino_mass_coeffs[soft_gaugino_masses[i]] = {0., 0.};
         soft_gaugino_mass_values[soft_gaugino_masses[i]] = 0.;
         soft_gaugino_mass_errors[soft_gaugino_masses[i]] = 0.;
      }

      std::map<CSE6SSM_info::Parameters, std::vector<double> > soft_trilinear_coeffs;
      std::map<CSE6SSM_info::Parameters, double> soft_trilinear_values;
      std::map<CSE6SSM_info::Parameters, double> soft_trilinear_errors;
      for (std::size_t i = 0; i < soft_trilinears.size(); ++i) {
         soft_trilinear_coeffs[soft_trilinears[i]] = {0., 0.};
         soft_trilinear_values[soft_trilinears[i]] = 0.;
         soft_trilinear_errors[soft_trilinears[i]] = 0.;
      }

      if (!error) {
         // calculate coefficients and percentage errors the original way
         for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_scalar_masses.begin(),
                 end = soft_scalar_masses.end(); it != end; ++it) {
            const Eigen::Array<double,4,1> coeffs = model.get_soft_scalar_mass_coeffs(*it, susy_scale, high_scale);
            soft_scalar_mass_coeffs[*it] = {coeffs(0), coeffs(1), coeffs(2), coeffs(3)};
            soft_scalar_mass_values[*it] = coeffs(0) * Sqr(model.get_input().m0) + coeffs(1) * Sqr(model.get_input().m12)
               + coeffs(2) * model.get_input().m12 * model.get_input().Azero + coeffs(3) * Sqr(model.get_input().Azero);
            soft_scalar_mass_errors[*it] = 100.0 * Abs((model.get_parameter(*it) - soft_scalar_mass_values[*it]) / 
                                                       (0.5 * (model.get_parameter(*it) + soft_scalar_mass_values[*it])));
            if (soft_scalar_mass_errors[*it] > coeffs_precision) coeffs_error = true; 
         }
         
         for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_gaugino_masses.begin(),
                 end = soft_gaugino_masses.end(); it != end; ++it) {
            const Eigen::Array<double,2,1> coeffs = model.get_soft_gaugino_mass_coeffs(*it, susy_scale, high_scale);
            soft_gaugino_mass_coeffs[*it] = {coeffs(0), coeffs(1)};
            soft_gaugino_mass_values[*it] = coeffs(0) * model.get_input().Azero + coeffs(1) * model.get_input().m12;
            soft_gaugino_mass_errors[*it] = 100.0 * Abs((model.get_parameter(*it) - soft_gaugino_mass_values[*it]) / 
                                                        (0.5 * (model.get_parameter(*it) + soft_gaugino_mass_values[*it])));
            if (soft_gaugino_mass_errors[*it] > coeffs_precision) coeffs_error = true; 
         }
         
         for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_trilinears.begin(),
                 end = soft_trilinears.end(); it != end; ++it) {
            const Eigen::Array<double,2,1> coeffs = model.get_soft_trilinear_coeffs(*it, susy_scale, high_scale);
            soft_trilinear_coeffs[*it] = {coeffs(0), coeffs(1)};
            soft_trilinear_values[*it] = coeffs(0) * model.get_input().Azero + coeffs(1) * model.get_input().m12;
            soft_trilinear_errors[*it] = 100.0 * Abs((model.get_parameter(*it) - soft_trilinear_values[*it]) /
                                                     (0.5 * (model.get_parameter(*it) + soft_trilinear_values[*it])));
            if (soft_trilinear_errors[*it] > coeffs_precision) coeffs_error = true; 
         }

         CSE6SSM<algorithm_type> running_model(model);

         running_model.run_to(high_scale);

         if (can_calculate_Lambdax_coeffs) {
            aLambdax = get_tree_level_Lambdax_soft_term(model, soft_scalar_mass_coeffs[CSE6SSM_info::mHd2][0], 
                                                        soft_scalar_mass_coeffs[CSE6SSM_info::mHu2][0]);
            bLambdax = get_tree_level_Lambdax_soft_term(model, soft_scalar_mass_coeffs[CSE6SSM_info::mHd2][1],
                                                        soft_scalar_mass_coeffs[CSE6SSM_info::mHu2][1]);
            cLambdax = get_tree_level_Lambdax_soft_term(model, soft_scalar_mass_coeffs[CSE6SSM_info::mHd2][2],
                                                        soft_scalar_mass_coeffs[CSE6SSM_info::mHu2][2]);
            dLambdax = get_tree_level_Lambdax_soft_term(model, soft_scalar_mass_coeffs[CSE6SSM_info::mHd2][3],
                                                        soft_scalar_mass_coeffs[CSE6SSM_info::mHu2][3]);
            lLambdax = get_tree_level_Lambdax_constant_term(model);

            running_model.run_to(susy_scale);
            running_model.solve_ewsb_tree_level();

            if (parameters.get_output_scale() > 0.) {
               running_model.run_to(output_scale);
            }
            
            tree_level_Lambdax = running_model.get_Lambdax();
         }

      }

      if (must_write_pole_masses)
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
         if (must_write_pole_masses)
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
      if (must_write_pole_masses)
         pole_mass_writer.write_pole_masses_line(pole_mass_out);

      if (must_write_drbar_masses)
         drbar_values_writer.write_drbar_masses_line(drbar_mass_out_stream);
      if (must_write_drbar_susy_pars)
         drbar_values_writer.write_drbar_susy_pars_line(drbar_susy_pars_out_stream);
      if (must_write_drbar_soft_pars)
         drbar_values_writer.write_drbar_soft_pars_line(drbar_soft_pars_out_stream);
      if (must_write_drbar_mixings)
         drbar_values_writer.write_drbar_mixings_line(drbar_mixings_out_stream);

      // write output
      std::cout << "  "
                << std::left << std::setw(width) << input.m0 << ' '
                << std::left << std::setw(width) << input.m12 << ' '
                << std::left << std::setw(width) << input.TanBeta << ' '
                << std::left << std::setw(width) << input.SignLambdax << ' '
                << std::left << std::setw(width) << input.Azero << ' '
                << std::left << std::setw(width) << input.sInput << ' '
                << std::left << std::setw(width) << input.QSInput << ' '
                << std::left << std::setw(width) << input.hEInput(0,0) << ' '
                << std::left << std::setw(width) << input.hEInput(1,0) << ' '
                << std::left << std::setw(width) << input.hEInput(2,0) << ' '
                << std::left << std::setw(width) << input.hEInput(0,1) << ' '
                << std::left << std::setw(width) << input.hEInput(1,1) << ' '
                << std::left << std::setw(width) << input.hEInput(2,1) << ' '
                << std::left << std::setw(width) << input.SigmaLInput << ' '
                << std::left << std::setw(width) << input.KappaPrInput << ' '
                << std::left << std::setw(width) << input.SigmaxInput << ' '
                << std::left << std::setw(width) << input.gDInput(0,0) << ' '
                << std::left << std::setw(width) << input.gDInput(1,0) << ' '
                << std::left << std::setw(width) << input.gDInput(2,0) << ' '
                << std::left << std::setw(width) << input.gDInput(0,1) << ' '
                << std::left << std::setw(width) << input.gDInput(1,1) << ' '
                << std::left << std::setw(width) << input.gDInput(2,1) << ' '
                << std::left << std::setw(width) << input.gDInput(0,2) << ' '
                << std::left << std::setw(width) << input.gDInput(1,2) << ' '
                << std::left << std::setw(width) << input.gDInput(2,2) << ' '
                << std::left << std::setw(width) << input.KappaInput(0,0) << ' '
                << std::left << std::setw(width) << input.KappaInput(1,0) << ' '
                << std::left << std::setw(width) << input.KappaInput(2,0) << ' '
                << std::left << std::setw(width) << input.KappaInput(0,1) << ' '
                << std::left << std::setw(width) << input.KappaInput(1,1) << ' '
                << std::left << std::setw(width) << input.KappaInput(2,1) << ' '
                << std::left << std::setw(width) << input.KappaInput(0,2) << ' '
                << std::left << std::setw(width) << input.KappaInput(1,2) << ' '
                << std::left << std::setw(width) << input.KappaInput(2,2) << ' '
                << std::left << std::setw(width) << input.Lambda12Input(0,0) << ' '
                << std::left << std::setw(width) << input.Lambda12Input(1,0) << ' '
                << std::left << std::setw(width) << input.Lambda12Input(0,1) << ' '
                << std::left << std::setw(width) << input.Lambda12Input(1,1) << ' '
                << std::left << std::setw(width) << input.fuInput(0,0) << ' '
                << std::left << std::setw(width) << input.fuInput(1,0) << ' '
                << std::left << std::setw(width) << input.fuInput(2,0) << ' '
                << std::left << std::setw(width) << input.fuInput(0,1) << ' '
                << std::left << std::setw(width) << input.fuInput(1,1) << ' '
                << std::left << std::setw(width) << input.fuInput(2,1) << ' '
                << std::left << std::setw(width) << input.fdInput(0,0) << ' '
                << std::left << std::setw(width) << input.fdInput(1,0) << ' '
                << std::left << std::setw(width) << input.fdInput(2,0) << ' '
                << std::left << std::setw(width) << input.fdInput(0,1) << ' '
                << std::left << std::setw(width) << input.fdInput(1,1) << ' '
                << std::left << std::setw(width) << input.fdInput(2,1) << ' '
                << std::left << std::setw(width) << input.MuPrInput << ' '
                << std::left << std::setw(width) << input.MuPhiInput << ' '
                << std::left << std::setw(width) << input.BMuPrInput << ' '
                << std::left << std::setw(width) << input.BMuPhiInput << ' '
                << std::left << std::setw(width) << output_scale << ' '
                << std::left << std::setw(width) << model.get_g1() << ' '
                << std::left << std::setw(width) << model.get_g2() << ' '
                << std::left << std::setw(width) << model.get_g3() << ' '
                << std::left << std::setw(width) << model.get_g1p() << ' '
                << std::left << std::setw(width) << model.get_Lambdax() << ' '
                << std::left << std::setw(width) << model.get_vd() << ' '
                << std::left << std::setw(width) << model.get_vu() << ' '
                << std::left << std::setw(width) << model.get_vs() << ' '
                << std::left << std::setw(width) << model.get_vsb() << ' '
                << std::left << std::setw(width) << model.get_vphi() << ' ';
      for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_scalar_masses.begin(),
              end = soft_scalar_masses.end(); it != end; ++it) {
         std::cout << std::left << std::setw(width) << soft_scalar_mass_coeffs[*it][0]<< ' '
                   << std::left << std::setw(width) << soft_scalar_mass_coeffs[*it][1]<< ' '
                   << std::left << std::setw(width) << soft_scalar_mass_coeffs[*it][2]<< ' '
                   << std::left << std::setw(width) << soft_scalar_mass_coeffs[*it][3]<< ' '
                   << std::left << std::setw(width) << model.get_parameter(*it) << ' '
                   << std::left << std::setw(width) << soft_scalar_mass_values[*it] << ' ';
      }
      for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_gaugino_masses.begin(),
              end = soft_gaugino_masses.end(); it != end; ++it) {
         std::cout << std::left << std::setw(width) << soft_gaugino_mass_coeffs[*it][0] << ' '
                   << std::left << std::setw(width) << soft_gaugino_mass_coeffs[*it][1] << ' '
                   << std::left << std::setw(width) << model.get_parameter(*it) << ' '
                   << std::left << std::setw(width) << soft_gaugino_mass_values[*it] << ' ';
      }
      for (std::vector<CSE6SSM_info::Parameters>::const_iterator it = soft_trilinears.begin(),
              end = soft_trilinears.end(); it != end; ++it) {
         std::cout << std::left << std::setw(width) << soft_trilinear_coeffs[*it][0] << ' '
                   << std::left << std::setw(width) << soft_trilinear_coeffs[*it][1] << ' '
                   << std::left << std::setw(width) << model.get_parameter(*it) << ' '
                   << std::left << std::setw(width) << soft_trilinear_values[*it] << ' ';
      }
      if (can_calculate_Lambdax_coeffs) {
         std::cout << std::left << std::setw(width) << aLambdax << ' '
                   << std::left << std::setw(width) << bLambdax << ' '
                   << std::left << std::setw(width) << cLambdax << ' '
                   << std::left << std::setw(width) << dLambdax << ' '
                   << std::left << std::setw(width) << lLambdax << ' '
                   << std::left << std::setw(width) << tree_level_Lambdax << ' '
                   << std::left << std::setw(width) << input.SignLambdax * 
            Sqrt(aLambdax * Sqr(model.get_input().m0) + bLambdax * Sqr(model.get_input().m12) 
                 + cLambdax * model.get_input().m12 * model.get_input().Azero + dLambdax * Sqr(model.get_input().Azero) 
                 + lLambdax) << ' ';
      }
      std::cout << std::left << std::setw(width) << (error || coeffs_error);
      
      if (error || coeffs_error) {
         if (error && !coeffs_error) {
            std::cout << "\t# " << problems;
         } else if (!error && coeffs_error) {
            std::cout << "\t# coefficients inaccurate";
         } else {
            std::cout << "\t# coefficients inaccurate, " << problems;
         }
      }

      std::cout << '\n';

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
