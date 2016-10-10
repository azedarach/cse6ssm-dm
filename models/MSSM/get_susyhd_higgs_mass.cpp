
#include "CMSSM_semi_two_scale_input_parameters.hpp"
#include "MSSM_mass_eigenstates.hpp"

#include "command_line_options.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "susyhd_call.hpp"
#include "wrappers.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: get_susyhd_higgs_mass.x [options]\n"
      "Options:\n"
      "  --susy-pars-file=<name>\n"
      "  --soft-pars-file=<name>\n"
      "  --help,-h                         print this help message"
             << std::endl;
}

void trim(std::string& str)
{
   std::size_t startpos = str.find_first_not_of(" \t\n\v\f\r");
   if (startpos != std::string::npos) str.erase(0, startpos);
   std::size_t endpos = str.find_last_not_of(" \t\n\v\f\r");
   if (endpos != std::string::npos) str.erase(endpos+1);
}

void strip_units(std::string& str)
{
   std::size_t unitspos = str.find_first_of("/");
   if (unitspos != std::string::npos) str.erase(unitspos);
}

// removes leading comment character
void strip_comment_character(std::string& str)
{
   std::size_t comment_start = str.find_first_of("#");
   if (comment_start != std::string::npos) str.erase(0, comment_start + 1);
}

void read_column_headers(const std::string& susy_pars_header, const std::string& soft_pars_header, std::map<std::string,std::size_t>& columns)
{
   std::vector<std::string> susy_columns;
   boost::split(susy_columns, susy_pars_header, boost::is_any_of(" \t\n"),
                boost::token_compress_on);
   for (std::size_t i = 0; i < susy_columns.size(); ++i) {
      trim(susy_columns[i]);
      strip_units(susy_columns[i]);
      columns[susy_columns[i]] = i;
   }
   std::vector<std::string> soft_columns;
   boost::split(soft_columns, soft_pars_header, boost::is_any_of(" \t\n"),
                boost::token_compress_on);
   for (std::size_t i = 0; i < soft_columns.size(); ++i) {
      trim(soft_columns[i]);
      strip_units(soft_columns[i]);
      columns[soft_columns[i]] = i;
   }
}

CMSSM_semianalytic_input_parameters<Two_scale> read_inputs_from_datapoint(const std::string& susy_pars_line,
                                                    const std::map<std::string,std::size_t>& columns)
{
   CMSSM_semianalytic_input_parameters<Two_scale> inputs;

   // read and set input parameters
   std::vector<std::string> input_parameters;
   boost::split(input_parameters, susy_pars_line, boost::is_any_of(" \t\n"),
                boost::token_compress_on);

   inputs.m12 = boost::lexical_cast<double>(input_parameters[columns.at("m12")]);;
   inputs.Azero = boost::lexical_cast<double>(input_parameters[columns.at("Azero")]);
   inputs.TanBeta = boost::lexical_cast<double>(input_parameters[columns.at("TanBeta")]);
   inputs.MuInput = boost::lexical_cast<double>(input_parameters[columns.at("MuInput")]);

   return inputs;
}

MSSM_mass_eigenstates initialize_model_from_datapoint(
   const std::string& susy_pars_line, const std::string& soft_pars_line,
   const std::map<std::string,std::size_t>& columns) {

   MSSM_mass_eigenstates model;

   // read and set SUSY parameters
   std::vector<std::string> susy_parameters;
   boost::split(susy_parameters, susy_pars_line, boost::is_any_of(" \t\n"),
                boost::token_compress_on);

   model.set_Yd(0, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(1,1)")]));
   model.set_Yd(0, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(1,2)")]));
   model.set_Yd(0, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(1,3)")]));
   model.set_Yd(1, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(2,1)")]));
   model.set_Yd(1, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(2,2)")]));
   model.set_Yd(1, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(2,3)")]));
   model.set_Yd(2, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(3,1)")]));
   model.set_Yd(2, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(3,2)")]));
   model.set_Yd(2, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Yd(3,3)")]));

   model.set_Ye(0, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(1,1)")]));
   model.set_Ye(0, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(1,2)")]));
   model.set_Ye(0, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(1,3)")]));
   model.set_Ye(1, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(2,1)")]));
   model.set_Ye(1, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(2,2)")]));
   model.set_Ye(1, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(2,3)")]));
   model.set_Ye(2, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(3,1)")]));
   model.set_Ye(2, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(3,2)")]));
   model.set_Ye(2, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Ye(3,3)")]));

   model.set_Yu(0, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(1,1)")]));
   model.set_Yu(0, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(1,2)")]));
   model.set_Yu(0, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(1,3)")]));
   model.set_Yu(1, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(2,1)")]));
   model.set_Yu(1, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(2,2)")]));
   model.set_Yu(1, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(2,3)")]));
   model.set_Yu(2, 0, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(3,1)")]));
   model.set_Yu(2, 1, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(3,2)")]));
   model.set_Yu(2, 2, boost::lexical_cast<double>(susy_parameters[columns.at("Yu(3,3)")]));

   model.set_Mu(boost::lexical_cast<double>(susy_parameters[columns.at("Mu")]));

   if (columns.count("g1") != 0) {
      model.set_g1(boost::lexical_cast<double>(susy_parameters[columns.at("g1")]));
   } else if (columns.count("gY") != 0) {
      model.set_g1(Sqrt(5.0 / 3.0) * boost::lexical_cast<double>(susy_parameters[columns.at("gY")]));
   }
   model.set_g2(boost::lexical_cast<double>(susy_parameters[columns.at("g2")]));
   model.set_g3(boost::lexical_cast<double>(susy_parameters[columns.at("g3")]));

   model.set_vd(boost::lexical_cast<double>(susy_parameters[columns.at("vd")]));
   model.set_vu(boost::lexical_cast<double>(susy_parameters[columns.at("vu")]));

   // read and set soft parameters
   std::vector<std::string> soft_parameters;
   boost::split(soft_parameters, soft_pars_line, boost::is_any_of(" \t\n"),
                boost::token_compress_on);

   model.set_TYd(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(1,1)")]));
   model.set_TYd(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(1,2)")]));
   model.set_TYd(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(1,3)")]));
   model.set_TYd(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(2,1)")]));
   model.set_TYd(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(2,2)")]));
   model.set_TYd(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(2,3)")]));
   model.set_TYd(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(3,1)")]));
   model.set_TYd(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(3,2)")]));
   model.set_TYd(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYd(3,3)")]));

   model.set_TYe(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(1,1)")]));
   model.set_TYe(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(1,2)")]));
   model.set_TYe(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(1,3)")]));
   model.set_TYe(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(2,1)")]));
   model.set_TYe(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(2,2)")]));
   model.set_TYe(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(2,3)")]));
   model.set_TYe(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(3,1)")]));
   model.set_TYe(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(3,2)")]));
   model.set_TYe(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYe(3,3)")]));

   model.set_TYu(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(1,1)")]));
   model.set_TYu(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(1,2)")]));
   model.set_TYu(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(1,3)")]));
   model.set_TYu(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(2,1)")]));
   model.set_TYu(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(2,2)")]));
   model.set_TYu(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(2,3)")]));
   model.set_TYu(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(3,1)")]));
   model.set_TYu(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(3,2)")]));
   model.set_TYu(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("TYu(3,3)")]));

   model.set_BMu(boost::lexical_cast<double>(soft_parameters[columns.at("BMu")]));

   model.set_mq2(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(1,1)")]));
   model.set_mq2(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(1,2)")]));
   model.set_mq2(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(1,3)")]));
   model.set_mq2(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(2,1)")]));
   model.set_mq2(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(2,2)")]));
   model.set_mq2(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(2,3)")]));
   model.set_mq2(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(3,1)")]));
   model.set_mq2(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(3,2)")]));
   model.set_mq2(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("mq2(3,3)")]));

   model.set_ml2(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(1,1)")]));
   model.set_ml2(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(1,2)")]));
   model.set_ml2(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(1,3)")]));
   model.set_ml2(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(2,1)")]));
   model.set_ml2(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(2,2)")]));
   model.set_ml2(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(2,3)")]));
   model.set_ml2(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(3,1)")]));
   model.set_ml2(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(3,2)")]));
   model.set_ml2(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("ml2(3,3)")]));

   model.set_mHd2(boost::lexical_cast<double>(soft_parameters[columns.at("mHd2")]));
   model.set_mHu2(boost::lexical_cast<double>(soft_parameters[columns.at("mHu2")]));

   model.set_md2(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("md2(1,1)")]));
   model.set_md2(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("md2(1,2)")]));
   model.set_md2(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("md2(1,3)")]));
   model.set_md2(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("md2(2,1)")]));
   model.set_md2(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("md2(2,2)")]));
   model.set_md2(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("md2(2,3)")]));
   model.set_md2(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("md2(3,1)")]));
   model.set_md2(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("md2(3,2)")]));
   model.set_md2(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("md2(3,3)")]));

   model.set_mu2(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(1,1)")]));
   model.set_mu2(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(1,2)")]));
   model.set_mu2(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(1,3)")]));
   model.set_mu2(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(2,1)")]));
   model.set_mu2(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(2,2)")]));
   model.set_mu2(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(2,3)")]));
   model.set_mu2(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(3,1)")]));
   model.set_mu2(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(3,2)")]));
   model.set_mu2(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("mu2(3,3)")]));

   model.set_me2(0, 0, boost::lexical_cast<double>(soft_parameters[columns.at("me2(1,1)")]));
   model.set_me2(0, 1, boost::lexical_cast<double>(soft_parameters[columns.at("me2(1,2)")]));
   model.set_me2(0, 2, boost::lexical_cast<double>(soft_parameters[columns.at("me2(1,3)")]));
   model.set_me2(1, 0, boost::lexical_cast<double>(soft_parameters[columns.at("me2(2,1)")]));
   model.set_me2(1, 1, boost::lexical_cast<double>(soft_parameters[columns.at("me2(2,2)")]));
   model.set_me2(1, 2, boost::lexical_cast<double>(soft_parameters[columns.at("me2(2,3)")]));
   model.set_me2(2, 0, boost::lexical_cast<double>(soft_parameters[columns.at("me2(3,1)")]));
   model.set_me2(2, 1, boost::lexical_cast<double>(soft_parameters[columns.at("me2(3,2)")]));
   model.set_me2(2, 2, boost::lexical_cast<double>(soft_parameters[columns.at("me2(3,3)")]));

   model.set_MassB(boost::lexical_cast<double>(soft_parameters[columns.at("MassB")]));
   model.set_MassWB(boost::lexical_cast<double>(soft_parameters[columns.at("MassWB")]));
   model.set_MassG(boost::lexical_cast<double>(soft_parameters[columns.at("MassG")]));

   return model;
}

SUSYHD::SUSYHD_input_parameters initialize_MSSM(MSSM_mass_eigenstates model)
{
   SUSYHD::SUSYHD_input_parameters input;

   input.TanBeta = model.get_vu() / model.get_vd();
   input.M1 = model.get_MassB();
   input.M2 = model.get_MassWB();
   input.M3 = model.get_MassG();
   input.Mu = model.get_Mu();
   input.At = model.get_TYu(2,2) / model.get_Yu(2,2);
   input.mQ3 = Sqrt(model.get_mq2(2,2));
   input.mU3 = Sqrt(model.get_mu2(2,2));
   input.mD3 = Sqrt(model.get_md2(2,2));
   input.mQ2 = Sqrt(model.get_mq2(1,1));
   input.mU2 = Sqrt(model.get_mu2(1,1));
   input.mD2 = Sqrt(model.get_md2(1,1));
   input.mQ1 = Sqrt(model.get_mq2(0,0));
   input.mU1 = Sqrt(model.get_mu2(0,0));
   input.mD1 = Sqrt(model.get_md2(0,0));
   input.mL3 = Sqrt(model.get_ml2(2,2));
   input.mE3 = Sqrt(model.get_me2(2,2));
   input.mL2 = Sqrt(model.get_ml2(1,1));
   input.mE2 = Sqrt(model.get_me2(1,1));
   input.mL1 = Sqrt(model.get_ml2(0,0));
   input.mE1 = Sqrt(model.get_me2(0,0));
   const double Sin2Beta = Sin(2.0 * ArcTan(input.TanBeta));
   input.mA = Sqrt(2.0 * model.get_BMu() / Sin2Beta);

   return input;
}

void set_command_line_args(int argc, const char* argv[], std::string& susy_pars_file, std::string& soft_pars_file)
{
   for (int i = 1; i < argc; ++i) {
      const char* option = argv[i];

      if (Command_line_options::starts_with(option, "--susy-pars-file")) {
         susy_pars_file = Command_line_options::get_option_value(option, "=");
         continue;
      }

      if (Command_line_options::starts_with(option, "--soft-pars-file")) {
         soft_pars_file = Command_line_options::get_option_value(option, "=");
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

void print_result_header(int width)
{
   std::cout << "#" << ' '
             << std::left << std::setw(width) << "m12/GeV" << ' '
             << std::left << std::setw(width) << "Azero/GeV" << ' '
             << std::left << std::setw(width) << "TanBeta" << ' '
             << std::left << std::setw(width) << "MuInput/GeV" << ' '
             << std::left << std::setw(width) << "TanBeta(MS)" << ' '
             << std::left << std::setw(width) << "MassB(MS)/GeV" << ' '
             << std::left << std::setw(width) << "MassWB(MS)/GeV" << ' '
             << std::left << std::setw(width) << "MassG(MS)/GeV" << ' '
             << std::left << std::setw(width) << "MuEff(MS)/GeV" << ' '
             << std::left << std::setw(width) << "At(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mQ3(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mU3(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mD3(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mQ2(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mU2(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mD2(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mQ1(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mU1(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mD1(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mL3(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mE3(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mL2(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mE2(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mL1(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mE1(MS)/GeV" << ' '
             << std::left << std::setw(width) << "mA(MS)/GeV" << ' '
             << std::left << std::setw(width) << "Rscale/GeV" << ' '
             << std::left << std::setw(width) << "MhhEFT/GeV" << ' '
             << std::left << std::setw(width) << "error" << '\n';
}

void print_result_line(const CMSSM_semianalytic_input_parameters<Two_scale>& inputs, const SUSYHD::SUSYHD_input_parameters& higgs_mass_inputs, double Rscale, double Mhh, int width)
{
   const bool error = (Mhh < 1.0 ? true : false);

   std::cout << " " 
             << std::left << std::setw(width) << inputs.m12 << ' '
             << std::left << std::setw(width) << inputs.Azero << ' '
             << std::left << std::setw(width) << inputs.TanBeta << ' '
             << std::left << std::setw(width) << inputs.MuInput << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.TanBeta << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.M1 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.M2 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.M3 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.Mu << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.At << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mQ3 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mU3 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mD3 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mQ2 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mU2 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mD2 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mQ1 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mU1 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mD1 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mL3 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mE3 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mL2 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mE2 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mL1 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mE1 << ' '
             << std::left << std::setw(width) << higgs_mass_inputs.mA << ' '
             << std::left << std::setw(width) << Rscale << ' '
             << std::left << std::setw(width) << Mhh << ' '
             << std::left << std::setw(width) << error << ' ';
   if (Mhh < 1.0) {
      std::cout << "\t# invalid Higgs mass\n";
   } else {
      std::cout << '\n';
   }
}

} // namespace flexiblesusy

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;

   std::string susy_pars_file;
   std::string soft_pars_file;

   set_command_line_args(argc, argv, susy_pars_file, soft_pars_file);

   if (susy_pars_file.empty()) {
      ERROR("No SUSY parameters file given!\n"
            "   Please provide one via the option --susy-pars-file=");
      return EXIT_FAILURE;
   }

   if (soft_pars_file.empty()) {
      ERROR("No soft parameters file given!\n"
            "   Please provide one via the option --soft-pars-file=");
      return EXIT_FAILURE;
   }

   try {
      // open files
      std::ifstream susy_ifs(susy_pars_file, std::ifstream::in);
      if (susy_ifs.fail())
         throw ReadError("unable to open file " + susy_pars_file);
      std::ifstream soft_ifs(soft_pars_file, std::ifstream::in);
      if (soft_ifs.fail())
         throw ReadError("unable to open file " + soft_pars_file);

      std::map<std::string,std::size_t> columns;

      std::string susy_header_line;
      std::string soft_header_line;

      std::getline(susy_ifs, susy_header_line);
      std::getline(soft_ifs, soft_header_line);

      strip_comment_character(susy_header_line);
      strip_comment_character(soft_header_line);

      trim(susy_header_line);
      trim(soft_header_line);

      read_column_headers(susy_header_line, soft_header_line, columns);

      // open kernel link
      mathematica::MathematicaLink link("-linkname \"math -mathlink\"");
      SUSYHD::SUSYHDLink susyhd(&link);

      const int width = 18;
      print_result_header(width);

      std::string susy_pars_line;
      std::string soft_pars_line;
      while (std::getline(susy_ifs, susy_pars_line)
             && std::getline(soft_ifs, soft_pars_line)) {
         if (susy_pars_line.empty() || soft_pars_line.empty())
            continue;
         if (susy_pars_line.find("#") != std::string::npos)
            susy_pars_line = susy_pars_line.substr(0, susy_pars_line.find("#"));
         if (soft_pars_line.find("#") != std::string::npos)
            soft_pars_line = soft_pars_line.substr(0, soft_pars_line.find("#"));

         trim(susy_pars_line);
         trim(soft_pars_line);

         CMSSM_semianalytic_input_parameters<Two_scale> inputs
            = read_inputs_from_datapoint(susy_pars_line, columns);

         MSSM_mass_eigenstates model(
            initialize_model_from_datapoint(susy_pars_line, soft_pars_line, columns));

         model.calculate_DRbar_masses();

         const Eigen::Array<double,6,1> MSu(model.get_MSu());
         const Eigen::Matrix<double,6,6> ZU(model.get_ZU());
         const double Rscale = Sqrt(Power(MSu(0),Sqr(Abs(ZU(0,2)))
            + Sqr(Abs(ZU(0,5))))*Power(MSu(1),Sqr(Abs(ZU(1,2))) +
            Sqr(Abs(ZU(1,5))))*Power(MSu(2),Sqr(Abs(ZU(2,2))) +
            Sqr(Abs(ZU(2,5))))*Power(MSu(3),Sqr(Abs(ZU(3,2))) +
            Sqr(Abs(ZU(3,5))))*Power(MSu(4),Sqr(Abs(ZU(4,2))) +
            Sqr(Abs(ZU(4,5))))*Power(MSu(5),Sqr(Abs(ZU(5,2))) +
            Sqr(Abs(ZU(5,5)))));

         SUSYHD::SUSYHD_input_parameters higgs_mass_inputs
            = initialize_MSSM(model);

         susyhd.set_Rscale(Rscale);

         const double eft_mhh = susyhd.calculate_MHiggs(higgs_mass_inputs);

         print_result_line(inputs, higgs_mass_inputs, Rscale, eft_mhh, width);
      }
   } catch (const mathematica::Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
