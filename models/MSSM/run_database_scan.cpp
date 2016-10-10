#include "CMSSM_semi_two_scale_input_parameters.hpp"
#include "CMSSM_semianalytic_spectrum_generator.hpp"
#include "CMSSM_scan_utilities.hpp"

#include "error.hpp"
#include "lowe.h"
#include "grid_scanner.hpp"
#include "scan_command_line_options.hpp"

#include <iostream>
#include <string>
#include <vector>

#ifdef ENABLE_SQLITE
#include <sqlite3.h>
#include <boost/lexical_cast.hpp>
#endif

namespace flexiblesusy {

int extract_callback(void* data, int argc, char** argv, char** col_name)
{
   std::vector<double>* values = static_cast<std::vector<double>*>(data);

   for (int i = 0; i < argc; i++) {
      values->push_back(boost::lexical_cast<double>(argv[i]));
   }

   return 0;
}

typedef int (*TCallback)(void*, int, char**, char**);

void execute(sqlite3* db, const std::string& cmd, TCallback callback, void* data)
{
   if (!db)
      return;

   char* zErrMsg = 0;
   const int rc = sqlite3_exec(db, cmd.c_str(), callback, data, &zErrMsg);

   if (rc != SQLITE_OK) {
      ERROR("SQL error while executing command \"" << cmd << "\": " << zErrMsg);
      sqlite3_free(zErrMsg);
   } else {
      VERBOSE_MSG("SQL command \"" << cmd << "\" executed successfully");
   }
}

void read_input_parameter_values(const std::string & scan_input_file,
                                 const std::string& table_name,
                                 std::vector<double>& m12_values,
                                 std::vector<double>& Azero_values,
                                 std::vector<double>& TanBeta_values,
                                 std::vector<double>& MuInput_values)
{
#ifndef ENABLE_SQLITE
   ERROR("SQLite support is not supported!\n");
   exit(EXIT_FAILURE);
#else
   sqlite3* db = 0;

   const int rc = sqlite3_open(scan_input_file.c_str(), &db);

   if (rc) {
      ERROR("Can't open database: " << sqlite3_errmsg(db));
      exit(EXIT_FAILURE);
   }

   std::string sql;
   sql = "SELECT m12 FROM " + table_name + ";";
   execute(db, sql, extract_callback, static_cast<void*>(&m12_values));

   sql = "SELECT Azero FROM " + table_name + ";";
   execute(db, sql, extract_callback, static_cast<void*>(&Azero_values));

   sql = "SELECT TanBeta FROM " + table_name + ";";
   execute(db, sql, extract_callback, static_cast<void*>(&TanBeta_values));

   sql = "SELECT MuInput FROM " + table_name + ";";
   execute(db, sql, extract_callback, static_cast<void*>(&MuInput_values));

   if (db)
      sqlite3_close(db);
#endif
}

void write_inputs_comment_line(std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("write_inputs_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# "
           << std::left << std::setw(width) << "m12/GeV" << ' '
           << std::left << std::setw(width) << "Azero/GeV" << ' '
           << std::left << std::setw(width) << "TanBeta/GeV" << ' '
           << std::left << std::setw(width) << "MuInput/GeV" << ' '
           << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' '
           << std::left << std::setw(width) << "AlphaSMZ" << ' '
           << std::left << std::setw(width) << "MtPole/GeV" << ' '
           << std::left << std::setw(width) << "error" << '\n';
}

void write_inputs_line(std::ostream & filestr,
                       double m12, double Azero, double TanBeta,
                       double MuInput, double m0Sqr,
                       double alphas, double mtp,
                       const Problems<MSSM_info::NUMBER_OF_PARTICLES>& problems,
                       std::size_t width)
{
   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   if (!filestr.good()) {
      ERROR("write_inputs_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "  "
           << std::left << std::setw(width) << m12 << ' '
           << std::left << std::setw(width) << Azero << ' '
           << std::left << std::setw(width) << TanBeta << ' '
           << std::left << std::setw(width) << MuInput << ' '
           << std::left << std::setw(width) << m0Sqr << ' '
           << std::left << std::setw(width) << alphas << ' '
           << std::left << std::setw(width) << mtp << ' '
           << std::left << std::setw(width)
           << problems.have_problem() << ' ';

   if (problems.have_problem() || problems.have_warning()) {
      filestr << "\t# " << problems << '\n';
   } else {
      filestr << '\n';
   }
   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

} // namespace flexiblesusy

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   typedef Two_scale algorithm_type;

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

   if (scan_input_file.empty()) {
      ERROR("No scan input database provided!\n"
            "   Please provide one via the option --scan-input-file=");
      return EXIT_FAILURE;
   }

   // output streams
   bool must_write_pole_masses = false;
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

   std::ofstream pole_mass_out_stream;
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

   if (!pole_mass_output_file.empty()) {
      must_write_pole_masses = true;
      pole_mass_out_stream.open(pole_mass_output_file, std::ofstream::out);
   }
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

   const std::string& table_name("slha_conventions");
   std::vector<double> m12_values;
   std::vector<double> Azero_values;
   std::vector<double> TanBeta_values;
   std::vector<double> MuInput_values;

   read_input_parameter_values(scan_input_file, table_name,
                               m12_values, Azero_values,
                               TanBeta_values, MuInput_values);

   const std::size_t n_seeds = m12_values.size();
   const std::size_t n_pts = 5;

   const double alphas_cent = 0.1181;
   const double alphas_uncert = 0.0013;
   const double mtp_cent = 173.21;
   const double mtp_uncert = 0.87;

   const double delta_alphas = 4.0 * alphas_uncert / (n_pts - 1.0);
   const double delta_mtp = 4.0 * mtp_uncert / (n_pts - 1.0);

   CMSSM_semianalytic_pole_mass_writer pole_mass_writer;
   CMSSM_semianalytic_drbar_values_writer drbar_values_writer;
   CMSSM_semianalytic_slha_values_writer slha_values_writer;
   bool must_write_comment_line = true;
   for (std::size_t i = 0; i < n_seeds; ++i) {
      CMSSM_semianalytic_input_parameters<algorithm_type> input;
      input.m12 = m12_values[i];
      input.TanBeta = TanBeta_values[i];
      input.Azero = Azero_values[i];
      input.MuInput = MuInput_values[i];

      for (std::size_t j = 0; j < n_pts; ++j) {
         const double alphas = alphas_cent - 2.0 * alphas_uncert
            + j * delta_alphas;
         for (std::size_t k = 0; k < n_pts; ++k) {
            const double mtp = mtp_cent - 2.0 * mtp_uncert
               + k * delta_mtp;

            softsusy::QedQcd oneset;
            oneset.setAlpha(ALPHAS, alphas);
            oneset.setPoleMt(mtp);
            oneset.toMz();

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
            const Problems<MSSM_info::NUMBER_OF_PARTICLES> problems
               = model.get_problems();

            slha_values_writer.set_high_scale(spectrum_generator.get_high_scale());
            slha_values_writer.set_susy_scale(spectrum_generator.get_susy_scale());
            slha_values_writer.set_low_scale(spectrum_generator.get_low_scale());

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

               write_inputs_comment_line(std::cout, 18);

               if (must_write_pole_masses)
                  pole_mass_writer.write_pole_masses_comment_line(pole_mass_out_stream);
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

            write_inputs_line(std::cout, m12_values[i],
                              Azero_values[i],
                              TanBeta_values[i], MuInput_values[i],
                              model.get_ewsb_output_parameter(0),
                              alphas, mtp, problems, 18);

            if (must_write_pole_masses)
               pole_mass_writer.write_pole_masses_line(pole_mass_out_stream);
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

         }
      }
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

   return 0;
}
