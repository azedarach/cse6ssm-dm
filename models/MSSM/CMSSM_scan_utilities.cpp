#include "CMSSM_scan_utilities.hpp"
#include "CMSSM_two_scale_input_parameters.hpp"
#include "CMSSM_semi_two_scale_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace softsusy;

namespace flexiblesusy {

void write_CMSSM_inputs_list(std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CMSSM_inputs_list: "
            "file stream is corrupted");
      return;
   }

   filestr << std::left << std::setw(width) << "m0/GeV" << ' '
           << std::left << std::setw(width) << "m12/GeV" << ' '
           << std::left << std::setw(width) << "TanBeta" << ' '
           << std::left << std::setw(width) << "SignMu" << ' '
           << std::left << std::setw(width) << "Azero/GeV" << ' ';
}

void write_CMSSM_semianalytic_inputs_list(std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CMSSM_semianalytic_inputs_list: "
            "file stream is corrupted");
      return;
   }

   filestr << std::left << std::setw(width) << "m12/GeV" << ' '
           << std::left << std::setw(width) << "Azero/GeV" << ' '
           << std::left << std::setw(width) << "TanBeta" << ' '
           << std::left << std::setw(width) << "MuInput/GeV" << ' '
           << std::left << std::setw(width) << "MuInput_at_MS" << ' ';
}

void write_CMSSM_inputs(const CMSSM_input_parameters<Two_scale>& inputs, std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CMSSM_inputs: "
            "file stream is corrupted");
      return;
   }

   filestr << "  "
           << std::left << std::setw(width) << inputs.m0 << ' '
           << std::left << std::setw(width) << inputs.m12 << ' '
           << std::left << std::setw(width) << inputs.TanBeta << ' '
           << std::left << std::setw(width) << inputs.SignMu << ' '
           << std::left << std::setw(width) << inputs.Azero << ' ';
}

void write_CMSSM_inputs(const CMSSM_semianalytic_input_parameters<Two_scale>& inputs, std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CMSSM_inputs: "
            "file stream is corrupted");
      return;
   }

   filestr << "  "
           << std::left << std::setw(width) << inputs.m12 << ' '
           << std::left << std::setw(width) << inputs.Azero << ' '
           << std::left << std::setw(width) << inputs.TanBeta << ' '
           << std::left << std::setw(width) << inputs.MuInput << ' '
           << std::left << std::setw(width) << inputs.MuInput_at_MS << ' ';
}

CMSSM_pole_mass_writer::CMSSM_pole_mass_writer()
   : pole_masses()
   , pole_masses_inputs()
   , pole_masses_problems(MSSM_info::particle_names)
   , pole_masses_scale(0.0)
   , width(18)
{
}

CMSSM_semianalytic_pole_mass_writer::CMSSM_semianalytic_pole_mass_writer()
   : pole_masses()
   , pole_masses_inputs()
   , pole_masses_problems(MSSM_info::particle_names)
   , pole_masses_scale(0.0)
   , m0Sqr(0.0)
   , width(18)
{
}

void CMSSM_pole_mass_writer::write_pole_masses_comment_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_pole_mass_writer::write_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_pole_mass_writer::write_pole_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = pole_masses[p].name;
      const std::size_t multiplicity = pole_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_pole_mass_writer::write_pole_masses_comment_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_pole_mass_writer::write_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   // additional write-out for m0
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_pole_mass_writer::write_pole_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = pole_masses[p].name;
      const std::size_t multiplicity = pole_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_pole_mass_writer::write_pole_masses_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   write_CMSSM_inputs(pole_masses_inputs, filestr, width);

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_pole_mass_writer::write_pole_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = pole_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << pole_masses_problems.have_problem() << ' ';

   if (pole_masses_problems.have_problem() || pole_masses_problems.have_warning()) {
      filestr << "\t# " << pole_masses_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_semianalytic_pole_mass_writer::write_pole_masses_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   write_CMSSM_inputs(pole_masses_inputs, filestr, width);

   filestr << std::left << std::setw(width) << m0Sqr << ' ';

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_pole_mass_writer::write_pole_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = pole_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << pole_masses_problems.have_problem() << ' ';

   if (pole_masses_problems.have_problem() || pole_masses_problems.have_warning()) {
      filestr << "\t# " << pole_masses_problems << '\n';
   } else {
      filestr << '\n';
   }
}

CMSSM_drbar_values_writer::CMSSM_drbar_values_writer()
   : drbar_masses()
   , drbar_susy_pars()
   , drbar_soft_pars()
   , drbar_mixings()
   , drbar_masses_inputs()
   , drbar_susy_pars_inputs()
   , drbar_soft_pars_inputs()
   , drbar_mixings_inputs()
   , drbar_masses_problems(MSSM_info::particle_names)
   , drbar_susy_pars_problems(MSSM_info::particle_names)
   , drbar_soft_pars_problems(MSSM_info::particle_names)
   , drbar_mixings_problems(MSSM_info::particle_names)
   , drbar_masses_scale()
   , drbar_susy_pars_scale()
   , drbar_soft_pars_scale()
   , drbar_mixings_scale()
   , width(18)
{
}

CMSSM_semianalytic_drbar_values_writer::CMSSM_semianalytic_drbar_values_writer()
   : drbar_masses()
   , drbar_susy_pars()
   , drbar_soft_pars()
   , drbar_mixings()
   , drbar_masses_inputs()
   , drbar_susy_pars_inputs()
   , drbar_soft_pars_inputs()
   , drbar_mixings_inputs()
   , drbar_masses_problems(MSSM_info::particle_names)
   , drbar_susy_pars_problems(MSSM_info::particle_names)
   , drbar_soft_pars_problems(MSSM_info::particle_names)
   , drbar_mixings_problems(MSSM_info::particle_names)
   , drbar_masses_scale()
   , drbar_susy_pars_scale()
   , drbar_soft_pars_scale()
   , drbar_mixings_scale()
   , width(18)
{
}

void CMSSM_drbar_values_writer::write_drbar_masses_comment_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_drbar_values_writer::write_drbar_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_masses[p].name;
      const std::size_t multiplicity = drbar_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_masses_comment_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_masses[p].name;
      const std::size_t multiplicity = drbar_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_drbar_values_writer::write_drbar_susy_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_drbar_values_writer::write_drbar_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_susy_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_susy_pars[p].name;
      const std::size_t mass_dimension = drbar_susy_pars[p].mass_dimension;
      const std::size_t rows = drbar_susy_pars[p].rows;
      const std::size_t cols = drbar_susy_pars[p].cols;
      const std::size_t npars = rows * cols;

      for (std::size_t i = 0; i < npars; ++i) {
         std::string entry(name);
         if (npars > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_susy_pars[p].name;
      const std::size_t mass_dimension = drbar_susy_pars[p].mass_dimension;
      const std::size_t rows = drbar_susy_pars[p].rows;
      const std::size_t cols = drbar_susy_pars[p].cols;
      const std::size_t npars = rows * cols;

      for (std::size_t i = 0; i < npars; ++i) {
         std::string entry(name);
         if (npars > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_drbar_values_writer::write_drbar_soft_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_drbar_values_writer::write_drbar_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_soft_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_soft_pars[p].name;
      const std::size_t mass_dimension = drbar_soft_pars[p].mass_dimension;
      const std::size_t rows = drbar_soft_pars[p].rows;
      const std::size_t cols = drbar_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         std::string entry(name);
         if (num_entries > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_soft_pars[p].name;
      const std::size_t mass_dimension = drbar_soft_pars[p].mass_dimension;
      const std::size_t rows = drbar_soft_pars[p].rows;
      const std::size_t cols = drbar_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         std::string entry(name);
         if (num_entries > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_drbar_values_writer::write_drbar_mixings_comment_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_drbar_values_writer::write_drbar_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_mixings_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_mixings[p].name;
      const std::size_t dimension = drbar_mixings[p].dimension;
      const std::size_t num_entries = dimension * dimension;
      const bool is_real = drbar_mixings[p].is_real;

      if (is_real) {
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if (num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << entry << ' ';
         }
      } else {
         // first num_entries values are the real parts
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Re(" + entry + ")" << ' ';
         }
         // second num_entries values are the imaginary parts
         for (std::size_t i = num_entries; i < 2 * num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = (i - num_entries) / dimension;
               std::size_t row_index = i - num_entries - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Im(" + entry + ")" << ' ';
         }
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_mixings_comment_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_mixings_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = drbar_mixings[p].name;
      const std::size_t dimension = drbar_mixings[p].dimension;
      const std::size_t num_entries = dimension * dimension;
      const bool is_real = drbar_mixings[p].is_real;

      if (is_real) {
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if (num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << entry << ' ';
         }
      } else {
         // first num_entries values are the real parts
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Re(" + entry + ")" << ' ';
         }
         // second num_entries values are the imaginary parts
         for (std::size_t i = num_entries; i < 2 * num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = (i - num_entries) / dimension;
               std::size_t row_index = i - num_entries - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Im(" + entry + ")" << ' ';
         }
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_drbar_values_writer::write_drbar_masses_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   write_CMSSM_inputs(drbar_masses_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = drbar_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_masses_problems.have_problem() << ' ';

   if (drbar_masses_problems.have_problem() || drbar_masses_problems.have_warning()) {
      filestr << "\t# " << drbar_masses_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_masses_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   write_CMSSM_inputs(drbar_masses_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = drbar_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_masses_problems.have_problem() << ' ';

   if (drbar_masses_problems.have_problem() || drbar_masses_problems.have_warning()) {
      filestr << "\t# " << drbar_masses_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_drbar_values_writer::write_drbar_susy_pars_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   write_CMSSM_inputs(drbar_susy_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_susy_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = drbar_susy_pars[p].values;
      const std::size_t rows = drbar_susy_pars[p].rows;
      const std::size_t cols = drbar_susy_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_susy_pars_problems.have_problem() << ' ';

   if (drbar_susy_pars_problems.have_problem() || drbar_susy_pars_problems.have_warning()) {
      filestr << "\t# " << drbar_susy_pars_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   write_CMSSM_inputs(drbar_susy_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = drbar_susy_pars[p].values;
      const std::size_t rows = drbar_susy_pars[p].rows;
      const std::size_t cols = drbar_susy_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_susy_pars_problems.have_problem() << ' ';

   if (drbar_susy_pars_problems.have_problem() || drbar_susy_pars_problems.have_warning()) {
      filestr << "\t# " << drbar_susy_pars_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_drbar_values_writer::write_drbar_soft_pars_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   write_CMSSM_inputs(drbar_soft_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_soft_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = drbar_soft_pars[p].values;
      const std::size_t rows = drbar_soft_pars[p].rows;
      const std::size_t cols = drbar_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_soft_pars_problems.have_problem() << ' ';

   if (drbar_soft_pars_problems.have_problem() || drbar_soft_pars_problems.have_warning()) {
      filestr << "\t# " << drbar_soft_pars_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   write_CMSSM_inputs(drbar_soft_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = drbar_soft_pars[p].values;
      const std::size_t rows = drbar_soft_pars[p].rows;
      const std::size_t cols = drbar_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_soft_pars_problems.have_problem() << ' ';

   if (drbar_soft_pars_problems.have_problem() || drbar_soft_pars_problems.have_warning()) {
      filestr << "\t# " << drbar_soft_pars_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_drbar_values_writer::write_drbar_mixings_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   write_CMSSM_inputs(drbar_mixings_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_drbar_values_writer::write_drbar_mixings_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& mixings = drbar_mixings[p].mixings;
      const std::size_t num_mixings = mixings.size();

      for (std::size_t i = 0; i < num_mixings; ++i) {
         filestr << std::left << std::setw(width) << mixings[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_mixings_problems.have_problem() << ' ';

   if (drbar_mixings_problems.have_problem() || drbar_mixings_problems.have_warning()) {
      filestr << "\t# " << drbar_mixings_problems << '\n';
   } else {
      filestr << '\n';
   }
}

void CMSSM_semianalytic_drbar_values_writer::write_drbar_mixings_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   write_CMSSM_inputs(drbar_mixings_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_drbar_values_writer::write_drbar_mixings_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& mixings = drbar_mixings[p].mixings;
      const std::size_t num_mixings = mixings.size();

      for (std::size_t i = 0; i < num_mixings; ++i) {
         filestr << std::left << std::setw(width) << mixings[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << drbar_mixings_problems.have_problem() << ' ';

   if (drbar_mixings_problems.have_problem() || drbar_mixings_problems.have_warning()) {
      filestr << "\t# " << drbar_mixings_problems << '\n';
   } else {
      filestr << '\n';
   }
}

std::valarray<double> to_valarray(double v)
{
   return std::valarray<double>(&v, 1);
}

CMSSM_slha_values_writer::CMSSM_slha_values_writer()
   : slha_pole_masses()
   , slha_running_masses()
   , slha_susy_pars()
   , slha_soft_pars()
   , slha_pole_mixings()
   , slha_running_mixings()
   , slha_pole_masses_inputs()
   , slha_running_masses_inputs()
   , slha_susy_pars_inputs()
   , slha_soft_pars_inputs()
   , slha_pole_mixings_inputs()
   , slha_running_mixings_inputs()
   , slha_pole_masses_problems(MSSM_info::particle_names)
   , slha_running_masses_problems(MSSM_info::particle_names)
   , slha_susy_pars_problems(MSSM_info::particle_names)
   , slha_soft_pars_problems(MSSM_info::particle_names)
   , slha_pole_mixings_problems(MSSM_info::particle_names)
   , slha_running_mixings_problems(MSSM_info::particle_names)
   , high_scale(0)
   , susy_scale(0)
   , low_scale(0)
   , width(18)
{
}

CMSSM_semianalytic_slha_values_writer::CMSSM_semianalytic_slha_values_writer()
   : slha_pole_masses()
   , slha_running_masses()
   , slha_susy_pars()
   , slha_soft_pars()
   , slha_pole_mixings()
   , slha_running_mixings()
   , slha_pole_masses_inputs()
   , slha_running_masses_inputs()
   , slha_susy_pars_inputs()
   , slha_soft_pars_inputs()
   , slha_pole_mixings_inputs()
   , slha_running_mixings_inputs()
   , slha_pole_masses_problems(MSSM_info::particle_names)
   , slha_running_masses_problems(MSSM_info::particle_names)
   , slha_susy_pars_problems(MSSM_info::particle_names)
   , slha_soft_pars_problems(MSSM_info::particle_names)
   , slha_pole_mixings_problems(MSSM_info::particle_names)
   , slha_running_mixings_problems(MSSM_info::particle_names)
   , high_scale(0)
   , susy_scale(0)
   , low_scale(0)
   , slha_pole_masses_m0Sqr(0)
   , slha_running_masses_m0Sqr(0)
   , slha_susy_pars_m0Sqr(0)
   , slha_soft_pars_m0Sqr(0)
   , slha_pole_mixings_m0Sqr(0)
   , slha_running_mixings_m0Sqr(0)
   , width(18)
{
}

void CMSSM_slha_values_writer::write_slha_pole_masses_comment_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_pole_masses[p].name;
      const std::size_t multiplicity = slha_pole_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_pole_masses[p].name;
      const std::size_t multiplicity = slha_pole_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_slha_values_writer::write_slha_running_masses_comment_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_running_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_running_masses[p].name;
      const std::size_t multiplicity = slha_running_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_running_masses[p].name;
      const std::size_t multiplicity = slha_running_masses[p].masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string entry(name);
         if (multiplicity > 1) {
            std::ostringstream index;
            index << "(" << i + 1 << ")";
            entry.append(index.str());
         }
         entry.append("/GeV");
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_slha_values_writer::write_slha_susy_pars_comment_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_susy_pars[p].name;
      const std::size_t mass_dimension = slha_susy_pars[p].mass_dimension;
      const std::size_t rows = slha_susy_pars[p].rows;
      const std::size_t cols = slha_susy_pars[p].cols;
      const std::size_t npars = rows * cols;

      for (std::size_t i = 0; i < npars; ++i) {
         std::string entry(name);
         if (npars > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_susy_pars[p].name;
      const std::size_t mass_dimension = slha_susy_pars[p].mass_dimension;
      const std::size_t rows = slha_susy_pars[p].rows;
      const std::size_t cols = slha_susy_pars[p].cols;
      const std::size_t npars = rows * cols;

      for (std::size_t i = 0; i < npars; ++i) {
         std::string entry(name);
         if (npars > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_slha_values_writer::write_slha_soft_pars_comment_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_soft_pars[p].name;
      const std::size_t mass_dimension = slha_soft_pars[p].mass_dimension;
      const std::size_t rows = slha_soft_pars[p].rows;
      const std::size_t cols = slha_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         std::string entry(name);
         if (num_entries > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_soft_pars[p].name;
      const std::size_t mass_dimension = slha_soft_pars[p].mass_dimension;
      const std::size_t rows = slha_soft_pars[p].rows;
      const std::size_t cols = slha_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         std::string entry(name);
         if (num_entries > 1) {
            std::size_t col_index = i / rows;
            std::size_t row_index = i - rows * col_index;
            std::ostringstream index;
            index << "(" << row_index + 1 << "," << col_index + 1 << ")";
            entry.append(index.str());
         }
         if (mass_dimension > 0) {
            entry.append("/GeV");
            if (mass_dimension > 1) {
               std::ostringstream dim;
               dim << "^" << mass_dimension;
               entry.append(dim.str());
            }
         }
         filestr << std::left << std::setw(width) << entry << ' ';
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_slha_values_writer::write_slha_pole_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_pole_mixings[p].name;
      const std::size_t dimension = slha_pole_mixings[p].dimension;
      const std::size_t num_entries = dimension * dimension;
      const bool is_real = slha_pole_mixings[p].is_real;

      if (is_real) {
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if (num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << entry << ' ';
         }
      } else {
         // first num_entries values are the real parts
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Re(" + entry + ")" << ' ';
         }
         // second num_entries values are the imaginary parts
         for (std::size_t i = num_entries; i < 2 * num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = (i - num_entries) / dimension;
               std::size_t row_index = i - num_entries - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Im(" + entry + ")" << ' ';
         }
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';


   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_pole_mixings[p].name;
      const std::size_t dimension = slha_pole_mixings[p].dimension;
      const std::size_t num_entries = dimension * dimension;
      const bool is_real = slha_pole_mixings[p].is_real;

      if (is_real) {
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if (num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << entry << ' ';
         }
      } else {
         // first num_entries values are the real parts
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Re(" + entry + ")" << ' ';
         }
         // second num_entries values are the imaginary parts
         for (std::size_t i = num_entries; i < 2 * num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = (i - num_entries) / dimension;
               std::size_t row_index = i - num_entries - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Im(" + entry + ")" << ' ';
         }
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_slha_values_writer::write_slha_running_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_running_mixings[p].name;
      const std::size_t dimension = slha_running_mixings[p].dimension;
      const std::size_t num_entries = dimension * dimension;
      const bool is_real = slha_running_mixings[p].is_real;

      if (is_real) {
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if (num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << entry << ' ';
         }
      } else {
         // first num_entries values are the real parts
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Re(" + entry + ")" << ' ';
         }
         // second num_entries values are the imaginary parts
         for (std::size_t i = num_entries; i < 2 * num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = (i - num_entries) / dimension;
               std::size_t row_index = i - num_entries - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Im(" + entry + ")" << ' ';
         }
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CMSSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV^2" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';


   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = slha_running_mixings[p].name;
      const std::size_t dimension = slha_running_mixings[p].dimension;
      const std::size_t num_entries = dimension * dimension;
      const bool is_real = slha_running_mixings[p].is_real;

      if (is_real) {
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if (num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << entry << ' ';
         }
      } else {
         // first num_entries values are the real parts
         for (std::size_t i = 0; i < num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = i / dimension;
               std::size_t row_index = i - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Re(" + entry + ")" << ' ';
         }
         // second num_entries values are the imaginary parts
         for (std::size_t i = num_entries; i < 2 * num_entries; ++i) {
            std::string entry(name);
            if(num_entries > 1) {
               std::size_t col_index = (i - num_entries) / dimension;
               std::size_t row_index = i - num_entries - dimension * col_index;
               std::ostringstream index;
               index << "(" << row_index + 1 << "," << col_index + 1 << ")";
               entry.append(index.str());
            }
            filestr << std::left << std::setw(width) << "Im(" + entry + ")" << ' ';
         }
      }
   }

   filestr << std::left << std::setw(width) << "error" << '\n';
}

void CMSSM_slha_values_writer::write_slha_pole_masses_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_pole_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_pole_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = slha_pole_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_pole_masses_problems.have_problem() << ' ';

   if (slha_pole_masses_problems.have_problem() || slha_pole_masses_problems.have_warning()) {
      filestr << "\t# " << slha_pole_masses_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_pole_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_pole_masses_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = slha_pole_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_pole_masses_problems.have_problem() << ' ';

   if (slha_pole_masses_problems.have_problem() || slha_pole_masses_problems.have_warning()) {
      filestr << "\t# " << slha_pole_masses_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_slha_values_writer::write_slha_running_masses_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_running_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_running_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = slha_running_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_running_masses_problems.have_problem() << ' ';

   if (slha_running_masses_problems.have_problem() || slha_running_masses_problems.have_warning()) {
      filestr << "\t# " << slha_running_masses_problems << '\n';
   } else {
      filestr << '\n';
   }
   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_running_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_running_masses_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& masses = slha_running_masses[p].masses;
      const std::size_t multiplicity = masses.size();

      for (std::size_t i = 0; i < multiplicity; ++i) {
         filestr << std::left << std::setw(width) << masses[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_running_masses_problems.have_problem() << ' ';

   if (slha_running_masses_problems.have_problem() || slha_running_masses_problems.have_warning()) {
      filestr << "\t# " << slha_running_masses_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_slha_values_writer::write_slha_susy_pars_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_susy_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_susy_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = slha_susy_pars[p].values;
      const std::size_t rows = slha_susy_pars[p].rows;
      const std::size_t cols = slha_susy_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_susy_pars_problems.have_problem() << ' ';

   if (slha_susy_pars_problems.have_problem() || slha_susy_pars_problems.have_warning()) {
      filestr << "\t# " << slha_susy_pars_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_susy_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_susy_pars_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = slha_susy_pars[p].values;
      const std::size_t rows = slha_susy_pars[p].rows;
      const std::size_t cols = slha_susy_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_susy_pars_problems.have_problem() << ' ';

   if (slha_susy_pars_problems.have_problem() || slha_susy_pars_problems.have_warning()) {
      filestr << "\t# " << slha_susy_pars_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_slha_values_writer::write_slha_soft_pars_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_soft_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_soft_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = slha_soft_pars[p].values;
      const std::size_t rows = slha_soft_pars[p].rows;
      const std::size_t cols = slha_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_soft_pars_problems.have_problem() << ' ';

   if (slha_soft_pars_problems.have_problem() || slha_soft_pars_problems.have_warning()) {
      filestr << "\t# " << slha_soft_pars_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_soft_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_soft_pars_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& values = slha_soft_pars[p].values;
      const std::size_t rows = slha_soft_pars[p].rows;
      const std::size_t cols = slha_soft_pars[p].cols;
      const std::size_t num_entries = rows * cols;

      for (std::size_t i = 0; i < num_entries; ++i) {
         filestr << std::left << std::setw(width) << values[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_soft_pars_problems.have_problem() << ' ';

   if (slha_soft_pars_problems.have_problem() || slha_soft_pars_problems.have_warning()) {
      filestr << "\t# " << slha_soft_pars_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_slha_values_writer::write_slha_pole_mixings_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_pole_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_pole_mixings_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& mixings = slha_pole_mixings[p].mixings;
      const std::size_t num_mixings = mixings.size();

      for (std::size_t i = 0; i < num_mixings; ++i) {
         filestr << std::left << std::setw(width) << mixings[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_pole_mixings_problems.have_problem() << ' ';

   if (slha_pole_mixings_problems.have_problem() || slha_pole_mixings_problems.have_warning()) {
      filestr << "\t# " << slha_pole_mixings_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_pole_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_pole_mixings_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& mixings = slha_pole_mixings[p].mixings;
      const std::size_t num_mixings = mixings.size();

      for (std::size_t i = 0; i < num_mixings; ++i) {
         filestr << std::left << std::setw(width) << mixings[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_pole_mixings_problems.have_problem() << ' ';

   if (slha_pole_mixings_problems.have_problem() || slha_pole_mixings_problems.have_warning()) {
      filestr << "\t# " << slha_pole_mixings_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_slha_values_writer::write_slha_running_mixings_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_running_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_slha_values_writer::write_slha_running_mixings_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& mixings = slha_running_mixings[p].mixings;
      const std::size_t num_mixings = mixings.size();

      for (std::size_t i = 0; i < num_mixings; ++i) {
         filestr << std::left << std::setw(width) << mixings[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_running_mixings_problems.have_problem() << ' ';

   if (slha_running_mixings_problems.have_problem() || slha_running_mixings_problems.have_warning()) {
      filestr << "\t# " << slha_running_mixings_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

void CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CMSSM_inputs(slha_running_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_running_mixings_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CMSSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& mixings = slha_running_mixings[p].mixings;
      const std::size_t num_mixings = mixings.size();

      for (std::size_t i = 0; i < num_mixings; ++i) {
         filestr << std::left << std::setw(width) << mixings[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << slha_running_mixings_problems.have_problem() << ' ';

   if (slha_running_mixings_problems.have_problem() || slha_running_mixings_problems.have_warning()) {
      filestr << "\t# " << slha_running_mixings_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

} // namespace flexiblesusy

