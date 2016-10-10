#include "CSE6SSM_scan_utilities.hpp"
#include "CSE6SSM_two_scale_input_parameters.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace softsusy;

namespace flexiblesusy {

void write_CSE6SSM_inputs_list(std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CSE6SSM_inputs_list: "
            "file stream is corrupted");
      return;
   }

   filestr << std::left << std::setw(width) << "m0/GeV" << ' '
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
           << std::left << std::setw(width) << "BMuPhiInput/GeV^2" << ' ';
}

void write_CSE6SSM_semianalytic_inputs_list(std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CSE6SSM_semianalytic_inputs_list: "
            "file stream is corrupted");
      return;
   }

   filestr << std::left << std::setw(width) << "m12/GeV" << ' '
           << std::left << std::setw(width) << "Azero/GeV" << ' '
           << std::left << std::setw(width) << "TanBeta" << ' '
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
           << std::left << std::setw(width) << "LambdaxInput" << ' '
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
           << std::left << std::setw(width) << "BMuPhiInput/GeV^2" << ' ';
}

void write_CSE6SSM_inputs(const CSE6SSM_input_parameters<Two_scale>& inputs, std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CSE6SSM_inputs: "
            "file stream is corrupted");
      return;
   }

   filestr << "  "
           << std::left << std::setw(width) << inputs.m0 << ' '
           << std::left << std::setw(width) << inputs.m12 << ' '
           << std::left << std::setw(width) << inputs.TanBeta << ' '
           << std::left << std::setw(width) << inputs.SignLambdax << ' '
           << std::left << std::setw(width) << inputs.Azero << ' '
           << std::left << std::setw(width) << inputs.sInput << ' '
           << std::left << std::setw(width) << inputs.QSInput << ' '
           << std::left << std::setw(width) << inputs.hEInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.hEInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.hEInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.hEInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.hEInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.hEInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.SigmaLInput << ' '
           << std::left << std::setw(width) << inputs.KappaPrInput << ' '
           << std::left << std::setw(width) << inputs.SigmaxInput << ' '
           << std::left << std::setw(width) << inputs.gDInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.gDInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.gDInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.gDInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.gDInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.gDInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.gDInput(0,2) << ' '
           << std::left << std::setw(width) << inputs.gDInput(1,2) << ' '
           << std::left << std::setw(width) << inputs.gDInput(2,2) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(0,2) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(1,2) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(2,2) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(0,0) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(1,0) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(0,1) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(1,1) << ' '
           << std::left << std::setw(width) << inputs.fuInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.fuInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.fuInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.fuInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.fuInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.fuInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.fdInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.fdInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.fdInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.fdInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.fdInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.fdInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.MuPrInput << ' '
           << std::left << std::setw(width) << inputs.MuPhiInput << ' '
           << std::left << std::setw(width) << inputs.BMuPrInput << ' '
           << std::left << std::setw(width) << inputs.BMuPhiInput << ' ';
}

void write_CSE6SSM_inputs(const CSE6SSM_semianalytic_input_parameters<Two_scale>& inputs, std::ostream & filestr, std::size_t width)
{
   if (!filestr.good()) {
      ERROR("flexiblesusy::write_CSE6SSM_inputs: "
            "file stream is corrupted");
      return;
   }

   filestr << "  "
           << std::left << std::setw(width) << inputs.m12 << ' '
           << std::left << std::setw(width) << inputs.Azero << ' '
           << std::left << std::setw(width) << inputs.TanBeta << ' '
           << std::left << std::setw(width) << inputs.sInput << ' '
           << std::left << std::setw(width) << inputs.QSInput << ' '
           << std::left << std::setw(width) << inputs.hEInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.hEInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.hEInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.hEInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.hEInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.hEInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.SigmaLInput << ' '
           << std::left << std::setw(width) << inputs.KappaPrInput << ' '
           << std::left << std::setw(width) << inputs.SigmaxInput << ' '
           << std::left << std::setw(width) << inputs.gDInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.gDInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.gDInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.gDInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.gDInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.gDInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.gDInput(0,2) << ' '
           << std::left << std::setw(width) << inputs.gDInput(1,2) << ' '
           << std::left << std::setw(width) << inputs.gDInput(2,2) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(0,2) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(1,2) << ' '
           << std::left << std::setw(width) << inputs.KappaInput(2,2) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(0,0) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(1,0) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(0,1) << ' '
           << std::left << std::setw(width) << inputs.Lambda12Input(1,1) << ' '
           << std::left << std::setw(width) << inputs.LambdaxInput << ' '
           << std::left << std::setw(width) << inputs.fuInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.fuInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.fuInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.fuInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.fuInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.fuInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.fdInput(0,0) << ' '
           << std::left << std::setw(width) << inputs.fdInput(1,0) << ' '
           << std::left << std::setw(width) << inputs.fdInput(2,0) << ' '
           << std::left << std::setw(width) << inputs.fdInput(0,1) << ' '
           << std::left << std::setw(width) << inputs.fdInput(1,1) << ' '
           << std::left << std::setw(width) << inputs.fdInput(2,1) << ' '
           << std::left << std::setw(width) << inputs.MuPrInput << ' '
           << std::left << std::setw(width) << inputs.MuPhiInput << ' '
           << std::left << std::setw(width) << inputs.BMuPrInput << ' '
           << std::left << std::setw(width) << inputs.BMuPhiInput << ' ';
}

CSE6SSM_pole_mass_writer::CSE6SSM_pole_mass_writer()
   : pole_masses()
   , pole_masses_inputs()
   , pole_masses_problems(CSE6SSM_info::particle_names)
   , pole_masses_scale(0.0)
   , width(18)
{
}

CSE6SSM_semianalytic_pole_mass_writer::CSE6SSM_semianalytic_pole_mass_writer()
   : pole_masses()
   , pole_masses_inputs()
   , pole_masses_problems(CSE6SSM_info::particle_names)
   , pole_masses_scale(0.0)
   , m0Sqr(0.0)
   , width(18)
{
}

void CSE6SSM_pole_mass_writer::write_pole_masses_comment_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_pole_mass_writer::write_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_pole_mass_writer::write_pole_masses_comment_line: "
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

void CSE6SSM_semianalytic_pole_mass_writer::write_pole_masses_comment_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_pole_mass_writer::write_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   // additional write-out for m0Sqr
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_pole_mass_writer::write_pole_masses_comment_line: "
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

void CSE6SSM_pole_mass_writer::write_pole_masses_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   write_CSE6SSM_inputs(pole_masses_inputs, filestr, width);

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_pole_mass_writer::write_pole_masses_line: "
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

void CSE6SSM_semianalytic_pole_mass_writer::write_pole_masses_line(std::ostream & filestr) const
{
   if (pole_masses.empty())
      return;

   write_CSE6SSM_inputs(pole_masses_inputs, filestr, width);

   filestr << std::left << std::setw(width) << m0Sqr << ' ';

   for (std::size_t p = 0; p < pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_pole_mass_writer::write_pole_masses_line: "
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

CSE6SSM_drbar_values_writer::CSE6SSM_drbar_values_writer()
   : drbar_masses()
   , drbar_susy_pars()
   , drbar_soft_pars()
   , drbar_mixings()
   , drbar_masses_inputs()
   , drbar_susy_pars_inputs()
   , drbar_soft_pars_inputs()
   , drbar_mixings_inputs()
   , drbar_masses_problems(CSE6SSM_info::particle_names)
   , drbar_susy_pars_problems(CSE6SSM_info::particle_names)
   , drbar_soft_pars_problems(CSE6SSM_info::particle_names)
   , drbar_mixings_problems(CSE6SSM_info::particle_names)
   , drbar_masses_scale()
   , drbar_susy_pars_scale()
   , drbar_soft_pars_scale()
   , drbar_mixings_scale()
   , width(18)
{
}

CSE6SSM_semianalytic_drbar_values_writer::CSE6SSM_semianalytic_drbar_values_writer()
   : drbar_masses()
   , drbar_susy_pars()
   , drbar_soft_pars()
   , drbar_mixings()
   , drbar_masses_inputs()
   , drbar_susy_pars_inputs()
   , drbar_soft_pars_inputs()
   , drbar_mixings_inputs()
   , drbar_masses_problems(CSE6SSM_info::particle_names)
   , drbar_susy_pars_problems(CSE6SSM_info::particle_names)
   , drbar_soft_pars_problems(CSE6SSM_info::particle_names)
   , drbar_mixings_problems(CSE6SSM_info::particle_names)
   , drbar_masses_scale()
   , drbar_susy_pars_scale()
   , drbar_soft_pars_scale()
   , drbar_mixings_scale()
   , width(18)
{
}

void CSE6SSM_drbar_values_writer::write_drbar_masses_comment_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_drbar_values_writer::write_drbar_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_masses_comment_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_masses_comment_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_masses_comment_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_susy_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_drbar_values_writer::write_drbar_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_susy_pars_comment_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_comment_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_soft_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_drbar_values_writer::write_drbar_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_soft_pars_comment_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_comment_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_comment_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_mixings_comment_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_drbar_values_writer::write_drbar_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_mixings_comment_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_mixings_comment_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_mixings_comment_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_masses_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   write_CSE6SSM_inputs(drbar_masses_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_masses_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_masses_line(std::ostream & filestr) const
{
   if (drbar_masses.empty())
      return;

   write_CSE6SSM_inputs(drbar_masses_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_masses_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_susy_pars_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   write_CSE6SSM_inputs(drbar_susy_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_susy_pars_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_line(std::ostream & filestr) const
{
   if (drbar_susy_pars.empty())
      return;

   write_CSE6SSM_inputs(drbar_susy_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_susy_pars_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_soft_pars_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   write_CSE6SSM_inputs(drbar_soft_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_soft_pars_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_line(std::ostream & filestr) const
{
   if (drbar_soft_pars.empty())
      return;

   write_CSE6SSM_inputs(drbar_soft_pars_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_soft_pars_line: "
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

void CSE6SSM_drbar_values_writer::write_drbar_mixings_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   write_CSE6SSM_inputs(drbar_mixings_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_drbar_values_writer::write_drbar_mixings_line: "
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

void CSE6SSM_semianalytic_drbar_values_writer::write_drbar_mixings_line(std::ostream & filestr) const
{
   if (drbar_mixings.empty())
      return;

   write_CSE6SSM_inputs(drbar_mixings_inputs, filestr, width);

   for (std::size_t p = 0; p < drbar_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_drbar_values_writer::write_drbar_mixings_line: "
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

CSE6SSM_slha_values_writer::CSE6SSM_slha_values_writer()
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
   , slha_pole_masses_problems(CSE6SSM_info::particle_names)
   , slha_running_masses_problems(CSE6SSM_info::particle_names)
   , slha_susy_pars_problems(CSE6SSM_info::particle_names)
   , slha_soft_pars_problems(CSE6SSM_info::particle_names)
   , slha_pole_mixings_problems(CSE6SSM_info::particle_names)
   , slha_running_mixings_problems(CSE6SSM_info::particle_names)
   , high_scale(0)
   , susy_scale(0)
   , low_scale(0)
   , width(18)
{
}

CSE6SSM_semianalytic_slha_values_writer::CSE6SSM_semianalytic_slha_values_writer()
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
   , slha_pole_masses_problems(CSE6SSM_info::particle_names)
   , slha_running_masses_problems(CSE6SSM_info::particle_names)
   , slha_susy_pars_problems(CSE6SSM_info::particle_names)
   , slha_soft_pars_problems(CSE6SSM_info::particle_names)
   , slha_pole_mixings_problems(CSE6SSM_info::particle_names)
   , slha_running_mixings_problems(CSE6SSM_info::particle_names)
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

void CSE6SSM_slha_values_writer::write_slha_pole_masses_comment_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_comment_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_comment_line: "
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

void CSE6SSM_slha_values_writer::write_slha_running_masses_comment_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_comment_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_comment_line: "
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

void CSE6SSM_slha_values_writer::write_slha_susy_pars_comment_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_comment_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_comment_line: "
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

void CSE6SSM_slha_values_writer::write_slha_soft_pars_comment_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_comment_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_comment_line: "
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

void CSE6SSM_slha_values_writer::write_slha_pole_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_comment_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';


   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_comment_line: "
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

void CSE6SSM_slha_values_writer::write_slha_running_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_comment_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';


   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_comment_line: "
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

void CSE6SSM_slha_values_writer::write_slha_pole_masses_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_pole_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_pole_masses_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_line(std::ostream & filestr) const
{
   if (slha_pole_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_pole_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_pole_masses_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_masses_line: "
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

void CSE6SSM_slha_values_writer::write_slha_running_masses_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_running_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_running_masses_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_line(std::ostream & filestr) const
{
   if (slha_running_masses.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_running_masses_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_running_masses_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_masses.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_masses_line: "
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

void CSE6SSM_slha_values_writer::write_slha_susy_pars_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_susy_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_susy_pars_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_line(std::ostream & filestr) const
{
   if (slha_susy_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_susy_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_susy_pars_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_susy_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_susy_pars_line: "
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

void CSE6SSM_slha_values_writer::write_slha_soft_pars_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_soft_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_soft_pars_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_line(std::ostream & filestr) const
{
   if (slha_soft_pars.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_soft_pars_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_soft_pars_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_soft_pars.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_soft_pars_line: "
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

void CSE6SSM_slha_values_writer::write_slha_pole_mixings_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_pole_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_pole_mixings_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line(std::ostream & filestr) const
{
   if (slha_pole_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_pole_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_pole_mixings_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_pole_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_pole_mixings_line: "
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

void CSE6SSM_slha_values_writer::write_slha_running_mixings_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_running_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_slha_values_writer::write_slha_running_mixings_line: "
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

void CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_line(std::ostream & filestr) const
{
   if (slha_running_mixings.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(slha_running_mixings_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << slha_running_mixings_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < slha_running_mixings.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_slha_values_writer::write_slha_running_mixings_line: "
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

CSE6SSM_semianalytic_coefficients_writer::CSE6SSM_semianalytic_coefficients_writer()
   : coefficients()
   , coefficients_inputs()
   , coefficients_problems(CSE6SSM_info::particle_names)
   , high_scale(0)
   , susy_scale(0)
   , low_scale(0)
   , width(18)
{
}

void CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line(std::ostream & filestr) const
{
   if (coefficients.empty())
      return;

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line: "
            "file stream is corrupted");
      return;
   }

   filestr << "# ";

   write_CSE6SSM_semianalytic_inputs_list(filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "m0Sqr/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "LowScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "SUSYScale/GeV" << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << "HighScale/GeV" << ' ';

   for (std::size_t p = 0; p < coefficients.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_comment_line: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = coefficients[p].name;
      const std::size_t mass_dimension = coefficients[p].mass_dimension;
      const std::size_t rows = coefficients[p].rows;
      const std::size_t cols = coefficients[p].cols;
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

void CSE6SSM_semianalytic_coefficients_writer::write_coefficients_line(std::ostream & filestr) const
{
   if (coefficients.empty())
      return;

   // ensure SLHA formatted output
   std::ios_base::fmtflags old_flags =
      filestr.setf(ios::scientific, ios::floatfield);
   int old_precision = filestr.precision(8);

   write_CSE6SSM_inputs(coefficients_inputs, filestr, width);

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << coefficients_m0Sqr << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << low_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << susy_scale << ' ';

   if (!filestr.good()) {
      ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_line: "
            "file stream is corrupted");
      return;
   }
   filestr << std::left << std::setw(width) << high_scale << ' ';

   for (std::size_t p = 0; p < coefficients.size(); ++p) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_semianalytic_coefficients_writer::write_coefficients_line: "
               "file stream is corrupted");
         break;
      }

      const std::valarray<double>& coeffs = coefficients[p].values;
      const std::size_t num_coeffs = coeffs.size();

      for (std::size_t i = 0; i < num_coeffs; ++i) {
         filestr << std::left << std::setw(width) << coeffs[i] << ' ';
      }
   }

   filestr << std::left << std::setw(width) << coefficients_problems.have_problem() << ' ';

   if (coefficients_problems.have_problem() || coefficients_problems.have_warning()) {
      filestr << "\t# " << coefficients_problems << '\n';
   } else {
      filestr << '\n';
   }

   filestr.setf(old_flags);
   filestr.precision(old_precision);
}

} // namespace flexiblesusy

