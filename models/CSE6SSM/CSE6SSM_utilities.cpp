// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Wed 3 Jun 2015 23:47:48

#include "CSE6SSM_utilities.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#define PHYSICAL(p) model.get_physical().p
#define MODELPARAMETER(p) model.get_##p()

namespace flexiblesusy {

CSE6SSM_spectrum_plotter::CSE6SSM_spectrum_plotter()
   : spectrum()
   , scale(0.0)
   , width(16)
{
}

void CSE6SSM_spectrum_plotter::write_to_file(const std::string& file_name) const
{
   if (spectrum.empty())
      return;

   std::ofstream filestr(file_name.c_str(), std::ios::out);
   VERBOSE_MSG("CSE6SSM_spectrum_plotter::write_to_file: opening file: "
               << file_name.c_str());
   if (filestr.fail()) {
      ERROR("CSE6SSM_spectrum_plotter::write_to_file: can't open file "
            << file_name);
      return;
   }

   filestr << "### one-loop pole masses (Q = " << scale << " GeV)\n";
   write_spectrum(spectrum, filestr);

   filestr.close();
   VERBOSE_MSG("CSE6SSM_spectrum_plotter::write_to_file: file written: "
               << file_name.c_str());
}

void CSE6SSM_spectrum_plotter::extract_spectrum(const CSE6SSM_mass_eigenstates& model)
{
   spectrum.clear();
   scale = model.get_scale();

   spectrum.push_back(TParticle("Glu", "\\tilde{g}", to_valarray(PHYSICAL(MGlu))));
   spectrum.push_back(TParticle("ChaP", "\\tilde{\\chi}^{'-}", to_valarray(PHYSICAL(MChaP))));
   spectrum.push_back(TParticle("VZp", "{Z'}", to_valarray(PHYSICAL(MVZp))));
   spectrum.push_back(TParticle("Sd", "\\tilde{d}", to_valarray(PHYSICAL(MSd))));
   spectrum.push_back(TParticle("Sv", "\\tilde{\\nu}", to_valarray(PHYSICAL(MSv))));
   spectrum.push_back(TParticle("Su", "\\tilde{u}", to_valarray(PHYSICAL(MSu))));
   spectrum.push_back(TParticle("Se", "\\tilde{e}", to_valarray(PHYSICAL(MSe))));
   spectrum.push_back(TParticle("SDX", "\\tilde{x}", to_valarray(PHYSICAL(MSDX))));
   spectrum.push_back(TParticle("hh", "h", to_valarray(PHYSICAL(Mhh))));
   spectrum.push_back(TParticle("Ah", "A^0", to_valarray(PHYSICAL(MAh))));
   spectrum.push_back(TParticle("Hpm", "H^-", to_valarray(PHYSICAL(MHpm))));
   spectrum.push_back(TParticle("Chi", "\\tilde{\\chi}^0", to_valarray(PHYSICAL(MChi))));
   spectrum.push_back(TParticle("Cha", "\\tilde{\\chi}^-", to_valarray(PHYSICAL(MCha))));
   spectrum.push_back(TParticle("FDX", "x", to_valarray(PHYSICAL(MFDX))));
   spectrum.push_back(TParticle("SHI0", "h^{0,Inert}", to_valarray(PHYSICAL(MSHI0))));
   spectrum.push_back(TParticle("SHIPM", "h^{-,Inert}", to_valarray(PHYSICAL(MSHIPM))));
   spectrum.push_back(TParticle("ChaI", "\\tilde{\\chi}^{-,Inert}", to_valarray(PHYSICAL(MChaI))));
   spectrum.push_back(TParticle("ChiI", "\\tilde{\\chi}^{0,Inert}", to_valarray(PHYSICAL(MChiI))));
   spectrum.push_back(TParticle("SHp0", "H^{'0}", to_valarray(PHYSICAL(MSHp0))));
   spectrum.push_back(TParticle("SHpp", "H^{'-}", to_valarray(PHYSICAL(MSHpp))));
   spectrum.push_back(TParticle("ChiP", "\\tilde{\\chi}^{'0}", to_valarray(PHYSICAL(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      spectrum.push_back(TParticle("Fd", "d", to_valarray(PHYSICAL(MFd))));
      spectrum.push_back(TParticle("Fe", "e", to_valarray(PHYSICAL(MFe))));
      spectrum.push_back(TParticle("Fu", "u", to_valarray(PHYSICAL(MFu))));
      spectrum.push_back(TParticle("Fv", "\\nu", to_valarray(PHYSICAL(MFv))));
      spectrum.push_back(TParticle("VG", "g", to_valarray(PHYSICAL(MVG))));
      spectrum.push_back(TParticle("VP", "\\gamma", to_valarray(PHYSICAL(MVP))));
      spectrum.push_back(TParticle("VWm", "W^-", to_valarray(PHYSICAL(MVWm))));
      spectrum.push_back(TParticle("VZ", "Z", to_valarray(PHYSICAL(MVZ))));

   }
}

void CSE6SSM_spectrum_plotter::write_spectrum(const TSpectrum& spectrum, std::ofstream& filestr) const
{
   for (std::size_t s = 0; s < spectrum.size(); ++s) {
      if (!filestr.good()) {
         ERROR("CSE6SSM_spectrum_plotter::write_spectrum: "
               "file stream is corrupted");
         break;
      }

      const std::string& name = spectrum[s].name;
      const std::string& latex_name = spectrum[s].latex_name;
      const std::valarray<double>& masses = spectrum[s].masses;
      const std::size_t multiplicity = masses.size();

      filestr << std::left << "# " << name << '\n';

      for (std::size_t i = 0; i < multiplicity; ++i) {
         std::string lname("$" + latex_name + "$");
         std::stringstream lname_with_index;
         lname_with_index << "$" << latex_name;
         if (multiplicity > 1)
            lname_with_index << "_{" << (i+1) << "}";
         lname_with_index  << "$";

         filestr << std::left << std::setw(width) << s
                 << std::left << std::setw(width) << masses[i]
                 << std::left << std::setw(width) << name
                 << std::left << std::setw(2*width) << lname
                 << std::left << std::setw(2*width) << lname_with_index.str()
                 << '\n';
      }
   }
}

std::valarray<double> CSE6SSM_spectrum_plotter::to_valarray(double v)
{
   return std::valarray<double>(&v, 1);
}

} // namespace flexiblesusy
