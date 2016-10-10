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

// File generated at Fri 10 Jul 2015 12:04:44

#ifndef CMSSM_TWO_SCALE_INPUT_PARAMETERS_H
#define CMSSM_TWO_SCALE_INPUT_PARAMETERS_H

#include "CMSSM_input_parameters.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template<>
struct CMSSM_input_parameters<Two_scale> {
   double m0;
   double m12;
   double TanBeta;
   int SignMu;
   double Azero;

   CMSSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), SignMu(1), Azero(0)

   {}
};

std::ostream& operator<<(std::ostream&, const CMSSM_input_parameters<Two_scale>&);

} // namespace flexiblesusy

#endif
