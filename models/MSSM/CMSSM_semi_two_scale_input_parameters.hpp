
#ifndef CMSSM_SEMI_TWO_SCALE_INPUT_PARAMETERS_H
#define CMSSM_SEMI_TWO_SCALE_INPUT_PARAMETERS_H

#include "CMSSM_semi_input_parameters.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template<>
struct CMSSM_semianalytic_input_parameters<Two_scale> {
   double m12;
   double Azero;
   double TanBeta;
   double MuInput;
   bool MuInput_at_MS;

   CMSSM_semianalytic_input_parameters()
      : m12(0), Azero(0), TanBeta(0), MuInput(0), MuInput_at_MS(false)
      {}
};

std::ostream& operator<<(std::ostream&, const CMSSM_semianalytic_input_parameters<Two_scale>&);

} // namespace flexiblesusy

#endif
