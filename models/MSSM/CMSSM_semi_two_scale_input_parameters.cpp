
#include "CMSSM_semi_two_scale_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const CMSSM_semianalytic_input_parameters<Two_scale>& input)
{
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "MuInput_at_MS = " << INPUT(MuInput_at_MS) << ", ";

   return ostr;
}

} // namespace flexiblesusy
