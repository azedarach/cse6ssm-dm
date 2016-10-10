// ====================================================================
// Specialisation of input parameters for models solved using the
// semianalytic version of the two-scale algorithm
// ====================================================================

#include "CSE6SSM_semi_two_scale_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const CSE6SSM_semianalytic_input_parameters<Two_scale>& input)
{
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "sInput = " << INPUT(sInput) << ", ";
   ostr << "QSInput = " << INPUT(QSInput) << ", ";
   ostr << "hEInput = " << INPUT(hEInput) << ", ";
   ostr << "SigmaLInput = " << INPUT(SigmaLInput) << ", ";
   ostr << "KappaPrInput = " << INPUT(KappaPrInput) << ", ";
   ostr << "SigmaxInput = " << INPUT(SigmaxInput) << ", ";
   ostr << "gDInput = " << INPUT(gDInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "Lambda12Input = " << INPUT(Lambda12Input) << ", ";
   ostr << "LambdaxInput = " << INPUT(LambdaxInput) << ", ";
   ostr << "fuInput = " << INPUT(fuInput) << ", ";
   ostr << "fdInput = " << INPUT(fdInput) << ", ";
   ostr << "MuPrInput = " << INPUT(MuPrInput) << ", ";
   ostr << "MuPhiInput = " << INPUT(MuPhiInput) << ", ";
   ostr << "BMuPrInput = " << INPUT(BMuPrInput) << ", ";
   ostr << "BMuPhiInput = " << INPUT(BMuPhiInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
