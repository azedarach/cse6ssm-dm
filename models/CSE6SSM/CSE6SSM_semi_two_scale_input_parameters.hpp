// ====================================================================
// Specialisation of input parameters for models solved using the
// semianalytic version of the two-scale algorithm
// ====================================================================

#ifndef CSE6SSM_SEMI_TWO_SCALE_INPUT_PARAMETERS_H
#define CSE6SSM_SEMI_TWO_SCALE_INPUT_PARAMETERS_H

#include "CSE6SSM_semi_input_parameters.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template<>
struct CSE6SSM_semianalytic_input_parameters<Two_scale> {
   double m12;
   double Azero;
   double TanBeta;
   double sInput;
   double QSInput;
   Eigen::Matrix<double,3,2> hEInput;
   double SigmaLInput;
   double KappaPrInput;
   double SigmaxInput;
   Eigen::Matrix<double,3,3> gDInput;
   Eigen::Matrix<double,3,3> KappaInput;
   Eigen::Matrix<double,2,2> Lambda12Input;
   double LambdaxInput;
   Eigen::Matrix<double,3,2> fuInput;
   Eigen::Matrix<double,3,2> fdInput;
   double MuPrInput;
   double MuPhiInput;
   double BMuPrInput;
   double BMuPhiInput;

   CSE6SSM_semianalytic_input_parameters()
      : m12(0), Azero(0), TanBeta(0), sInput(0), QSInput(0)
      , hEInput(Eigen::Matrix<double,3,2>::Zero()), SigmaLInput(0)
      , KappaPrInput(0), SigmaxInput(0)
      , gDInput(Eigen::Matrix<double,3,3>::Zero())
      , KappaInput(Eigen::Matrix<double,3,3>::Zero())
      , Lambda12Input(Eigen::Matrix<double,2,2>::Zero()), LambdaxInput(0)
      , fuInput(Eigen::Matrix<double,3,2>::Zero())
      , fdInput(Eigen::Matrix<double,3,2>::Zero()), MuPrInput(0)
      , MuPhiInput(0), BMuPrInput(0), BMuPhiInput(0)
      {}
};

std::ostream& operator<<(std::ostream&, const CSE6SSM_semianalytic_input_parameters<Two_scale>&);

} // namespace flexiblesusy

#endif
