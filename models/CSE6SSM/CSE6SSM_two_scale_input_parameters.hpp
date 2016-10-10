// ====================================================================
// Specialisation of input parameters for models solved using the
// two scale algorithm
// ====================================================================

#ifndef CSE6SSM_TWO_SCALE_INPUT_PARAMETERS_H
#define CSE6SSM_TWO_SCALE_INPUT_PARAMETERS_H

#include "CSE6SSM_input_parameters.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

template<>
struct CSE6SSM_input_parameters<Two_scale> {
   double m0;
   double m12;
   double TanBeta;
   int SignLambdax;
   double Azero;
   double sInput;
   double QSInput;
   Eigen::Matrix<double,3,2> hEInput;
   double SigmaLInput;
   double KappaPrInput;
   double SigmaxInput;
   Eigen::Matrix<double,3,3> gDInput;
   Eigen::Matrix<double,3,3> KappaInput;
   Eigen::Matrix<double,2,2> Lambda12Input;
   Eigen::Matrix<double,3,2> fuInput;
   Eigen::Matrix<double,3,2> fdInput;
   double MuPrInput;
   double MuPhiInput;
   double BMuPrInput;
   double BMuPhiInput;

   CSE6SSM_input_parameters()
      : m0(0), m12(0), TanBeta(0), SignLambdax(1), Azero(0), sInput(0), QSInput(0),
   hEInput(Eigen::Matrix<double,3,2>::Zero()), SigmaLInput(0), KappaPrInput(0),
   SigmaxInput(0), gDInput(Eigen::Matrix<double,3,3>::Zero()), KappaInput(
   Eigen::Matrix<double,3,3>::Zero()), Lambda12Input(Eigen::Matrix<double,2,2>
   ::Zero()), fuInput(Eigen::Matrix<double,3,2>::Zero()), fdInput(Eigen::Matrix
   <double,3,2>::Zero()), MuPrInput(0), MuPhiInput(0), BMuPrInput(0),
   BMuPhiInput(0)

   {}
};

std::ostream& operator<<(std::ostream&, const CSE6SSM_input_parameters<Two_scale>&);

} // namespace flexiblesusy

#endif
