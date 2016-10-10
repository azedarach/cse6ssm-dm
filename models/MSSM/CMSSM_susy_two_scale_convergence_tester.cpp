
#include "CMSSM_susy_two_scale_convergence_tester.hpp"
#include "wrappers.hpp"

#include <cmath>
#include <algorithm>

namespace flexiblesusy {

#define OLD(p) ol.get_##p()
#define NEW(p) ne.get_##p()

#define OLD1(p,i) ol.get_##p()(i)
#define NEW1(p,i) ne.get_##p()(i)

#define OLD2(p,i,j) ol.get_##p(i,j)
#define NEW2(p,i,j) ne.get_##p(i,j)

#define OLD3(p,i,j,k) ol.get_##p(i,j,k)
#define NEW3(p,i,j,k) ne.get_##p(i,j,k)

#define OLD4(p,i,j,k,l) ol.get_##p(i,j,k,l)
#define NEW4(p,i,j,k,l) ne.get_##p(i,j,k,l)

CMSSM_susy_convergence_tester<Two_scale>::CMSSM_susy_convergence_tester(CMSSM_semianalytic<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CMSSM_semianalytic<Two_scale> >(model, accuracy_goal)
{
}

CMSSM_susy_convergence_tester<Two_scale>::~CMSSM_susy_convergence_tester()
{
}

double CMSSM_susy_convergence_tester<Two_scale>::max_rel_diff() const
{
   const CMSSM_semianalytic<Two_scale>& ol = get_last_iteration_model();
   const CMSSM_semianalytic<Two_scale>& ne = get_model();

   double diff[33] = { 0. };

   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i] = MaxRelDiff(OLD2(Yd,i,j),NEW2(Yd,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 9] = MaxRelDiff(OLD2(Ye,i,j),NEW2(Ye,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 18] = MaxRelDiff(OLD2(Yu,i,j),NEW2(Yu,i,j));
      }
   }
   diff[27] = MaxRelDiff(OLD(Mu),NEW(Mu));
   diff[28] = MaxRelDiff(OLD(g1),NEW(g1));
   diff[29] = MaxRelDiff(OLD(g2),NEW(g2));
   diff[30] = MaxRelDiff(OLD(g3),NEW(g3));
   diff[31] = MaxRelDiff(OLD(vd),NEW(vd));
   diff[32] = MaxRelDiff(OLD(vu),NEW(vu));

   return *std::max_element(diff, diff + 33);
}

} // namespace flexiblesusy
