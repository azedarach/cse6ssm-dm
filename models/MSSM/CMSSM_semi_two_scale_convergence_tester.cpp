
#include "CMSSM_semi_two_scale_convergence_tester.hpp"
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

CMSSM_semianalytic_convergence_tester<Two_scale>::CMSSM_semianalytic_convergence_tester(CMSSM_semianalytic<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CMSSM_semianalytic<Two_scale> >(model, accuracy_goal)
{
}

CMSSM_semianalytic_convergence_tester<Two_scale>::~CMSSM_semianalytic_convergence_tester()
{
}

double CMSSM_semianalytic_convergence_tester<Two_scale>::max_rel_diff() const
{
   const CMSSM_semianalytic<Two_scale>& ol = get_last_iteration_model();
   const CMSSM_semianalytic<Two_scale>& ne = get_model();

   double diff[34] = { 0 };

   diff[0] = MaxRelDiff(OLD(MGlu),NEW(MGlu));
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 1] = MaxRelDiff(OLD1(MSd,i),NEW1(MSd,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 7] = MaxRelDiff(OLD1(MSv,i),NEW1(MSv,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 10] = MaxRelDiff(OLD1(MSu,i),NEW1(MSu,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 16] = MaxRelDiff(OLD1(MSe,i),NEW1(MSe,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 22] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (unsigned i = 1; i < 2; i++) {
      diff[i + 24] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (unsigned i = 1; i < 2; i++) {
      diff[i + 26] = MaxRelDiff(OLD1(MHpm,i),NEW1(MHpm,i));
   }
   for (unsigned i = 0; i < 4; i++) {
      diff[i + 28] = MaxRelDiff(OLD1(MChi,i),NEW1(MChi,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 32] = MaxRelDiff(OLD1(MCha,i),NEW1(MCha,i));
   }

   return *std::max_element(diff, diff + 34);

}

} // namespace flexiblesusy
