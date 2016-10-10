// ====================================================================
// Implementation of convergence tester for semianalytic version of
// the two-scale algorithm, checking the SUSY parameters for
// convergence
// ====================================================================

#include "CSE6SSM_susy_two_scale_convergence_tester.hpp"
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

CSE6SSM_susy_convergence_tester<Two_scale>::CSE6SSM_susy_convergence_tester(CSE6SSM_semianalytic<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CSE6SSM_semianalytic<Two_scale> >(model, accuracy_goal)
{
}

CSE6SSM_susy_convergence_tester<Two_scale>::~CSE6SSM_susy_convergence_tester()
{
}

double CSE6SSM_susy_convergence_tester<Two_scale>::max_rel_diff() const
{
   const CSE6SSM_semianalytic<Two_scale>& ol = get_last_iteration_model();
   const CSE6SSM_semianalytic<Two_scale>& ne = get_model();

   double diff[84] = { 0. };

   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i] = MaxRelDiff(OLD2(Yd,i,j),NEW2(Yd,i,j));
      }
   }
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 9] = MaxRelDiff(OLD2(hE,i,j),NEW2(hE,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 15] = MaxRelDiff(OLD2(Ye,i,j),NEW2(Ye,i,j));
      }
   }
   diff[24] = MaxRelDiff(OLD(SigmaL),NEW(SigmaL));
   diff[25] = MaxRelDiff(OLD(KappaPr),NEW(KappaPr));
   diff[26] = MaxRelDiff(OLD(Sigmax),NEW(Sigmax));
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 27] = MaxRelDiff(OLD2(gD,i,j),NEW2(gD,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 36] = MaxRelDiff(OLD2(Kappa,i,j),NEW2(Kappa,i,j));
      }
   }
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 2; ++i) {
         diff[2 * j + i + 45] = MaxRelDiff(OLD2(Lambda12,i,j),NEW2(Lambda12,i,j));
      }
   }
   diff[49] = MaxRelDiff(OLD(Lambdax),NEW(Lambdax));
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 50] = MaxRelDiff(OLD2(fu,i,j),NEW2(fu,i,j));
      }
   }
   for (unsigned j = 0; j < 2; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 56] = MaxRelDiff(OLD2(fd,i,j),NEW2(fd,i,j));
      }
   }
   for (unsigned j = 0; j < 3; ++j) {
      for (unsigned i = 0; i < 3; ++i) {
         diff[3 * j + i + 62] = MaxRelDiff(OLD2(Yu,i,j),NEW2(Yu,i,j));
      }
   }
   diff[71] = MaxRelDiff(OLD(MuPr),NEW(MuPr));
   diff[72] = MaxRelDiff(OLD(MuPhi),NEW(MuPhi));
   diff[73] = MaxRelDiff(OLD(XiF),NEW(XiF));
   diff[74] = MaxRelDiff(OLD(g1),NEW(g1));
   diff[75] = MaxRelDiff(OLD(g2),NEW(g2));
   diff[76] = MaxRelDiff(OLD(g3),NEW(g3));
   diff[77] = MaxRelDiff(OLD(g1p),NEW(g1p));
   diff[78] = MaxRelDiff(OLD(vd),NEW(vd));
   diff[79] = MaxRelDiff(OLD(vu),NEW(vu));
   diff[80] = MaxRelDiff(OLD(vs),NEW(vs));
   diff[81] = MaxRelDiff(OLD(vsb),NEW(vsb));
   diff[82] = MaxRelDiff(OLD(vphi),NEW(vphi));
   diff[83] = MaxRelDiff(OLD(QS),NEW(QS));

   return *std::max_element(diff, diff + 84);

}

} // namespace flexiblesusy
