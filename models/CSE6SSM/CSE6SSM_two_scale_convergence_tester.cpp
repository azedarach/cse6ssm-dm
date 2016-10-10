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

// File generated at Wed 3 Jun 2015 23:47:46

#include "CSE6SSM_two_scale_convergence_tester.hpp"
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

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

CSE6SSM_convergence_tester<Two_scale>::CSE6SSM_convergence_tester(CSE6SSM<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<CSE6SSM<Two_scale> >(model, accuracy_goal)
{
}

CSE6SSM_convergence_tester<Two_scale>::~CSE6SSM_convergence_tester()
{
}

double CSE6SSM_convergence_tester<Two_scale>::max_rel_diff() const
{
   const CSE6SSM<Two_scale>& ol = get_last_iteration_model();
   const CSE6SSM<Two_scale>& ne = get_model();

   double diff[81] = { 0 };

   diff[0] = MaxRelDiff(OLD(MGlu),NEW(MGlu));
   diff[1] = MaxRelDiff(OLD(MChaP),NEW(MChaP));
   diff[2] = MaxRelDiff(OLD(MVZp),NEW(MVZp));
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 3] = MaxRelDiff(OLD1(MSd,i),NEW1(MSd,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 9] = MaxRelDiff(OLD1(MSv,i),NEW1(MSv,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 12] = MaxRelDiff(OLD1(MSu,i),NEW1(MSu,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 18] = MaxRelDiff(OLD1(MSe,i),NEW1(MSe,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 24] = MaxRelDiff(OLD1(MSDX,i),NEW1(MSDX,i));
   }
   for (unsigned i = 0; i < 5; i++) {
      diff[i + 30] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (unsigned i = 2; i < 5; i++) {
      diff[i + 35] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (unsigned i = 1; i < 2; i++) {
      diff[i + 40] = MaxRelDiff(OLD1(MHpm,i),NEW1(MHpm,i));
   }
   for (unsigned i = 0; i < 8; i++) {
      diff[i + 42] = MaxRelDiff(OLD1(MChi,i),NEW1(MChi,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 50] = MaxRelDiff(OLD1(MCha,i),NEW1(MCha,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 52] = MaxRelDiff(OLD1(MFDX,i),NEW1(MFDX,i));
   }
   for (unsigned i = 0; i < 7; i++) {
      diff[i + 55] = MaxRelDiff(OLD1(MSHI0,i),NEW1(MSHI0,i));
   }
   for (unsigned i = 0; i < 4; i++) {
      diff[i + 62] = MaxRelDiff(OLD1(MSHIPM,i),NEW1(MSHIPM,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 66] = MaxRelDiff(OLD1(MChaI,i),NEW1(MChaI,i));
   }
   for (unsigned i = 0; i < 7; i++) {
      diff[i + 68] = MaxRelDiff(OLD1(MChiI,i),NEW1(MChiI,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 75] = MaxRelDiff(OLD1(MSHp0,i),NEW1(MSHp0,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 77] = MaxRelDiff(OLD1(MSHpp,i),NEW1(MSHpp,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 79] = MaxRelDiff(OLD1(MChiP,i),NEW1(MChiP,i));
   }

   return *std::max_element(diff, diff + 81);

}

} // namespace flexiblesusy
