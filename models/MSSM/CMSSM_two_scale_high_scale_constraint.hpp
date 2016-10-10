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

// File generated at Fri 10 Jul 2015 12:04:45

#ifndef CMSSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define CMSSM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "CMSSM_high_scale_constraint.hpp"
#include "CMSSM_two_scale_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class CMSSM;

class Two_scale;

template<>
class CMSSM_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CMSSM_high_scale_constraint();
   CMSSM_high_scale_constraint(CMSSM<Two_scale>*);
   virtual ~CMSSM_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const CMSSM_input_parameters<Two_scale>& get_input_parameters() const;
   CMSSM<Two_scale>* get_model() const;
   void initialize();
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();

private:
   double scale;
   double initial_scale_guess;
   CMSSM<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
