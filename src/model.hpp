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

#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <ostream>
#include <typeinfo>

namespace flexiblesusy {

class Model {
public:
   virtual ~Model() {};
   virtual void calculate_spectrum() = 0;
   virtual std::string name() const = 0;
   virtual void run_to(double, double eps = -1.0) = 0;
   virtual void print(std::ostream& out) const = 0;

   // general parameter getter and setter
   virtual double get_parameter(unsigned) const = 0;
   virtual void set_parameter(unsigned, double) = 0;
};

template <class TargetModel, class AbstractModel>
TargetModel cast_model(AbstractModel abstract_model)
{
#ifdef ENABLE_DEBUG
   TargetModel tmp = dynamic_cast<TargetModel>(abstract_model);
   if (!tmp) {
      FATAL("model " << abstract_model << " is not of type "
            << typeid(TargetModel).name());
   }
   return tmp;
#else
   return static_cast<TargetModel>(abstract_model);
#endif
}

} // namespace flexiblesusy

#endif
