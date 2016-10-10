// ====================================================================
// Implementation of convergence tester for semianalytic version of
// the two-scale algorithm, checking the SUSY parameters for
// convergence
// ====================================================================

#ifndef CSE6SSM_SUSY_TWO_SCALE_CONVERGENCE_TESTER_H
#define CSE6SSM_SUSY_TWO_SCALE_CONVERGENCE_TESTER_H

#include "CSE6SSM_susy_convergence_tester.hpp"
#include "CSE6SSM_semi_two_scale_model.hpp"
#include "two_scale_convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class CSE6SSM_susy_convergence_tester<Two_scale> : public Convergence_tester_DRbar<CSE6SSM_semianalytic<Two_scale> > {
public:
   CSE6SSM_susy_convergence_tester(CSE6SSM_semianalytic<Two_scale>*, double);
   virtual ~CSE6SSM_susy_convergence_tester();

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
