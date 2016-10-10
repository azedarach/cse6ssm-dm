// ====================================================================
// Test implementation of low scale constraint to be used in
// semianalytic version of the two-scale algorithm
// ====================================================================

#ifndef CSE6SSM_SEMI_TWO_SCALE_LOW_SCALE_CONSTRAINT_H
#define CSE6SSM_SEMI_TWO_SCALE_LOW_SCALE_CONSTRAINT_H

#include "CSE6SSM_semi_low_scale_constraint.hpp"
#include "two_scale_constraint.hpp"
#include "lowe.h"
#include <Eigen/Core>

namespace flexiblesusy {

template <class T>
class CSE6SSM_semianalytic;

class Two_scale;

template<>
class CSE6SSM_semianalytic_low_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   CSE6SSM_semianalytic_low_scale_constraint();
   CSE6SSM_semianalytic_low_scale_constraint(CSE6SSM_semianalytic<Two_scale>*, const QedQcd&);
   virtual ~CSE6SSM_semianalytic_low_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm() const;
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns() const;
   double get_initial_scale_guess() const;
   void initialize();
   const QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const QedQcd&);
   void set_threshold_corrections_loop_order(unsigned); ///< threshold corrections loop order

private:
   double scale;
   double initial_scale_guess;
   CSE6SSM_semianalytic<Two_scale>* model;
   QedQcd oneset;
   Eigen::Matrix<std::complex<double>,3,3> ckm;
   Eigen::Matrix<std::complex<double>,3,3> pmns;
   Eigen::Matrix<double,3,3> neutrinoDRbar;
   double MWDRbar;
   double MZDRbar;
   double AlphaS;
   double EDRbar;
   double ThetaWDRbar;
   double new_g1, new_g2, new_g3;
   double self_energy_w_at_mw;
   unsigned threshold_corrections_loop_order; ///< threshold corrections loop order

   double calculate_theta_w(double);
   void calculate_threshold_corrections();
   void calculate_DRbar_gauge_couplings();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   void calculate_MNeutrino_DRbar();
   double calculate_delta_alpha_em(double) const;
   double calculate_delta_alpha_s(double) const;
   void recalculate_mw_pole();
   void update_scale();
};

} // namespace flexiblesusy

#endif
