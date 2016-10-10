// ====================================================================
// Computes an approximate 1-loop upper bound on the lightest CP-even
// Higgs
// ====================================================================

#ifndef CSE6SSM_HIGGS_UPPER_BOUND_H
#define CSE6SSM_HIGGS_UPPER_BOUND_H

#include "CSE6SSM_mass_eigenstates.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

class CSE6SSM_higgs_upper_bound {
public:
   CSE6SSM_higgs_upper_bound(const CSE6SSM_mass_eigenstates&);
   ~CSE6SSM_higgs_upper_bound();

   // if false (default), include only the third
   // generation contributions
   void set_include_all_SM_generations(bool flag) { include_all_gens = flag; }
   void set_include_up_tadpoles(bool flag) { include_ups = flag; }
   void set_include_down_tadpoles(bool flag) { include_downs = flag; }
   void set_include_exotic_tadpoles(bool flag) { include_exotics = flag; }
   void set_include_inert_singlet_tadpoles(bool flag) { include_inert_singlets = flag; }
   void set_include_inert_neutral_higgs_tadpoles(bool flag) {
      include_inert_neutral_higgs = flag;
   }
   void set_include_inert_charged_higgs_tadpoles(bool flag) {
      include_inert_charged_higgs = flag;
   }

   Eigen::Array<double,2,1> calculate_MSu2(unsigned) const;
   double calculate_Sin2ThetaSu(unsigned) const;
   double calculate_Cos2ThetaSu(unsigned) const;
   double calculate_MFu2(unsigned) const;

   Eigen::Array<double,2,1> calculate_MSd2(unsigned) const;
   double calculate_Sin2ThetaSd(unsigned) const;
   double calculate_Cos2ThetaSd(unsigned) const;
   double calculate_MFd2(unsigned) const;

   Eigen::Array<double,2,1> calculate_MSDX2(unsigned) const;
   double calculate_Sin2ThetaSDX(unsigned) const;
   double calculate_Cos2ThetaSDX(unsigned) const;
   double calculate_MFDX2(unsigned) const;

   Eigen::Array<double,2,1> calculate_MHI02(unsigned) const;
   double calculate_Sin2ThetaHI0(unsigned) const;
   double calculate_Cos2ThetaHI0(unsigned) const;
   double calculate_MFHI02(unsigned) const;

   double calculate_MSI02(unsigned) const;

   Eigen::Array<double,2,1> calculate_MHIPM2(unsigned) const;
   double calculate_Sin2ThetaHIPM(unsigned) const;
   double calculate_Cos2ThetaHIPM(unsigned) const;
   double calculate_MFHIPM2(unsigned) const;

   double get_upper_bound() const { return upper_bound; }

   double calculate_tree_level_upper_bound();
   double calculate_one_loop_upper_bound();

   double get_tadpole_vd() const;
   double get_tadpole_vu() const;
   double get_up_contribution(unsigned) const;
   double get_down_contribution(unsigned) const;
   double get_exotic_contribution(unsigned) const;
   double get_inert_singlet_contribution(unsigned) const;
   double get_inert_neutral_higgs_contribution(unsigned) const;
   double get_inert_charged_higgs_contribution(unsigned) const;

   double get_unrotated_up_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_down_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_exotic_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_inert_singlet_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_inert_neutral_higgs_contribution(unsigned,unsigned,unsigned) const;
   double get_unrotated_inert_charged_higgs_contribution(unsigned,unsigned,unsigned) const;

   // methods for comparing against 1-loop self-energies
   std::complex<double> get_full_self_energy(double, unsigned, unsigned) const;
   std::complex<double> get_unrotated_full_self_energy(double, unsigned, unsigned) const;

   std::complex<double> get_self_energy(double, unsigned, unsigned) const;
   std::complex<double> get_unrotated_self_energy(double, unsigned, unsigned) const;

   std::complex<double> get_up_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_down_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_exotic_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_inert_higgs_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_inert_charged_higgs_self_energy(double,unsigned,unsigned) const;

   std::complex<double> get_unrotated_up_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_unrotated_down_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_unrotated_exotic_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_unrotated_inert_higgs_self_energy(double,unsigned,unsigned) const;
   std::complex<double> get_unrotated_inert_charged_higgs_self_energy(double,unsigned,unsigned) const;

private:
   CSE6SSM_mass_eigenstates model;
   double upper_bound;
   bool include_all_gens;
   bool include_ups;
   bool include_downs;
   bool include_exotics;
   bool include_inert_singlets;
   bool include_inert_neutral_higgs;
   bool include_inert_charged_higgs;

   Eigen::Matrix<double,2,2> get_mass_matrix_Su(unsigned) const;
   Eigen::Matrix<double,2,2> get_mass_matrix_Sd(unsigned) const;
   Eigen::Matrix<double,2,2> get_mass_matrix_SDX(unsigned) const;
   Eigen::Matrix<double,2,2> get_mass_matrix_HI0(unsigned) const;
   Eigen::Matrix<double,2,2> get_mass_matrix_HIPM(unsigned) const;

   // wrappers for Passarino-Veltman functions
   double A0(double) const;
   double B0(double,double,double) const;
   double G0(double,double,double) const;

   // convenient helper functions for derivatives
   Eigen::Matrix<double,2,2> get_dmass_matrix_Su_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_Su_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_Su_dvd_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_Su_dvd_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_Su_dvu_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_Sd_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_Sd_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_Sd_dvd_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_Sd_dvd_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_Sd_dvu_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_SDX_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_SDX_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_SDX_dvd_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_SDX_dvd_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_SDX_dvu_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_HI0_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_HI0_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_HI0_dvd_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_HI0_dvd_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_HI0_dvu_dvu(unsigned) const;

   Eigen::Matrix<double,2,2> get_dmass_matrix_HIPM_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_dmass_matrix_HIPM_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_HIPM_dvd_dvd(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_HIPM_dvd_dvu(unsigned) const;
   Eigen::Matrix<double,2,2> get_d2mass_matrix_HIPM_dvu_dvu(unsigned) const;

   double get_dV1lp_up_dvd(unsigned) const;
   double get_dV1lp_up_dvu(unsigned) const;
   double get_d2V1lp_up_dvd_dvd(unsigned) const;
   double get_d2V1lp_up_dvd_dvu(unsigned) const;
   double get_d2V1lp_up_dvu_dvu(unsigned) const;

   double get_dV1lp_down_dvd(unsigned) const;
   double get_dV1lp_down_dvu(unsigned) const;
   double get_d2V1lp_down_dvd_dvd(unsigned) const;
   double get_d2V1lp_down_dvd_dvu(unsigned) const;
   double get_d2V1lp_down_dvu_dvu(unsigned) const;

   double get_dV1lp_exotic_dvd(unsigned) const;
   double get_dV1lp_exotic_dvu(unsigned) const;
   double get_d2V1lp_exotic_dvd_dvd(unsigned) const;
   double get_d2V1lp_exotic_dvd_dvu(unsigned) const;
   double get_d2V1lp_exotic_dvu_dvu(unsigned) const;

   double get_dV1lp_inert_singlet_dvd(unsigned) const;
   double get_dV1lp_inert_singlet_dvu(unsigned) const;
   double get_d2V1lp_inert_singlet_dvd_dvd(unsigned) const;
   double get_d2V1lp_inert_singlet_dvd_dvu(unsigned) const;
   double get_d2V1lp_inert_singlet_dvu_dvu(unsigned) const;

   double get_dV1lp_inert_neutral_higgs_dvd(unsigned) const;
   double get_dV1lp_inert_neutral_higgs_dvu(unsigned) const;
   double get_d2V1lp_inert_neutral_higgs_dvd_dvd(unsigned) const;
   double get_d2V1lp_inert_neutral_higgs_dvd_dvu(unsigned) const;
   double get_d2V1lp_inert_neutral_higgs_dvu_dvu(unsigned) const;

   double get_dV1lp_inert_charged_higgs_dvd(unsigned) const;
   double get_dV1lp_inert_charged_higgs_dvu(unsigned) const;
   double get_d2V1lp_inert_charged_higgs_dvd_dvd(unsigned) const;
   double get_d2V1lp_inert_charged_higgs_dvd_dvu(unsigned) const;
   double get_d2V1lp_inert_charged_higgs_dvu_dvu(unsigned) const;
};

} // namespace flexiblesusy

#endif
