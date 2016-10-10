// ====================================================================
// Test implementation of wrapper class for model class using
// semianalytic version of the two-scale algorithm
// ====================================================================

/**
 * @file CSE6SSM_semi_two_scale_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

#ifndef CSE6SSM_SEMI_TWO_SCALE_SLHA_H
#define CSE6SSM_SEMI_TWO_SCALE_SLHA_H

#include "CSE6SSM_semi_two_scale_model.hpp"
#include "CSE6SSM_physical.hpp"
#include "CSE6SSM_semi_model_slha.hpp"

namespace flexiblesusy {

class Two_scale;

/**
 * @class CSE6SSM_semianalytic_slha<Two_scale>
 * @brief model class wrapper in SLHA convention
 */

template<>
class CSE6SSM_semianalytic_slha<Two_scale> : public CSE6SSM_semianalytic<Two_scale> {
public:
   explicit CSE6SSM_semianalytic_slha(const CSE6SSM_semianalytic_input_parameters<Two_scale>& input_ = CSE6SSM_semianalytic_input_parameters<Two_scale>());
   explicit CSE6SSM_semianalytic_slha(const CSE6SSM_semianalytic<Two_scale>&);
   virtual ~CSE6SSM_semianalytic_slha();

   virtual void clear();
   void convert_to_slha(); ///< converts pole masses and couplings to SLHA convention
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm_matrix() const { return ckm; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns_matrix() const { return pmns; }
   const CSE6SSM_physical& get_physical_slha() const; ///< returns pole masses to SLHA convention
   CSE6SSM_physical& get_physical_slha(); ///< return pole masses to SLHA convention
   const CSE6SSM_physical& get_drbar_slha() const; ///< returns DR-bar masses to SLHA convention
   CSE6SSM_physical& get_drbar_slha(); ///< return DR-bar masses to SLHA convention

   // interface functions
   virtual void calculate_spectrum();
   virtual void print(std::ostream&) const;

   double get_MVG_pole_slha() const { return physical_slha.MVG; }
   double get_MGlu_pole_slha() const { return physical_slha.MGlu; }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return physical_slha.MFv; }
   double get_MFv_pole_slha(int i) const { return physical_slha.MFv(i); }
   double get_MChaP_pole_slha() const { return physical_slha.MChaP; }
   double get_MVP_pole_slha() const { return physical_slha.MVP; }
   double get_MVZ_pole_slha() const { return physical_slha.MVZ; }
   double get_MVZp_pole_slha() const { return physical_slha.MVZp; }
   const Eigen::Array<double,6,1>& get_MSd_pole_slha() const { return physical_slha.MSd; }
   double get_MSd_pole_slha(int i) const { return physical_slha.MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv_pole_slha() const { return physical_slha.MSv; }
   double get_MSv_pole_slha(int i) const { return physical_slha.MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu_pole_slha() const { return physical_slha.MSu; }
   double get_MSu_pole_slha(int i) const { return physical_slha.MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe_pole_slha() const { return physical_slha.MSe; }
   double get_MSe_pole_slha(int i) const { return physical_slha.MSe(i); }
   const Eigen::Array<double,6,1>& get_MSDX_pole_slha() const { return physical_slha.MSDX; }
   double get_MSDX_pole_slha(int i) const { return physical_slha.MSDX(i); }
   const Eigen::Array<double,5,1>& get_Mhh_pole_slha() const { return physical_slha.Mhh; }
   double get_Mhh_pole_slha(int i) const { return physical_slha.Mhh(i); }
   const Eigen::Array<double,5,1>& get_MAh_pole_slha() const { return physical_slha.MAh; }
   double get_MAh_pole_slha(int i) const { return physical_slha.MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm_pole_slha() const { return physical_slha.MHpm; }
   double get_MHpm_pole_slha(int i) const { return physical_slha.MHpm(i); }
   const Eigen::Array<double,8,1>& get_MChi_pole_slha() const { return physical_slha.MChi; }
   double get_MChi_pole_slha(int i) const { return physical_slha.MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha_pole_slha() const { return physical_slha.MCha; }
   double get_MCha_pole_slha(int i) const { return physical_slha.MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return physical_slha.MFe; }
   double get_MFe_pole_slha(int i) const { return physical_slha.MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return physical_slha.MFd; }
   double get_MFd_pole_slha(int i) const { return physical_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return physical_slha.MFu; }
   double get_MFu_pole_slha(int i) const { return physical_slha.MFu(i); }
   const Eigen::Array<double,3,1>& get_MFDX_pole_slha() const { return physical_slha.MFDX; }
   double get_MFDX_pole_slha(int i) const { return physical_slha.MFDX(i); }
   const Eigen::Array<double,7,1>& get_MSHI0_pole_slha() const { return physical_slha.MSHI0; }
   double get_MSHI0_pole_slha(int i) const { return physical_slha.MSHI0(i); }
   const Eigen::Array<double,4,1>& get_MSHIPM_pole_slha() const { return physical_slha.MSHIPM; }
   double get_MSHIPM_pole_slha(int i) const { return physical_slha.MSHIPM(i); }
   const Eigen::Array<double,2,1>& get_MChaI_pole_slha() const { return physical_slha.MChaI; }
   double get_MChaI_pole_slha(int i) const { return physical_slha.MChaI(i); }
   const Eigen::Array<double,7,1>& get_MChiI_pole_slha() const { return physical_slha.MChiI; }
   double get_MChiI_pole_slha(int i) const { return physical_slha.MChiI(i); }
   const Eigen::Array<double,2,1>& get_MSHp0_pole_slha() const { return physical_slha.MSHp0; }
   double get_MSHp0_pole_slha(int i) const { return physical_slha.MSHp0(i); }
   const Eigen::Array<double,2,1>& get_MSHpp_pole_slha() const { return physical_slha.MSHpp; }
   double get_MSHpp_pole_slha(int i) const { return physical_slha.MSHpp(i); }
   const Eigen::Array<double,2,1>& get_MChiP_pole_slha() const { return physical_slha.MChiP; }
   double get_MChiP_pole_slha(int i) const { return physical_slha.MChiP(i); }
   double get_MVWm_pole_slha() const { return physical_slha.MVWm; }

   const Eigen::Matrix<double,6,6>& get_ZD_pole_slha() const { return physical_slha.ZD; }
   double get_ZD_pole_slha(int i, int k) const { return physical_slha.ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV_pole_slha() const { return physical_slha.ZV; }
   double get_ZV_pole_slha(int i, int k) const { return physical_slha.ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU_pole_slha() const { return physical_slha.ZU; }
   double get_ZU_pole_slha(int i, int k) const { return physical_slha.ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE_pole_slha() const { return physical_slha.ZE; }
   double get_ZE_pole_slha(int i, int k) const { return physical_slha.ZE(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZDX_pole_slha() const { return physical_slha.ZDX; }
   double get_ZDX_pole_slha(int i, int k) const { return physical_slha.ZDX(i,k); }
   const Eigen::Matrix<double,5,5>& get_ZH_pole_slha() const { return physical_slha.ZH; }
   double get_ZH_pole_slha(int i, int k) const { return physical_slha.ZH(i,k); }
   const Eigen::Matrix<double,5,5>& get_ZA_pole_slha() const { return physical_slha.ZA; }
   double get_ZA_pole_slha(int i, int k) const { return physical_slha.ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP_pole_slha() const { return physical_slha.ZP; }
   double get_ZP_pole_slha(int i, int k) const { return physical_slha.ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,8,8>& get_ZN_pole_slha() const { return physical_slha.ZN; }
   const std::complex<double>& get_ZN_pole_slha(int i, int k) const { return physical_slha.ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM_pole_slha() const { return physical_slha.UM; }
   const std::complex<double>& get_UM_pole_slha(int i, int k) const { return physical_slha.UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP_pole_slha() const { return physical_slha.UP; }
   const std::complex<double>& get_UP_pole_slha(int i, int k) const { return physical_slha.UP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL_pole_slha() const { return physical_slha.ZEL; }
   const std::complex<double>& get_ZEL_pole_slha(int i, int k) const { return physical_slha.ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER_pole_slha() const { return physical_slha.ZER; }
   const std::complex<double>& get_ZER_pole_slha(int i, int k) const { return physical_slha.ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL_pole_slha() const { return physical_slha.ZDL; }
   const std::complex<double>& get_ZDL_pole_slha(int i, int k) const { return physical_slha.ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR_pole_slha() const { return physical_slha.ZDR; }
   const std::complex<double>& get_ZDR_pole_slha(int i, int k) const { return physical_slha.ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL_pole_slha() const { return physical_slha.ZUL; }
   const std::complex<double>& get_ZUL_pole_slha(int i, int k) const { return physical_slha.ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR_pole_slha() const { return physical_slha.ZUR; }
   const std::complex<double>& get_ZUR_pole_slha(int i, int k) const { return physical_slha.ZUR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXL_pole_slha() const { return physical_slha.ZDXL; }
   const std::complex<double>& get_ZDXL_pole_slha(int i, int k) const { return physical_slha.ZDXL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXR_pole_slha() const { return physical_slha.ZDXR; }
   const std::complex<double>& get_ZDXR_pole_slha(int i, int k) const { return physical_slha.ZDXR(i,k); }
   const Eigen::Matrix<double,7,7>& get_UHI0_pole_slha() const { return physical_slha.UHI0; }
   double get_UHI0_pole_slha(int i, int k) const { return physical_slha.UHI0(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHIPM_pole_slha() const { return physical_slha.UHIPM; }
   double get_UHIPM_pole_slha(int i, int k) const { return physical_slha.UHIPM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZMI_pole_slha() const { return physical_slha.ZMI; }
   const std::complex<double>& get_ZMI_pole_slha(int i, int k) const { return physical_slha.ZMI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZPI_pole_slha() const { return physical_slha.ZPI; }
   const std::complex<double>& get_ZPI_pole_slha(int i, int k) const { return physical_slha.ZPI(i,k); }
   const Eigen::Matrix<std::complex<double>,7,7>& get_ZNI_pole_slha() const { return physical_slha.ZNI; }
   const std::complex<double>& get_ZNI_pole_slha(int i, int k) const { return physical_slha.ZNI(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHp0_pole_slha() const { return physical_slha.UHp0; }
   double get_UHp0_pole_slha(int i, int k) const { return physical_slha.UHp0(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHpp_pole_slha() const { return physical_slha.UHpp; }
   double get_UHpp_pole_slha(int i, int k) const { return physical_slha.UHpp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZNp_pole_slha() const { return physical_slha.ZNp; }
   const std::complex<double>& get_ZNp_pole_slha(int i, int k) const { return physical_slha.ZNp(i,k); }

   double get_MVG_drbar_slha() const { return drbar_slha.MVG; }
   double get_MGlu_drbar_slha() const { return drbar_slha.MGlu; }
   const Eigen::Array<double,3,1>& get_MFv_drbar_slha() const { return drbar_slha.MFv; }
   double get_MFv_drbar_slha(int i) const { return drbar_slha.MFv(i); }
   double get_MChaP_drbar_slha() const { return drbar_slha.MChaP; }
   double get_MVP_drbar_slha() const { return drbar_slha.MVP; }
   double get_MVZ_drbar_slha() const { return drbar_slha.MVZ; }
   double get_MVZp_drbar_slha() const { return drbar_slha.MVZp; }
   const Eigen::Array<double,6,1>& get_MSd_drbar_slha() const { return drbar_slha.MSd; }
   double get_MSd_drbar_slha(int i) const { return drbar_slha.MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv_drbar_slha() const { return drbar_slha.MSv; }
   double get_MSv_drbar_slha(int i) const { return drbar_slha.MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu_drbar_slha() const { return drbar_slha.MSu; }
   double get_MSu_drbar_slha(int i) const { return drbar_slha.MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe_drbar_slha() const { return drbar_slha.MSe; }
   double get_MSe_drbar_slha(int i) const { return drbar_slha.MSe(i); }
   const Eigen::Array<double,6,1>& get_MSDX_drbar_slha() const { return drbar_slha.MSDX; }
   double get_MSDX_drbar_slha(int i) const { return drbar_slha.MSDX(i); }
   const Eigen::Array<double,5,1>& get_Mhh_drbar_slha() const { return drbar_slha.Mhh; }
   double get_Mhh_drbar_slha(int i) const { return drbar_slha.Mhh(i); }
   const Eigen::Array<double,5,1>& get_MAh_drbar_slha() const { return drbar_slha.MAh; }
   double get_MAh_drbar_slha(int i) const { return drbar_slha.MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm_drbar_slha() const { return drbar_slha.MHpm; }
   double get_MHpm_drbar_slha(int i) const { return drbar_slha.MHpm(i); }
   const Eigen::Array<double,8,1>& get_MChi_drbar_slha() const { return drbar_slha.MChi; }
   double get_MChi_drbar_slha(int i) const { return drbar_slha.MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha_drbar_slha() const { return drbar_slha.MCha; }
   double get_MCha_drbar_slha(int i) const { return drbar_slha.MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe_drbar_slha() const { return drbar_slha.MFe; }
   double get_MFe_drbar_slha(int i) const { return drbar_slha.MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd_drbar_slha() const { return drbar_slha.MFd; }
   double get_MFd_drbar_slha(int i) const { return drbar_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_drbar_slha() const { return drbar_slha.MFu; }
   double get_MFu_drbar_slha(int i) const { return drbar_slha.MFu(i); }
   const Eigen::Array<double,3,1>& get_MFDX_drbar_slha() const { return drbar_slha.MFDX; }
   double get_MFDX_drbar_slha(int i) const { return drbar_slha.MFDX(i); }
   const Eigen::Array<double,7,1>& get_MSHI0_drbar_slha() const { return drbar_slha.MSHI0; }
   double get_MSHI0_drbar_slha(int i) const { return drbar_slha.MSHI0(i); }
   const Eigen::Array<double,4,1>& get_MSHIPM_drbar_slha() const { return drbar_slha.MSHIPM; }
   double get_MSHIPM_drbar_slha(int i) const { return drbar_slha.MSHIPM(i); }
   const Eigen::Array<double,2,1>& get_MChaI_drbar_slha() const { return drbar_slha.MChaI; }
   double get_MChaI_drbar_slha(int i) const { return drbar_slha.MChaI(i); }
   const Eigen::Array<double,7,1>& get_MChiI_drbar_slha() const { return drbar_slha.MChiI; }
   double get_MChiI_drbar_slha(int i) const { return drbar_slha.MChiI(i); }
   const Eigen::Array<double,2,1>& get_MSHp0_drbar_slha() const { return drbar_slha.MSHp0; }
   double get_MSHp0_drbar_slha(int i) const { return drbar_slha.MSHp0(i); }
   const Eigen::Array<double,2,1>& get_MSHpp_drbar_slha() const { return drbar_slha.MSHpp; }
   double get_MSHpp_drbar_slha(int i) const { return drbar_slha.MSHpp(i); }
   const Eigen::Array<double,2,1>& get_MChiP_drbar_slha() const { return drbar_slha.MChiP; }
   double get_MChiP_drbar_slha(int i) const { return drbar_slha.MChiP(i); }
   double get_MVWm_drbar_slha() const { return drbar_slha.MVWm; }

   const Eigen::Matrix<double,6,6>& get_ZD_drbar_slha() const { return drbar_slha.ZD; }
   double get_ZD_drbar_slha(int i, int k) const { return drbar_slha.ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV_drbar_slha() const { return drbar_slha.ZV; }
   double get_ZV_drbar_slha(int i, int k) const { return drbar_slha.ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU_drbar_slha() const { return drbar_slha.ZU; }
   double get_ZU_drbar_slha(int i, int k) const { return drbar_slha.ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE_drbar_slha() const { return drbar_slha.ZE; }
   double get_ZE_drbar_slha(int i, int k) const { return drbar_slha.ZE(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZDX_drbar_slha() const { return drbar_slha.ZDX; }
   double get_ZDX_drbar_slha(int i, int k) const { return drbar_slha.ZDX(i,k); }
   const Eigen::Matrix<double,5,5>& get_ZH_drbar_slha() const { return drbar_slha.ZH; }
   double get_ZH_drbar_slha(int i, int k) const { return drbar_slha.ZH(i,k); }
   const Eigen::Matrix<double,5,5>& get_ZA_drbar_slha() const { return drbar_slha.ZA; }
   double get_ZA_drbar_slha(int i, int k) const { return drbar_slha.ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP_drbar_slha() const { return drbar_slha.ZP; }
   double get_ZP_drbar_slha(int i, int k) const { return drbar_slha.ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,8,8>& get_ZN_drbar_slha() const { return drbar_slha.ZN; }
   const std::complex<double>& get_ZN_drbar_slha(int i, int k) const { return drbar_slha.ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM_drbar_slha() const { return drbar_slha.UM; }
   const std::complex<double>& get_UM_drbar_slha(int i, int k) const { return drbar_slha.UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP_drbar_slha() const { return drbar_slha.UP; }
   const std::complex<double>& get_UP_drbar_slha(int i, int k) const { return drbar_slha.UP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL_drbar_slha() const { return drbar_slha.ZEL; }
   const std::complex<double>& get_ZEL_drbar_slha(int i, int k) const { return drbar_slha.ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER_drbar_slha() const { return drbar_slha.ZER; }
   const std::complex<double>& get_ZER_drbar_slha(int i, int k) const { return drbar_slha.ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL_drbar_slha() const { return drbar_slha.ZDL; }
   const std::complex<double>& get_ZDL_drbar_slha(int i, int k) const { return drbar_slha.ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR_drbar_slha() const { return drbar_slha.ZDR; }
   const std::complex<double>& get_ZDR_drbar_slha(int i, int k) const { return drbar_slha.ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL_drbar_slha() const { return drbar_slha.ZUL; }
   const std::complex<double>& get_ZUL_drbar_slha(int i, int k) const { return drbar_slha.ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR_drbar_slha() const { return drbar_slha.ZUR; }
   const std::complex<double>& get_ZUR_drbar_slha(int i, int k) const { return drbar_slha.ZUR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXL_drbar_slha() const { return drbar_slha.ZDXL; }
   const std::complex<double>& get_ZDXL_drbar_slha(int i, int k) const { return drbar_slha.ZDXL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXR_drbar_slha() const { return drbar_slha.ZDXR; }
   const std::complex<double>& get_ZDXR_drbar_slha(int i, int k) const { return drbar_slha.ZDXR(i,k); }
   const Eigen::Matrix<double,7,7>& get_UHI0_drbar_slha() const { return drbar_slha.UHI0; }
   double get_UHI0_drbar_slha(int i, int k) const { return drbar_slha.UHI0(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHIPM_drbar_slha() const { return drbar_slha.UHIPM; }
   double get_UHIPM_drbar_slha(int i, int k) const { return drbar_slha.UHIPM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZMI_drbar_slha() const { return drbar_slha.ZMI; }
   const std::complex<double>& get_ZMI_drbar_slha(int i, int k) const { return drbar_slha.ZMI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZPI_drbar_slha() const { return drbar_slha.ZPI; }
   const std::complex<double>& get_ZPI_drbar_slha(int i, int k) const { return drbar_slha.ZPI(i,k); }
   const Eigen::Matrix<std::complex<double>,7,7>& get_ZNI_drbar_slha() const { return drbar_slha.ZNI; }
   const std::complex<double>& get_ZNI_drbar_slha(int i, int k) const { return drbar_slha.ZNI(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHp0_drbar_slha() const { return drbar_slha.UHp0; }
   double get_UHp0_drbar_slha(int i, int k) const { return drbar_slha.UHp0(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHpp_drbar_slha() const { return drbar_slha.UHpp; }
   double get_UHpp_drbar_slha(int i, int k) const { return drbar_slha.UHpp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZNp_drbar_slha() const { return drbar_slha.ZNp; }
   const std::complex<double>& get_ZNp_drbar_slha(int i, int k) const { return drbar_slha.ZNp(i,k); }

   const Eigen::Array<double,3,1>& get_Yu_slha() const { return Yu_slha; }
   double get_Yu_slha(int i) const { return Yu_slha(i); }
   const Eigen::Array<double,3,1>& get_Yd_slha() const { return Yd_slha; }
   double get_Yd_slha(int i) const { return Yd_slha(i); }
   const Eigen::Array<double,3,1>& get_Ye_slha() const { return Ye_slha; }
   double get_Ye_slha(int i) const { return Ye_slha(i); }

   const Eigen::Matrix<double,3,3>& get_TYu_slha() const { return TYu_slha; }
   double get_TYu_slha(int i, int k) const { return TYu_slha(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYd_slha() const { return TYd_slha; }
   double get_TYd_slha(int i, int k) const { return TYd_slha(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYe_slha() const { return TYe_slha; }
   double get_TYe_slha(int i, int k) const { return TYe_slha(i,k); }

   const Eigen::Matrix<double,3,3>& get_mq2_slha() const { return mq2_slha; }
   double get_mq2_slha(int i, int k) const { return mq2_slha(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2_slha() const { return mu2_slha; }
   double get_mu2_slha(int i, int k) const { return mu2_slha(i,k); }
   const Eigen::Matrix<double,3,3>& get_md2_slha() const { return md2_slha; }
   double get_md2_slha(int i, int k) const { return md2_slha(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2_slha() const { return ml2_slha; }
   double get_ml2_slha(int i, int k) const { return ml2_slha(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2_slha() const { return me2_slha; }
   double get_me2_slha(int i, int k) const { return me2_slha(i,k); }

   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL_slha() const { return ZDL_slha; }
   const std::complex<double>& get_ZDL_slha(int i, int k) const { return ZDL_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL_slha() const { return ZUL_slha; }
   const std::complex<double>& get_ZUL_slha(int i, int k) const { return ZUL_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR_slha() const { return ZDR_slha; }
   const std::complex<double>& get_ZDR_slha(int i, int k) const { return ZDR_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR_slha() const { return ZUR_slha; }
   const std::complex<double>& get_ZUR_slha(int i, int k) const { return ZUR_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL_slha() const { return ZEL_slha; }
   const std::complex<double>& get_ZEL_slha(int i, int k) const { return ZEL_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER_slha() const { return ZER_slha; }
   const std::complex<double>& get_ZER_slha(int i, int k) const { return ZER_slha(i,k); }


private:
   CSE6SSM_physical physical_slha; ///< contains the pole masses and mixings in slha convention
   CSE6SSM_physical drbar_slha;    ///< contains the DR-bar masses and mixings in slha convention
   Eigen::Matrix<std::complex<double>,3,3> ckm;
   Eigen::Matrix<std::complex<double>,3,3> pmns;
   Eigen::Array<double,3,1> Yu_slha;
   Eigen::Array<double,3,1> Yd_slha;
   Eigen::Array<double,3,1> Ye_slha;

   Eigen::Matrix<std::complex<double>,3,3> ZDL_slha;
   Eigen::Matrix<std::complex<double>,3,3> ZUL_slha;
   Eigen::Matrix<std::complex<double>,3,3> ZDR_slha;
   Eigen::Matrix<std::complex<double>,3,3> ZUR_slha;
   Eigen::Matrix<std::complex<double>,3,3> ZEL_slha;
   Eigen::Matrix<std::complex<double>,3,3> ZER_slha;

   Eigen::Matrix<double,3,3> TYu_slha;
   Eigen::Matrix<double,3,3> TYd_slha;
   Eigen::Matrix<double,3,3> TYe_slha;

   Eigen::Matrix<double,3,3> mq2_slha;
   Eigen::Matrix<double,3,3> mu2_slha;
   Eigen::Matrix<double,3,3> md2_slha;
   Eigen::Matrix<double,3,3> ml2_slha;
   Eigen::Matrix<double,3,3> me2_slha;


   void calculate_ckm_matrix();
   void calculate_pmns_matrix();
   void convert_yukawa_couplings_to_slha();
   void convert_trilinear_couplings_to_slha();
   void convert_soft_squared_masses_to_slha();
};

} // namespace flexiblesusy

#endif
