
/**
 * @file CMSSM_semi_two_scale_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

#ifndef CMSSM_SEMI_TWO_SCALE_SLHA_H
#define CMSSM_SEMI_TWO_SCALE_SLHA_H

#include "CMSSM_semi_two_scale_model.hpp"
#include "MSSM_physical.hpp"
#include "CMSSM_semi_model_slha.hpp"

namespace flexiblesusy {

class Two_scale;

/**
 * @class CMSSM_semianalytic_slha<Two_scale>
 * @brief model class wrapper in SLHA convention
 */

template<>
class CMSSM_semianalytic_slha<Two_scale> : public CMSSM_semianalytic<Two_scale> {
public:
   explicit CMSSM_semianalytic_slha(const CMSSM_semianalytic_input_parameters<Two_scale>& input_ = CMSSM_semianalytic_input_parameters<Two_scale>());
   explicit CMSSM_semianalytic_slha(const CMSSM_semianalytic<Two_scale>&);
   virtual ~CMSSM_semianalytic_slha();

   virtual void clear();
   void convert_to_slha(); ///< converts pole masses and couplings to SLHA convention
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm_matrix() const { return ckm; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns_matrix() const { return pmns; }
   const MSSM_physical& get_physical_slha() const; ///< return pole masses to SLHA convention
   MSSM_physical& get_physical_slha(); ///< return pole masses to SLHA convention
   const MSSM_physical& get_drbar_slha() const; ///< returns DR-bar masses to SLHA convention
   MSSM_physical& get_drbar_slha(); ///< returns DR-bar masses to SLHA convention

   // interface functions
   virtual void calculate_spectrum();
   virtual void print(std::ostream&) const;

   double get_MVG_pole_slha() const { return physical_slha.MVG; }
   double get_MGlu_pole_slha() const { return physical_slha.MGlu; }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return physical_slha.MFv; }
   double get_MFv_pole_slha(int i) const { return physical_slha.MFv(i); }
   double get_MVP_pole_slha() const { return physical_slha.MVP; }
   double get_MVZ_pole_slha() const { return physical_slha.MVZ; }
   const Eigen::Array<double,6,1>& get_MSd_pole_slha() const { return physical_slha.MSd; }
   double get_MSd_pole_slha(int i) const { return physical_slha.MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv_pole_slha() const { return physical_slha.MSv; }
   double get_MSv_pole_slha(int i) const { return physical_slha.MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu_pole_slha() const { return physical_slha.MSu; }
   double get_MSu_pole_slha(int i) const { return physical_slha.MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe_pole_slha() const { return physical_slha.MSe; }
   double get_MSe_pole_slha(int i) const { return physical_slha.MSe(i); }
   const Eigen::Array<double,2,1>& get_Mhh_pole_slha() const { return physical_slha.Mhh; }
   double get_Mhh_pole_slha(int i) const { return physical_slha.Mhh(i); }
   const Eigen::Array<double,2,1>& get_MAh_pole_slha() const { return physical_slha.MAh; }
   double get_MAh_pole_slha(int i) const { return physical_slha.MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm_pole_slha() const { return physical_slha.MHpm; }
   double get_MHpm_pole_slha(int i) const { return physical_slha.MHpm(i); }
   const Eigen::Array<double,4,1>& get_MChi_pole_slha() const { return physical_slha.MChi; }
   double get_MChi_pole_slha(int i) const { return physical_slha.MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha_pole_slha() const { return physical_slha.MCha; }
   double get_MCha_pole_slha(int i) const { return physical_slha.MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return physical_slha.MFe; }
   double get_MFe_pole_slha(int i) const { return physical_slha.MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return physical_slha.MFd; }
   double get_MFd_pole_slha(int i) const { return physical_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return physical_slha.MFu; }
   double get_MFu_pole_slha(int i) const { return physical_slha.MFu(i); }
   double get_MVWm_pole_slha() const { return physical_slha.MVWm; }

   const Eigen::Matrix<double,6,6>& get_ZD_pole_slha() const { return physical_slha.ZD; }
   double get_ZD_pole_slha(int i, int k) const { return physical_slha.ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV_pole_slha() const { return physical_slha.ZV; }
   double get_ZV_pole_slha(int i, int k) const { return physical_slha.ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU_pole_slha() const { return physical_slha.ZU; }
   double get_ZU_pole_slha(int i, int k) const { return physical_slha.ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE_pole_slha() const { return physical_slha.ZE; }
   double get_ZE_pole_slha(int i, int k) const { return physical_slha.ZE(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZH_pole_slha() const { return physical_slha.ZH; }
   double get_ZH_pole_slha(int i, int k) const { return physical_slha.ZH(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZA_pole_slha() const { return physical_slha.ZA; }
   double get_ZA_pole_slha(int i, int k) const { return physical_slha.ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP_pole_slha() const { return physical_slha.ZP; }
   double get_ZP_pole_slha(int i, int k) const { return physical_slha.ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN_pole_slha() const { return physical_slha.ZN; }
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

   double get_MVG_drbar_slha() const { return drbar_slha.MVG; }
   double get_MGlu_drbar_slha() const { return drbar_slha.MGlu; }
   const Eigen::Array<double,3,1>& get_MFv_drbar_slha() const { return drbar_slha.MFv; }
   double get_MFv_drbar_slha(int i) const { return drbar_slha.MFv(i); }
   double get_MVP_drbar_slha() const { return drbar_slha.MVP; }
   double get_MVZ_drbar_slha() const { return drbar_slha.MVZ; }
   const Eigen::Array<double,6,1>& get_MSd_drbar_slha() const { return drbar_slha.MSd; }
   double get_MSd_drbar_slha(int i) const { return drbar_slha.MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv_drbar_slha() const { return drbar_slha.MSv; }
   double get_MSv_drbar_slha(int i) const { return drbar_slha.MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu_drbar_slha() const { return drbar_slha.MSu; }
   double get_MSu_drbar_slha(int i) const { return drbar_slha.MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe_drbar_slha() const { return drbar_slha.MSe; }
   double get_MSe_drbar_slha(int i) const { return drbar_slha.MSe(i); }
   const Eigen::Array<double,2,1>& get_Mhh_drbar_slha() const { return drbar_slha.Mhh; }
   double get_Mhh_drbar_slha(int i) const { return drbar_slha.Mhh(i); }
   const Eigen::Array<double,2,1>& get_MAh_drbar_slha() const { return drbar_slha.MAh; }
   double get_MAh_drbar_slha(int i) const { return drbar_slha.MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm_drbar_slha() const { return drbar_slha.MHpm; }
   double get_MHpm_drbar_slha(int i) const { return drbar_slha.MHpm(i); }
   const Eigen::Array<double,4,1>& get_MChi_drbar_slha() const { return drbar_slha.MChi; }
   double get_MChi_drbar_slha(int i) const { return drbar_slha.MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha_drbar_slha() const { return drbar_slha.MCha; }
   double get_MCha_drbar_slha(int i) const { return drbar_slha.MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe_drbar_slha() const { return drbar_slha.MFe; }
   double get_MFe_drbar_slha(int i) const { return drbar_slha.MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd_drbar_slha() const { return drbar_slha.MFd; }
   double get_MFd_drbar_slha(int i) const { return drbar_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_drbar_slha() const { return drbar_slha.MFu; }
   double get_MFu_drbar_slha(int i) const { return drbar_slha.MFu(i); }
   double get_MVWm_drbar_slha() const { return drbar_slha.MVWm; }

   const Eigen::Matrix<double,6,6>& get_ZD_drbar_slha() const { return drbar_slha.ZD; }
   double get_ZD_drbar_slha(int i, int k) const { return drbar_slha.ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV_drbar_slha() const { return drbar_slha.ZV; }
   double get_ZV_drbar_slha(int i, int k) const { return drbar_slha.ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU_drbar_slha() const { return drbar_slha.ZU; }
   double get_ZU_drbar_slha(int i, int k) const { return drbar_slha.ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE_drbar_slha() const { return drbar_slha.ZE; }
   double get_ZE_drbar_slha(int i, int k) const { return drbar_slha.ZE(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZH_drbar_slha() const { return drbar_slha.ZH; }
   double get_ZH_drbar_slha(int i, int k) const { return drbar_slha.ZH(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZA_drbar_slha() const { return drbar_slha.ZA; }
   double get_ZA_drbar_slha(int i, int k) const { return drbar_slha.ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP_drbar_slha() const { return drbar_slha.ZP; }
   double get_ZP_drbar_slha(int i, int k) const { return drbar_slha.ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN_drbar_slha() const { return drbar_slha.ZN; }
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
   MSSM_physical physical_slha; ///< contains the pole masses and mixings in slha convention
   MSSM_physical drbar_slha;    ///< contains the DR-bar masses and mixings in slha convention
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
