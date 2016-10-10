
/**
 * @file CMSSM_semi_two_scale_model_slha.cpp
 * @brief CMSSM model class wrapper for SLHA convention
 */

#include "CMSSM_semi_two_scale_model_slha.hpp"
#include "slha_io.hpp"
#include "ckm.hpp"
#include "pmns.hpp"

namespace flexiblesusy {

#define CLASSNAME CMSSM_semianalytic_slha<Two_scale>
#define LOCALPHYSICAL(p) physical.p

CLASSNAME::CMSSM_semianalytic_slha(const CMSSM_semianalytic_input_parameters<Two_scale>& input_)
   : CMSSM_semianalytic<Two_scale>(input_)
   , physical_slha()
   , drbar_slha()
   , ckm(Eigen::Matrix<std::complex<double>,3,3>::Identity())
   , pmns(Eigen::Matrix<std::complex<double>,3,3>::Identity())
{
}

/**
 * Copy constructor.  Copies from base class (two-scale model class in
 * BPMZ convention) and converts parameters to SLHA.
 *
 * @param model_ model class in BPMZ convention
 */
CLASSNAME::CMSSM_semianalytic_slha(const CMSSM_semianalytic<Two_scale>& model_)
   : CMSSM_semianalytic<Two_scale>(model_)
{
   convert_to_slha();
}

CLASSNAME::~CMSSM_semianalytic_slha()
{
}

void CLASSNAME::clear()
{
   CMSSM_semianalytic<Two_scale>::clear();
   physical_slha.clear();
   drbar_slha.clear();
}

void CLASSNAME::calculate_spectrum()
{
   CMSSM_semianalytic<Two_scale>::calculate_spectrum();
   convert_to_slha();
}

void CLASSNAME::convert_to_slha()
{
   physical_slha = get_physical();
   physical_slha.convert_to_slha();

   drbar_slha = get_drbar_masses();
   drbar_slha.convert_to_slha();

   convert_yukawa_couplings_to_slha();
   calculate_ckm_matrix();
   calculate_pmns_matrix();
   convert_trilinear_couplings_to_slha();
   convert_soft_squared_masses_to_slha();
}

void CLASSNAME::calculate_ckm_matrix()
{
   ckm = ZUL_slha * ZDL_slha.adjoint();
   CKM_parameters::to_pdg_convention(ckm, ZUL_slha, ZDL_slha, ZUR_slha, ZDR_slha);

}

void CLASSNAME::calculate_pmns_matrix()
{
   pmns << 1, 0, 0, 0, 1, 0, 0, 0, 1;

}

/**
 * Convert Yukawa couplings to SLHA convention
 */
void CLASSNAME::convert_yukawa_couplings_to_slha()
{
   fs_svd(Yu, Yu_slha, ZUR_slha, ZUL_slha);
   fs_svd(Yd, Yd_slha, ZDR_slha, ZDL_slha);
   fs_svd(Ye, Ye_slha, ZER_slha, ZEL_slha);

}

/**
 * Convert trilinear couplings to SLHA convention
 */
void CLASSNAME::convert_trilinear_couplings_to_slha()
{
   TYu_slha = (ZUR_slha.conjugate() * TYu * ZUL_slha.adjoint()).real();
   TYd_slha = (ZDR_slha.conjugate() * TYd * ZDL_slha.adjoint()).real();
   TYe_slha = (ZER_slha.conjugate() * TYe * ZEL_slha.adjoint()).real();

}

/**
 * Convert soft squared masses to SLHA convention
 */
void CLASSNAME::convert_soft_squared_masses_to_slha()
{
   mq2_slha = (ZDL_slha * mq2 * ZDL_slha.adjoint()).real();
   mu2_slha = (ZUR_slha.conjugate() * mu2 * ZUR_slha.transpose()).real();
   md2_slha = (ZDR_slha.conjugate() * md2 * ZDR_slha.transpose()).real();
   ml2_slha = (ZEL_slha * ml2 * ZEL_slha.adjoint()).real();
   me2_slha = (ZER_slha.conjugate() * me2 * ZER_slha.transpose()).real();

}

const MSSM_physical& CLASSNAME::get_physical_slha() const
{
   return physical_slha;
}

MSSM_physical& CLASSNAME::get_physical_slha()
{
   return physical_slha;
}

const MSSM_physical& CLASSNAME::get_drbar_slha() const
{
   return drbar_slha;
}

MSSM_physical& CLASSNAME::get_drbar_slha()
{
   return drbar_slha;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   CMSSM_semianalytic<Two_scale>::print(ostr);

   ostr << "----------------------------------------\n"
           "SLHA convention:\n"
           "----------------------------------------\n";
   ostr << "----------------------------------------\n"
           "Pole masses:\n"
           "----------------------------------------\n";
   physical_slha.print(ostr);
   ostr << "----------------------------------------\n"
           "DR-bar masses:\n"
           "----------------------------------------\n";
   drbar_slha.print(ostr);
}

} // namespace flexiblesusy
