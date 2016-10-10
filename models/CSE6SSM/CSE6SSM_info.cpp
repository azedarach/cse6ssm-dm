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

// File generated at Wed 3 Jun 2015 23:47:48

#include "CSE6SSM_info.hpp"

#include <iostream>

namespace flexiblesusy {

namespace CSE6SSM_info {
   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;
   const double normalization_g1p = 0.15811388300841897;

   const unsigned particle_multiplicities[NUMBER_OF_PARTICLES] = {1, 1, 3, 1, 1
      , 1, 1, 6, 3, 6, 6, 6, 5, 5, 2, 8, 2, 3, 3, 3, 3, 7, 4, 2, 7, 2, 2, 2, 1};

   const char* particle_names[NUMBER_OF_PARTICLES] = {"VG", "Glu", "Fv", "ChaP"
      , "VP", "VZ", "VZp", "Sd", "Sv", "Su", "Se", "SDX", "hh", "Ah", "Hpm", "Chi"
      , "Cha", "Fe", "Fd", "Fu", "FDX", "SHI0", "SHIPM", "ChaI", "ChiI", "SHp0",
      "SHpp", "ChiP", "VWm"};

   const char* particle_latex_names[NUMBER_OF_PARTICLES] = {   "g",
      "\\tilde{g}", "\\nu", "\\tilde{\\chi'}^{-}", "\\gamma", "Z", "{Z\\ {}^\\prime}",
      "\\tilde{d}", "\\tilde{\\nu}", "\\tilde{u}", "\\tilde{e}", "\\tilde{D}", "h"
      , "A^0", "H^{\\pm}", "\\tilde{\\chi}^0", "\\tilde{\\chi}^-", "e", "d", "u", "D",
      "h_I^{0}", "h_I^{-}", "\\tilde{\\chi}_I^{-}",
      "\\tilde{\\chi}_I^{0}", "L_4^{0}", "L_4^{-}", "\\tilde{\\chi'}^{0}", "W^-"
      };

   const char* parameter_names[NUMBER_OF_PARAMETERS] = {"Yd(0,0)", "Yd(0,1)",
      "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)",
      "hE(0,0)", "hE(0,1)", "hE(1,0)", "hE(1,1)", "hE(2,0)", "hE(2,1)", "Ye(0,0)"
      , "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)",
      "Ye(2,1)", "Ye(2,2)", "SigmaL", "KappaPr", "Sigmax", "gD(0,0)", "gD(0,1)",
      "gD(0,2)", "gD(1,0)", "gD(1,1)", "gD(1,2)", "gD(2,0)", "gD(2,1)", "gD(2,2)",
      "Kappa(0,0)", "Kappa(0,1)", "Kappa(0,2)", "Kappa(1,0)", "Kappa(1,1)",
      "Kappa(1,2)", "Kappa(2,0)", "Kappa(2,1)", "Kappa(2,2)", "Lambda12(0,0)",
      "Lambda12(0,1)", "Lambda12(1,0)", "Lambda12(1,1)", "Lambdax", "fu(0,0)",
      "fu(0,1)", "fu(1,0)", "fu(1,1)", "fu(2,0)", "fu(2,1)", "fd(0,0)", "fd(0,1)",
      "fd(1,0)", "fd(1,1)", "fd(2,0)", "fd(2,1)", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)"
      , "Yu(1,0)", "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "MuPr",
      "MuPhi", "XiF", "g1", "g2", "g3", "g1p", "vd", "vu", "vs", "vsb", "vphi",
      "QS", "TYd(0,0)", "TYd(0,1)", "TYd(0,2)", "TYd(1,0)", "TYd(1,1)", "TYd(1,2)",
      "TYd(2,0)", "TYd(2,1)", "TYd(2,2)", "ThE(0,0)", "ThE(0,1)", "ThE(1,0)",
      "ThE(1,1)", "ThE(2,0)", "ThE(2,1)", "TYe(0,0)", "TYe(0,1)", "TYe(0,2)",
      "TYe(1,0)", "TYe(1,1)", "TYe(1,2)", "TYe(2,0)", "TYe(2,1)", "TYe(2,2)",
      "TSigmaL", "TKappaPr", "TSigmax", "TgD(0,0)", "TgD(0,1)", "TgD(0,2)",
      "TgD(1,0)", "TgD(1,1)", "TgD(1,2)", "TgD(2,0)", "TgD(2,1)", "TgD(2,2)",
      "TKappa(0,0)", "TKappa(0,1)", "TKappa(0,2)", "TKappa(1,0)", "TKappa(1,1)",
      "TKappa(1,2)", "TKappa(2,0)", "TKappa(2,1)", "TKappa(2,2)", "TLambda12(0,0)"
      , "TLambda12(0,1)", "TLambda12(1,0)", "TLambda12(1,1)", "TLambdax",
      "Tfu(0,0)", "Tfu(0,1)", "Tfu(1,0)", "Tfu(1,1)", "Tfu(2,0)", "Tfu(2,1)",
      "Tfd(0,0)", "Tfd(0,1)", "Tfd(1,0)", "Tfd(1,1)", "Tfd(2,0)", "Tfd(2,1)",
      "TYu(0,0)", "TYu(0,1)", "TYu(0,2)", "TYu(1,0)", "TYu(1,1)", "TYu(1,2)",
      "TYu(2,0)", "TYu(2,1)", "TYu(2,2)", "BMuPr", "BMuPhi", "LXiF", "mq2(0,0)",
      "mq2(0,1)", "mq2(0,2)", "mq2(1,0)", "mq2(1,1)", "mq2(1,2)", "mq2(2,0)",
      "mq2(2,1)", "mq2(2,2)", "ml2(0,0)", "ml2(0,1)", "ml2(0,2)", "ml2(1,0)",
      "ml2(1,1)", "ml2(1,2)", "ml2(2,0)", "ml2(2,1)", "ml2(2,2)", "mHd2", "mHu2",
      "md2(0,0)", "md2(0,1)", "md2(0,2)", "md2(1,0)", "md2(1,1)", "md2(1,2)",
      "md2(2,0)", "md2(2,1)", "md2(2,2)", "mu2(0,0)", "mu2(0,1)", "mu2(0,2)",
      "mu2(1,0)", "mu2(1,1)", "mu2(1,2)", "mu2(2,0)", "mu2(2,1)", "mu2(2,2)",
      "me2(0,0)", "me2(0,1)", "me2(0,2)", "me2(1,0)", "me2(1,1)", "me2(1,2)",
      "me2(2,0)", "me2(2,1)", "me2(2,2)", "ms2", "msbar2", "mH1I2(0,0)",
      "mH1I2(0,1)", "mH1I2(1,0)", "mH1I2(1,1)", "mH2I2(0,0)", "mH2I2(0,1)",
      "mH2I2(1,0)", "mH2I2(1,1)", "mSI2(0,0)", "mSI2(0,1)", "mSI2(0,2)",
      "mSI2(1,0)", "mSI2(1,1)", "mSI2(1,2)", "mSI2(2,0)", "mSI2(2,1)", "mSI2(2,2)"
      , "mDx2(0,0)", "mDx2(0,1)", "mDx2(0,2)", "mDx2(1,0)", "mDx2(1,1)",
      "mDx2(1,2)", "mDx2(2,0)", "mDx2(2,1)", "mDx2(2,2)", "mDxbar2(0,0)",
      "mDxbar2(0,1)", "mDxbar2(0,2)", "mDxbar2(1,0)", "mDxbar2(1,1)",
      "mDxbar2(1,2)", "mDxbar2(2,0)", "mDxbar2(2,1)", "mDxbar2(2,2)", "mHp2",
      "mHpbar2", "mphi2", "MassB", "MassWB", "MassG", "MassBp"};

   const int parameter_mass_dimensions[NUMBER_OF_PARAMETERS] = {0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 1, 1, 2, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
      , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 1, 1, 1, 1};

   const char* parameter_latex_names[NUMBER_OF_PARAMETERS] = {"y_{11}^D",
      "y_{12}^D", "y_{13}^D", "y_{21}^D", "y_{22}^D", "y_{23}^D", "y_{31}^D",
      "y_{32}^D", "y_{33}^D", "h_{11}^E", "h_{12}^E", "h_{21}^E", "h_{22}^E",
      "h_{31}^E", "h_{32}^E", "y_{11}^E", "y_{12}^E", "y_{13}^E", "y_{21}^E",
      "y_{22}^E", "y_{23}^E", "y_{31}^E", "y_{32}^E", "y_{33}^E", "\\sigma_L",
      "\\kappa", "\\sigma", "g_{11}^D", "g_{12}^D", "g_{13}^D", "g_{21}^D",
      "g_{22}^D", "g_{23}^D", "g_{31}^D", "g_{32}^D", "g_{33}^D", "\\kappa_{11}",
      "\\kappa_{12}", "\\kappa_{13}", "\\kappa_{21}", "\\kappa_{22}", 
      "\\kappa_{23}", "\\kappa_{31}", "\\kappa_{32}", "\\kappa_{33}", "\\lambda_{11}",
      "\\lambda_{12}", "\\lambda_{21}", "\\lambda_{22}", "\\lambda", "\\tilde{f}_{11}",
      "\\tilde{f}_{12}", "\\tilde{f}_{21}", "\\tilde{f}_{22}", "\\tilde{f}_{31}",
      "\\tilde{f}_{32}", "f_{11}", "f_{12}", "f_{21}", "f_{22}", "f_{31}", "f_{32}",
      "y_{11}^U", "y_{12}^U", "y_{13}^U", "y_{21}^U", "y_{22}^U", "y_{23}^U",
      "y_{31}^U", "y_{32}^U", "y_{33}^U", "\\mu_L", "\\mu", "\\Lambda", "g_1", 
      "g_2", "g_3", "g_1'", "v_1", "v_2", "s_1", "s_2", "\\varphi", "Q_S", "T_{11}^D",
      "T_{12}^D", "T_{13}^D", "T_{21}^D", "T_{22}^D", "T_{23}^D", "T_{31}^D", 
      "T_{32}^D", "T_{33}^D", "T_{11}^h", "T_{12}^h", "T_{21}^h", "T_{22}^h",
      "T_{31}^h", "T_{32}^h", "T_{11}^E", "T_{12}^E", "T_{13}^E", "T_{21}^E",
      "T_{22}^E", "T_{23}^E", "T_{31}^E", "T_{32}^E", "T_{33}^E", "T_{\\sigma_L}",
      "T_{\\kappa}", "T_{\\sigma}", "T_{11}^g", "T_{12}^g", "T_{13}^g", "T_{21}^g",
      "T_{22}^g", "T_{23}^g", "T_{31}^g", "T_{32}^g", "T_{33}^g", "T_{11}^{\\kappa}",
      "T_{12}^{\\kappa}", "T_{13}^{\\kappa}", "T_{21}^{\\kappa}", "T_{22}^{\\kappa}",
      "T_{23}^{\\kappa}", "T_{31}^{\\kappa}", "T_{32}^{\\kappa}", "T_{33}^{\\kappa}",
      "T_{11}^{\\lambda}", "T_{12}^{\\lambda}", "T_{21}^{\\lambda}", "T_{22}^{\\lambda}",
      "T_{\\lambda}", "T_{11}^{\\tilde{f}}", "T_{12}^{\\tilde{f}}", "T_{21}^{\\tilde{f}}",
      "T_{22}^{\\tilde{f}}", "T_{31}^{\\tilde{f}}", "T_{32}^{\\tilde{f}}", "T_{11}^f",
      "T_{12}^f", "T_{21}^f", "T_{22}^f", "T_{31}^f", "T_{32}^f", "T_{11}^U", "T_{12}^U",
      "T_{13}^U", "T_{21}^U", "T_{22}^U", "T_{23}^U", "T_{31}^U", "T_{32}^U", "T_{33}^U",
      "b_{\\mu_L}", "b_{\\mu}", "\\xi\\Lambda", "m_{Q,11}^2", "m_{Q,12}^2", "m_{Q,13}^2",
      "m_{Q,21}^2", "m_{Q,22}^2", "m_{Q,23}^2", "m_{Q,31}^2", "m_{Q,32}^2", "m_{Q,33}^2",
      "m_{L,11}^2", "m_{L,12}^2", "m_{L,13}^2", "m_{L,21}^2", "m_{L,22}^2", "m_{L,23}^2",
      "m_{L,31}^2", "m_{L,32}^2", "m_{L,33}^2", "m_{H_d}^2", "m_{H_u}^2", "m_{d,11}^2",
      "m_{d,12}^2", "m_{d,13}^2", "m_{d,21}^2", "m_{d,22}^2", "m_{d,23}^2", "m_{d,31}^2",
      "m_{d,32}^2", "m_{d,33}^2", "m_{u,11}^2", "m_{u,12}^2", "m_{u,13}^2", "m_{u,21}^2",
      "m_{u,22}^2", "m_{u,23}^2", "m_{u,31}^2", "m_{u,32}^2", "m_{u,33}^2", "m_{e,11}^2",
      "m_{e,12}^2", "m_{e,13}^2", "m_{e,21}^2", "m_{e,22}^2", "m_{e,23}^2", "m_{e,31}^2",
      "m_{e,32}^2", "m_{e,33}^2", "m_S^2", "m_{\\bar{S}}^2", "m_{H_1,11}^2", "m_{H_1,12}^2",
      "m_{H_1,21}^2", "m_{H_1,22}^2", "m_{H_2,11}^2", "m_{H_2,12}^2", "m_{H_2,21}^2",
      "m_{H_2,22}^2", "m_{S,11}^2", "m_{S,12}^2", "m_{S,13}^2", "m_{S,21}^2", "m_{S,22}^2",
      "m_{S,23}^2", "m_{S,31}^2", "m_{S,32}^2", "m_{S,33}^2", "m_{D,11}^2", "m_{D,12}^2",
      "m_{D,13}^2", "m_{D,21}^2", "m_{D,22}^2", "m_{D,23}^2", "m_{D,31}^2", "m_{D,32}^2",
      "m_{D,33}^2", "m_{\\bar{D},11}^2", "m_{\\bar{D},12}^2", "m_{\\bar{D},13}^2",
      "m_{\\bar{D},21}^2", "m_{\\bar{D},22}^2", "m_{\\bar{D},23}^2", "m_{\\bar{D},31}^2",
      "m_{\\bar{D},32}^2", "m_{\\bar{D},33}^2", "m_{L_4}^2", "m_{\\bar{L_4}}^2", 
      "m_{\\phi}^2", "M_1", "M_2", "M_3", "M_1'"};

   namespace two_scale {

      const char* input_names[NUMBER_OF_INPUTS] = {"m0", "m12", "TanBeta",
         "SignLambdax", "Azero", "sInput", "QSInput", "hEInput(0,0)", "hEInput(0,1)",
         "hEInput(1,0)", "hEInput(1,1)", "hEInput(2,0)", "hEInput(2,1)", "SigmaLInput",
         "KappaPrInput", "SigmaxInput", "gDInput(0,0)", "gDInput(0,1)", "gDInput(0,2)",
         "gDInput(1,0)", "gDInput(1,1)", "gDInput(1,2)", "gDInput(2,0)", "gDInput(2,1)",
         "gDInput(2,2)", "KappaInput(0,0)", "KappaInput(0,1)", "KappaInput(0,2)",
         "KappaInput(1,0)", "KappaInput(1,1)", "KappaInput(1,2)", "KappaInput(2,0)",
         "KappaInput(2,1)", "KappaInput(2,2)", "Lambda12Input(0,0)", "Lambda12Input(0,1)",
         "Lambda12Input(1,0)", "Lambda12Input(1,1)", "fuInput(0,0)", "fuInput(0,1)",
         "fuInput(1,0)", "fuInput(1,1)", "fuInput(2,0)", "fuInput(2,1)", "fdInput(0,0)",
         "fdInput(0,1)", "fdInput(1,0)", "fdInput(1,1)", "fdInput(2,0)", "fdInput(2,1)",
         "MuPrInput", "MuPhiInput", "BMuPrInput", "BMuPhiInput"};
      
      const int input_mass_dimensions[NUMBER_OF_INPUTS] = {1, 1, 0, 0, 
         1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2};
      
      const char* input_latex_names[NUMBER_OF_INPUTS] = {"m_0", "m_{1/2}",
         "\\tan\\beta (M_Z)", "\\textrm{sgn}(\\lambda)", "A_0", "s (M_{\\mathrm{SUSY}})",
         "\\tilde{Q}_S", "h_{11}^E(M_{\\mathrm{GUT}})", "h_{12}^E(M_{\\mathrm{GUT}})",
         "h_{21}^E(M_{\\mathrm{GUT}})", "h_{22}^E(M_{\\mathrm{GUT}})",
         "h_{31}^E(M_{\\mathrm{GUT}})", "h_{32}^E(M_{\\mathrm{GUT}})",
         "\\sigma_L(M_{\\mathrm{GUT}})", "\\kappa(M_{\\mathrm{GUT}})",
         "\\sigma(M_{\\mathrm{GUT}})", "g_{11}^D(M_{\\mathrm{GUT}})",
         "g_{12}^D(M_{\\mathrm{GUT}})", "g_{13}^D(M_{\\mathrm{GUT}})",
         "g_{21}^D(M_{\\mathrm{GUT}})", "g_{22}^D(M_{\\mathrm{GUT}})",
         "g_{23}^D(M_{\\mathrm{GUT}})", "g_{31}^D(M_{\\mathrm{GUT}})",
         "g_{32}^D(M_{\\mathrm{GUT}})", "g_{33}^D(M_{\\mathrm{GUT}})",
         "\\kappa_{11}(M_{\\mathrm{GUT}})", "\\kappa_{12}(M_{\\mathrm{GUT}})",
         "\\kappa_{13}(M_{\\mathrm{GUT}})", "\\kappa_{21}(M_{\\mathrm{GUT}})",
         "\\kappa_{22}(M_{\\mathrm{GUT}})", "\\kappa_{23}(M_{\\mathrm{GUT}})",
         "\\kappa_{31}(M_{\\mathrm{GUT}})", "\\kappa_{32}(M_{\\mathrm{GUT}})",
         "\\kappa_{33}(M_{\\mathrm{GUT}})", "\\lambda_{11}(M_{\\mathrm{GUT}})",
         "\\lambda_{12}(M_{\\mathrm{GUT}})", "\\lambda_{21}(M_{\\mathrm{GUT}})",
         "\\lambda_{22}(M_{\\mathrm{GUT}})", "\\tilde{f}_{11}(M_{\\mathrm{GUT}})",
         "\\tilde{f}_{12}(M_{\\mathrm{GUT}})", "\\tilde{f}_{21}(M_{\\mathrm{GUT}})",
         "\\tilde{f}_{22}(M_{\\mathrm{GUT}})", "\\tilde{f}_{31}(M_{\\mathrm{GUT}})",
         "\\tilde{f}_{32}(M_{\\mathrm{GUT}})", "f_{11}(M_{\\mathrm{GUT}})",
         "f_{12}(M_{\\mathrm{GUT}})", "f_{21}(M_{\\mathrm{GUT}})",
         "f_{22}(M_{\\mathrm{GUT}})", "f_{31}(M_{\\mathrm{GUT}})",
         "f_{32}(M_{\\mathrm{GUT}})", "\\mu_L(M_{\\mathrm{GUT}})",
         "\\mu(M_{\\mathrm{GUT}})", "B_{\\mu_L}(M_{\\mathrm{GUT}})",
         "B_\\mu(M_{\\mathrm{GUT}})"};

   } // namespace two_scale

   namespace semianalytic {

      const char* input_names[NUMBER_OF_INPUTS] = {"m12", "Azero",
         "TanBeta", "sInput", "QSInput", "hEInput(0,0)", "hEInput(0,1)", "hEInput(1,0)",
         "hEInput(1,1)", "hEInput(2,0)", "hEInput(2,1)", "SigmaLInput", "KappaPrInput",
         "SigmaxInput", "gDInput(0,0)", "gDInput(0,1)", "gDInput(0,2)", "gDInput(1,0)",
         "gDInput(1,1)", "gDInput(1,2)", "gDInput(2,0)", "gDInput(2,1)", "gDInput(2,2)",
         "KappaInput(0,0)", "KappaInput(0,1)", "KappaInput(0,2)", "KappaInput(1,0)",
         "KappaInput(1,1)", "KappaInput(1,2)", "KappaInput(2,0)", "KappaInput(2,1)",
         "KappaInput(2,2)", "Lambda12Input(0,0)", "Lambda12Input(0,1)", "Lambda12Input(1,0)",
         "Lambda12Input(1,1)", "LambdaxInput", "fuInput(0,0)", "fuInput(0,1)", "fuInput(1,0)",
         "fuInput(1,1)", "fuInput(2,0)", "fuInput(2,1)", "fdInput(0,0)", "fdInput(0,1)",
         "fdInput(1,0)", "fdInput(1,1)", "fdInput(2,0)", "fdInput(2,1)", "MuPrInput",
         "MuPhiInput", "BMuPrInput", "BMuPhiInput"};

      const int input_mass_dimensions[NUMBER_OF_INPUTS] = {1, 1,
         0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2};
      
      const char* input_latex_names[NUMBER_OF_INPUTS] = {"m_{1/2}",
         "A_0", "\\tan\\beta (M_Z)", "s(M_{\\mathrm{SUSY}})", "\\tilde{Q}_S",
         "h_{11}^E(M_{\\mathrm{GUT}})", "h_{12}^E(M_{\\mathrm{GUT}})", "h_{21}^E(M_{\\mathrm{GUT}})",
         "h_{22}^E(M_{\\mathrm{GUT}})", "h_{31}^E(M_{\\mathrm{GUT}})", "h_{32}^E(M_{\\mathrm{GUT}})",
         "\\sigma_L(M_{\\mathrm{GUT}})", "\\kappa(M_{\\mathrm{GUT}})", "\\sigma(M_{\\mathrm{GUT}})",
         "g_{11}^D(M_{\\mathrm{GUT}})", "g_{12}^D(M_{\\mathrm{GUT}})", "g_{13}^D(M_{\\mathrm{GUT}})",
         "g_{21}^D(M_{\\mathrm{GUT}})", "g_{22}^D(M_{\\mathrm{GUT}})", "g_{23}^D(M_{\\mathrm{GUT}})",
         "g_{31}^D(M_{\\mathrm{GUT}})", "g_{32}^D(M_{\\mathrm{GUT}})", "g_{33}^D(M_{\\mathrm{GUT}})",
         "\\kappa_{11}(M_{\\mathrm{GUT}})", "\\kappa_{12}(M_{\\mathrm{GUT}})",
         "\\kappa_{13}(M_{\\mathrm{GUT}})", "\\kappa_{21}(M_{\\mathrm{GUT}})",
         "\\kappa_{22}(M_{\\mathrm{GUT}})", "\\kappa_{23}(M_{\\mathrm{GUT}})",
         "\\kappa_{31}(M_{\\mathrm{GUT}})", "\\kappa_{32}(M_{\\mathrm{GUT}})",
         "\\kappa_{33}(M_{\\mathrm{GUT}})", "\\lambda_{11}(M_{\\mathrm{GUT}})",
         "\\lambda_{12}(M_{\\mathrm{GUT}})", "\\lambda_{21}(M_{\\mathrm{GUT}})",
         "\\lambda_{22}(M_{\\mathrm{GUT}})", "\\lambda (M_{\\mathrm{GUT}})",
         "\\tilde{f}_{11}(M_{\\mathrm{GUT}})", "\\tilde{f}_{12}(M_{\\mathrm{GUT}})",
         "\\tilde{f}_{21}(M_{\\mathrm{GUT}})", "\\tilde{f}_{22}(M_{\\mathrm{GUT}})",
         "\\tilde{f}_{31}(M_{\\mathrm{GUT}})", "\\tilde{f}_{32}(M_{\\mathrm{GUT}})",
         "f_{11}(M_{\\mathrm{GUT}})", "f_{12}(M_{\\mathrm{GUT}})", "f_{21}(M_{\\mathrm{GUT}})",
         "f_{22}(M_{\\mathrm{GUT}})", "f_{31}(M_{\\mathrm{GUT}})", "f_{32}(M_{\\mathrm{GUT}})",
         "\\mu_L(M_{\\mathrm{GUT}})", "\\mu(M_{\\mathrm{GUT}})", "B_{\\mu_L}(M_{\\mathrm{GUT}})",
         "B_\\mu(M_{\\mathrm{GUT}})"};

   } // namespace semianalytic

   const unsigned mixing_dimensions[NUMBER_OF_MIXINGS] = {6, 3, 6, 6
      , 6, 5, 5, 2, 8, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 7, 4, 
      2, 2, 7, 2, 2, 2};

   const char* mixing_names[NUMBER_OF_MIXINGS] = {"ZD", "ZV", "ZU" 
      , "ZE", "ZDX", "ZH", "ZA", "ZP", "ZN", "UM", "UP", "ZEL"
      , "ZER", "ZDL", "ZDR", "ZUL", "ZUR", "ZDXL", "ZDXR", "UHI0"
      , "UHIPM", "ZMI", "ZPI", "ZNI", "UHp0", "UHpp", "ZNp"};

   const char* mixing_latex_names[NUMBER_OF_MIXINGS] = {"Z_d", "Z_{\\nu}", "Z_u"
      , "Z_e", "Z_D", "U_h", "U_A", "U_{H^\\pm}", "N", "U", "V"
      , "U_{e_L}", "U_{e_R}", "U_{d_L}", "U_{d_R}", "U_{u_L}"
      , "U_{u_R}", "U_{D_L}", "U_{D_R}", "U_{h_I}", "U_{H^\\pm_I}"
      , "U_I", "V_I", "N_I", "U_{L_4^0}", "U_{L_4^\\pm}", "N_{L_4}"};

   const char* model_name = "CSE6SSM";
   const bool is_low_energy_model = false;
   const bool is_supersymmetric_model = true;

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                " << model_name << '\n'
      << "Is a low-energy model:     "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model: "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Number of multiplets:      " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:      " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                ";
   for (unsigned i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                ";
   for (unsigned i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }
   ostr << '\n';
}

} // namespace CSE6SSM_info

} // namespace flexiblesusy

