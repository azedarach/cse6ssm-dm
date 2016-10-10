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

#ifndef CSE6SSM_INFO_H
#define CSE6SSM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace CSE6SSM_info {
   enum Particles : unsigned {VG, Glu, Fv, ChaP, VP, VZ, VZp, Sd, Sv, Su, Se,
      SDX, hh, Ah, Hpm, Chi, Cha, Fe, Fd, Fu, FDX, SHI0, SHIPM, ChaI, ChiI, SHp0,
      SHpp, ChiP, VWm, NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, hE00, hE01, hE10, hE11, hE20, hE21, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12
      , Ye20, Ye21, Ye22, SigmaL, KappaPr, Sigmax, gD00, gD01, gD02, gD10, gD11,
      gD12, gD20, gD21, gD22, Kappa00, Kappa01, Kappa02, Kappa10, Kappa11, Kappa12
      , Kappa20, Kappa21, Kappa22, Lambda1200, Lambda1201, Lambda1210, Lambda1211,
      Lambdax, fu00, fu01, fu10, fu11, fu20, fu21, fd00, fd01, fd10, fd11, fd20,
      fd21, Yu00, Yu01, Yu02, Yu10, Yu11, Yu12, Yu20, Yu21, Yu22, MuPr, MuPhi, XiF
      , g1, g2, g3, g1p, vd, vu, vs, vsb, vphi, QS, TYd00, TYd01, TYd02, TYd10, TYd11,
      TYd12, TYd20, TYd21, TYd22, ThE00, ThE01, ThE10, ThE11, ThE20, ThE21, TYe00
      , TYe01, TYe02, TYe10, TYe11, TYe12, TYe20, TYe21, TYe22, TSigmaL, TKappaPr,
      TSigmax, TgD00, TgD01, TgD02, TgD10, TgD11, TgD12, TgD20, TgD21, TgD22,
      TKappa00, TKappa01, TKappa02, TKappa10, TKappa11, TKappa12, TKappa20,
      TKappa21, TKappa22, TLambda1200, TLambda1201, TLambda1210, TLambda1211,
      TLambdax, Tfu00, Tfu01, Tfu10, Tfu11, Tfu20, Tfu21, Tfd00, Tfd01, Tfd10,
      Tfd11, Tfd20, Tfd21, TYu00, TYu01, TYu02, TYu10, TYu11, TYu12, TYu20, TYu21,
      TYu22, BMuPr, BMuPhi, LXiF, mq200, mq201, mq202, mq210, mq211, mq212, mq220
      , mq221, mq222, ml200, ml201, ml202, ml210, ml211, ml212, ml220, ml221,
      ml222, mHd2, mHu2, md200, md201, md202, md210, md211, md212, md220, md221,
      md222, mu200, mu201, mu202, mu210, mu211, mu212, mu220, mu221, mu222, me200,
      me201, me202, me210, me211, me212, me220, me221, me222, ms2, msbar2,
      mH1I200, mH1I201, mH1I210, mH1I211, mH2I200, mH2I201, mH2I210, mH2I211,
      mSI200, mSI201, mSI202, mSI210, mSI211, mSI212, mSI220, mSI221, mSI222,
      mDx200, mDx201, mDx202, mDx210, mDx211, mDx212, mDx220, mDx221, mDx222,
      mDxbar200, mDxbar201, mDxbar202, mDxbar210, mDxbar211, mDxbar212, mDxbar220,
      mDxbar221, mDxbar222, mHp2, mHpbar2, mphi2, MassB, MassWB, MassG, MassBp,
      NUMBER_OF_PARAMETERS};

   // temp solution to allow for both sets of inputs without changing names
   namespace two_scale {

      // DH:: added enum for input parameters as well
      enum Inputs : unsigned {m0, m12, TanBeta, SignLambdax, Azero, ssumInput,
            QSInput, hEInput00, hEInput01, hEInput10, hEInput11, hEInput20,
            hEInput21, SigmaLInput, KappaPrInput, SigmaxInput, gDInput00,
            gDInput01, gDInput02, gDInput10, gDInput11, gDInput12, gDInput20,
            gDInput21, gDInput22, KappaInput00, KappaInput01, KappaInput02,
            KappaInput10, KappaInput11, KappaInput12, KappaInput20, KappaInput21,
            KappaInput22, Lambda12Input00, Lambda12Input01, Lambda12Input10,
            Lambda12Input11, fuInput00, fuInput01, fuInput10, fuInput11,
            fuInput20, fuInput21, fdInput00, fdInput01, fdInput10, fdInput11,
            fdInput20, fdInput21, MuPrInput, MuPhiInput, BMuPrInput, BMuPhiInput,
            NUMBER_OF_INPUTS};

      extern const char* input_names[NUMBER_OF_INPUTS];
      extern const int input_mass_dimensions[NUMBER_OF_INPUTS];
      extern const char* input_latex_names[NUMBER_OF_INPUTS];

   } // namespace two_scale

   namespace semianalytic {

      enum Inputs : unsigned {m12, Azero, TanBeta, sInput,
            QSInput, hEInput00, hEInput01, hEInput10, hEInput11, hEInput20,
            hEInput21, SigmaLInput, KappaPrInput, SigmaxInput, gDInput00,
            gDInput01, gDInput02, gDInput10, gDInput11, gDInput12, gDInput20,
            gDInput21, gDInput22, KappaInput00, KappaInput01, KappaInput02,
            KappaInput10, KappaInput11, KappaInput12, KappaInput20, KappaInput21,
            KappaInput22, Lambda12Input00, Lambda12Input01, Lambda12Input10,
            Lambda12Input11, LambdaxInput, fuInput00, fuInput01, fuInput10, fuInput11,
            fuInput20, fuInput21, fdInput00, fdInput01, fdInput10, fdInput11,
            fdInput20, fdInput21, MuPrInput, MuPhiInput, BMuPrInput, BMuPhiInput,
            NUMBER_OF_INPUTS};

      extern const char* input_names[NUMBER_OF_INPUTS];
      extern const int input_mass_dimensions[NUMBER_OF_INPUTS];
      extern const char* input_latex_names[NUMBER_OF_INPUTS];

   } // namespace semianalytic

   // DH:: and enum for mixings
   enum Mixings : unsigned {ZD, ZV, ZU, ZE, ZDX, ZH, ZA, ZP, ZN, UM, UP, ZEL,
         ZER, ZDL, ZDR, ZUL, ZUR, ZDXL, ZDXR, UHI0, UHIPM, ZMI, ZPI, ZNI, 
         UHp0, UHpp, ZNp, NUMBER_OF_MIXINGS};

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;
   extern const double normalization_g1p;

   extern const unsigned particle_multiplicities[NUMBER_OF_PARTICLES];
   extern const char* particle_names[NUMBER_OF_PARTICLES];
   extern const char* particle_latex_names[NUMBER_OF_PARTICLES];
   extern const char* parameter_names[NUMBER_OF_PARAMETERS];
   extern const int parameter_mass_dimensions[NUMBER_OF_PARAMETERS];
   extern const char* parameter_latex_names[NUMBER_OF_PARAMETERS];

   extern const unsigned mixing_dimensions[NUMBER_OF_MIXINGS];
   extern const char* mixing_names[NUMBER_OF_MIXINGS];
   extern const char* mixing_latex_names[NUMBER_OF_MIXINGS];

   extern const char* model_name;
   extern const bool is_low_energy_model;
   extern const bool is_supersymmetric_model;

   void print(std::ostream&);
}

} // namespace flexiblesusy

#endif
