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

/**
 * @file sfermions.cpp
 * @brief finding mass eigenstates and mixing of sfermions in absence of 
 *        family mixing, where we have a 2 by 2 mass matrix.
 */

#include "sfermions.hpp"
#include "linalg2.hpp"
#include "wrappers.hpp"
#include "logger.hpp"

namespace flexiblesusy {
namespace sfermions {

static const double oneOverRoot2 = 1./sqrt(2.); // 0.7071067811865475

const double Isospin[NUMBER_OF_MSSM_SPARTICLES] = {
   0.5, -0.5, 0.5, -0.5
};

const double Hypercharge_left[NUMBER_OF_MSSM_SPARTICLES] = {
   1./3., 1./3., -1., -1.
};

const double Hypercharge_right[NUMBER_OF_MSSM_SPARTICLES] = {
   -4./3., 2./3., 0., 2.
};

const double U1prime_charge_left[NUMBER_OF_MSSM_SPARTICLES] = {
   1., 1., 2., 2.
};

const double U1prime_charge_right[NUMBER_OF_MSSM_SPARTICLES] = {
   1., 2., 0., 1.
};

/**
 * Obtains 2 x 2 mass matrix using input parameters in first argument 
 * and diagonalises it.  Fills the second argument with the eigenvalues
 * and returns the mixing angle.
 */ 
double diagonalize_sfermions_2x2(const Mass_data& pars,
                                 Eigen::Array<double,2,1>& msf)
{
   const double ml2    = pars.ml2;
   const double mr2    = pars.mr2;
   const double yf     = pars.yf;
   const double vd     = pars.vd;
   const double vu     = pars.vu;
   const double gY     = pars.gY;
   const double g2     = pars.g2;
   const double Tyf    = pars.Tyf;
   const double mu     = pars.mu;
   const double T3     = pars.T3;
   const double Yl     = pars.Yl;
   const double Yr     = pars.Yr;
   const double vev2   = 0.25 * (Sqr(vd) - Sqr(vu));
   Eigen::Matrix<double,2,2> mass_matrix;
   /// fill sfermion phi in mass matix in basis (phi_L phi_R)
   if (Sign(T3) > 0) {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2;
      mass_matrix(0,1) = oneOverRoot2 * (vu*Conj(Tyf) - vd*Conj(yf)*mu);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         - 0.5 * Yr * Sqr(gY) * vev2;
   } else {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2;
      mass_matrix(0,1) = oneOverRoot2 * (vd*Conj(Tyf) - vu*Conj(yf)*mu);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         - 0.5 * Yr * Sqr(gY) * vev2;
   }

   Eigen::Matrix<double, 2, 2> Zf;
   diagonalize_hermitian(mass_matrix, msf, Zf);

#ifdef ENABLE_VERBOSE
   if (msf.minCoeff() < 0.)
      WARNING("diagonalize_sfermions_2x2: sfermion tachyon");
#endif

   msf = AbsSqrt(msf);

   double theta;

   if (Sign(Zf(0,0)) == Sign(Zf(1,1))) {
      theta = ArcCos(Abs(Zf(0,0)));
   } else {
      theta = ArcCos(Abs(Zf(0,1)));
      Zf.col(0).swap(Zf.col(1));
      std::swap(msf(0), msf(1));
   }

   theta = Sign(mass_matrix(0,1) / (mass_matrix(0,0) - mass_matrix(1,1)))
      * Abs(theta);

   return theta;
}

/**
 * Obtains 2 x 2 mass matrix using input parameters in first argument 
 * and diagonalises it, for a U(1)-extended model with one SM singlet.  
 * Fills the second argument with the eigenvalues and returns the mixing angle.
 */ 
double diagonalize_sfermions_2x2(const U1_extended_mass_data_one_singlet& pars,
                                 Eigen::Array<double,2,1>& msf)
{
   const double ml2    = pars.ml2;
   const double mr2    = pars.mr2;
   const double yf     = pars.yf;
   const double vd     = pars.vd;
   const double vu     = pars.vu;
   const double vs     = pars.vs;
   const double QHd    = pars.QHd;
   const double QHu    = pars.QHu;
   const double QS     = pars.QS;
   const double gY     = pars.gY;
   const double g2     = pars.g2;
   const double gN     = pars.gN;
   const double Tyf    = pars.Tyf;
   const double mueff  = pars.mueff;
   const double T3     = pars.T3;
   const double Yl     = pars.Yl;
   const double Yr     = pars.Yr;
   const double Ql     = pars.Ql;
   const double Qr     = pars.Qr;
   const double vev2   = 0.25 * (Sqr(vd) - Sqr(vu));
   const double deltap = QHd * Sqr(vd) + QHu * Sqr(vu) + QS * Sqr(vs);
   Eigen::Matrix<double,2,2> mass_matrix;
   /// fill sfermion phi in mass matix in basis (phi_L phi_R)
   if (Sign(T3) > 0) {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2
         + 0.5 * Sqr(gN) * Ql * deltap;
      mass_matrix(0,1) = oneOverRoot2 * (vu*Conj(Tyf) - vd*Conj(yf)*mueff);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         - 0.5 * Yr * Sqr(gY) * vev2 + 0.5 * Sqr(gN) * Qr * deltap;
   } else {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2
         + 0.5 * Sqr(gN) * Ql * deltap;
      mass_matrix(0,1) = oneOverRoot2 * (vd*Conj(Tyf) - vu*Conj(yf)*mueff);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         - 0.5 * Yr * Sqr(gY) * vev2 + 0.5 * Sqr(gN) * Qr * deltap;
   }

   Eigen::Matrix<double, 2, 2> Zf;
   diagonalize_hermitian(mass_matrix, msf, Zf);

#ifdef ENABLE_VERBOSE
   if (msf.minCoeff() < 0.)
      WARNING("diagonalize_sfermions_2x2: sfermion tachyon");
#endif

   msf = AbsSqrt(msf);

   double theta;

   if (Sign(Zf(0,0)) == Sign(Zf(1,1))) {
      theta = ArcCos(Abs(Zf(0,0)));
   } else {
      theta = ArcCos(Abs(Zf(0,1)));
      Zf.col(0).swap(Zf.col(1));
      std::swap(msf(0), msf(1));
   }

   theta = Sign(mass_matrix(0,1) / (mass_matrix(0,0) - mass_matrix(1,1)))
      * Abs(theta);

   return theta;
}

/**
 * Obtains 2 x 2 mass matrix using input parameters in first argument 
 * and diagonalises it, for a U(1)-extended model with two SM singlets.  
 * Fills the second argument with the eigenvalues and returns the mixing angle.
 */ 
double diagonalize_sfermions_2x2(const U1_extended_mass_data_two_singlets& pars,
                                 Eigen::Array<double,2,1>& msf)
{
   const double ml2    = pars.ml2;
   const double mr2    = pars.mr2;
   const double yf     = pars.yf;
   const double vd     = pars.vd;
   const double vu     = pars.vu;
   const double vs     = pars.vs;
   const double vsb    = pars.vsb;
   const double QHd    = pars.QHd;
   const double QHu    = pars.QHu;
   const double QS     = pars.QS;
   const double QSb    = pars.QSb;
   const double gY     = pars.gY;
   const double g2     = pars.g2;
   const double gN     = pars.gN;
   const double Tyf    = pars.Tyf;
   const double mueff  = pars.mueff;
   const double T3     = pars.T3;
   const double Yl     = pars.Yl;
   const double Yr     = pars.Yr;
   const double Ql     = pars.Ql;
   const double Qr     = pars.Qr;
   const double vev2   = 0.25 * (Sqr(vd) - Sqr(vu));
   const double deltap = QHd * Sqr(vd) + QHu * Sqr(vu) + QS * Sqr(vs) + QSb * Sqr(vsb);
   Eigen::Matrix<double,2,2> mass_matrix;
   /// fill sfermion phi in mass matix in basis (phi_L phi_R)
   if (Sign(T3) > 0) {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2
         + 0.5 * Sqr(gN) * Ql * deltap;
      mass_matrix(0,1) = oneOverRoot2 * (vu*Conj(Tyf) - vd*Conj(yf)*mueff);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         - 0.5 * Yr * Sqr(gY) * vev2 + 0.5 * Sqr(gN) * Qr * deltap;
   } else {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2
         + 0.5 * Sqr(gN) * Ql * deltap;
      mass_matrix(0,1) = oneOverRoot2 * (vd*Conj(Tyf) - vu*Conj(yf)*mueff);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         - 0.5 * Yr * Sqr(gY) * vev2 + 0.5 * Sqr(gN) * Qr * deltap;
   }

   Eigen::Matrix<double, 2, 2> Zf;
   diagonalize_hermitian(mass_matrix, msf, Zf);

#ifdef ENABLE_VERBOSE
   if (msf.minCoeff() < 0.)
      WARNING("diagonalize_sfermions_2x2: sfermion tachyon");
#endif

   msf = AbsSqrt(msf);

   double theta;

   if (Sign(Zf(0,0)) == Sign(Zf(1,1))) {
      theta = ArcCos(Abs(Zf(0,0)));
   } else {
      theta = ArcCos(Abs(Zf(0,1)));
      Zf.col(0).swap(Zf.col(1));
      std::swap(msf(0), msf(1));
   }

   theta = Sign(mass_matrix(0,1) / (mass_matrix(0,0) - mass_matrix(1,1)))
      * Abs(theta);

   return theta;
}

} // namespace sfermions
} // namespace flexiblesusy
