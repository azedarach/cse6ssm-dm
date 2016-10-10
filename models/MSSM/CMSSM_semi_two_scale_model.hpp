
/**
 * @file CMSSM_semi_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the semianalytic two-scale solver by
 *        solving EWSB and determine the pole masses and mixings
 */

#ifndef CMSSM_SEMI_TWO_SCALE_MODEL_H
#define CMSSM_SEMI_TWO_SCALE_MODEL_H

#include "CMSSM_semi_model.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "MSSM_info.hpp"
#include "CMSSM_semi_two_scale_input_parameters.hpp"
#include "two_scale_model.hpp"
#include "config.h"

#include <iosfwd>
#include <string>

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
class Two_scale;
/**
 * @class CMSSM_semianalytic<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class CMSSM_semianalytic<Two_scale> : public Two_scale_model, public MSSM_mass_eigenstates {
public:
   explicit CMSSM_semianalytic(const CMSSM_semianalytic_input_parameters<Two_scale>& input_ = CMSSM_semianalytic_input_parameters<Two_scale>());
   virtual ~CMSSM_semianalytic();

   /// number of ewsb equations
   static const std::size_t number_of_tadpole_equations = 2;

   virtual void clear();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(unsigned);
   void set_number_of_ewsb_iterations(std::size_t);
   const CMSSM_semianalytic_input_parameters<Two_scale>& get_input() const;
   CMSSM_semianalytic_input_parameters<Two_scale>& get_input();
   void set_input_parameters(const CMSSM_semianalytic_input_parameters<Two_scale>&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level
   double get_ewsb_output_parameter(unsigned) const;

   /// calculate the tadpoles at current loop order
   void tadpole_equations(double[number_of_tadpole_equations]) const;

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);

   // model interface
   double get_parameter(unsigned) const;
   void set_parameter(unsigned, double);

   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TYd() const { return TYd_Azero_coeff; }
   double get_Azero_coeff_TYd(int i, int k) const { return TYd_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TYd() const { return TYd_m12_coeff; }
   double get_m12_coeff_TYd(int i, int k) const { return TYd_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TYe() const { return TYe_Azero_coeff; }
   double get_Azero_coeff_TYe(int i, int k) const { return TYe_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TYe() const { return TYe_m12_coeff; }
   double get_m12_coeff_TYe(int i, int k) const { return TYe_m12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero_coeff_TYu() const { return TYu_Azero_coeff; }
   double get_Azero_coeff_TYu(int i, int k) const { return TYu_Azero_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m12_coeff_TYu() const { return TYu_m12_coeff; }
   double get_m12_coeff_TYu(int i, int k) const { return TYu_m12_coeff(i,k); }
   double get_BMu0_coeff_BMu() const { return BMu_BMu0_coeff; }
   double get_Azero_coeff_BMu() const { return BMu_Azero_coeff; }
   double get_m12_coeff_BMu() const { return BMu_m12_coeff; }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mq2() const { return mq2_m02_coeff; }
   double get_m02_coeff_mq2(int i, int k) const { return mq2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mq2() const { return mq2_m122_coeff; }
   double get_m122_coeff_mq2(int i, int k) const { return mq2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mq2() const { return mq2_Azerom12_coeff; }
   double get_Azerom12_coeff_mq2(int i, int k) const { return mq2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mq2() const { return mq2_Azero2_coeff; }
   double get_Azero2_coeff_mq2(int i, int k) const { return mq2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_ml2() const { return ml2_m02_coeff; }
   double get_m02_coeff_ml2(int i, int k) const { return ml2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_ml2() const { return ml2_m122_coeff; }
   double get_m122_coeff_ml2(int i, int k) const { return ml2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_ml2() const { return ml2_Azerom12_coeff; }
   double get_Azerom12_coeff_ml2(int i, int k) const { return ml2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_ml2() const { return ml2_Azero2_coeff; }
   double get_Azero2_coeff_ml2(int i, int k) const { return ml2_Azero2_coeff(i,k); }
   double get_m02_coeff_mHd2() const { return mHd2_m02_coeff; }
   double get_m122_coeff_mHd2() const { return mHd2_m122_coeff; }
   double get_Azerom12_coeff_mHd2() const { return mHd2_Azerom12_coeff; }
   double get_Azero2_coeff_mHd2() const { return mHd2_Azero2_coeff; }
   double get_m02_coeff_mHu2() const { return mHu2_m02_coeff; }
   double get_m122_coeff_mHu2() const { return mHu2_m122_coeff; }
   double get_Azerom12_coeff_mHu2() const { return mHu2_Azerom12_coeff; }
   double get_Azero2_coeff_mHu2() const { return mHu2_Azero2_coeff; }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_mu2() const { return mu2_m02_coeff; }
   double get_m02_coeff_mu2(int i, int k) const { return mu2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_mu2() const { return mu2_m122_coeff; }
   double get_m122_coeff_mu2(int i, int k) const { return mu2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_mu2() const { return mu2_Azerom12_coeff; }
   double get_Azerom12_coeff_mu2(int i, int k) const { return mu2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_mu2() const { return mu2_Azero2_coeff; }
   double get_Azero2_coeff_mu2(int i, int k) const { return mu2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_md2() const { return md2_m02_coeff; }
   double get_m02_coeff_md2(int i, int k) const { return md2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_md2() const { return md2_m122_coeff; }
   double get_m122_coeff_md2(int i, int k) const { return md2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_md2() const { return md2_Azerom12_coeff; }
   double get_Azerom12_coeff_md2(int i, int k) const { return md2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_md2() const { return md2_Azero2_coeff; }
   double get_Azero2_coeff_md2(int i, int k) const { return md2_Azero2_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m02_coeff_me2() const { return me2_m02_coeff; }
   double get_m02_coeff_me2(int i, int k) const { return me2_m02_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_m122_coeff_me2() const { return me2_m122_coeff; }
   double get_m122_coeff_me2(int i, int k) const { return me2_m122_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azerom12_coeff_me2() const { return me2_Azerom12_coeff; }
   double get_Azerom12_coeff_me2(int i, int k) const { return me2_Azerom12_coeff(i,k); }
   const Eigen::Matrix<double,3,3>& get_Azero2_coeff_me2() const { return me2_Azero2_coeff; }
   double get_Azero2_coeff_me2(int i, int k) const { return me2_Azero2_coeff(i,k); }
   double get_Azero_coeff_MassB() const { return MassB_Azero_coeff; }
   double get_m12_coeff_MassB() const { return MassB_m12_coeff; }
   double get_Azero_coeff_MassWB() const { return MassWB_Azero_coeff; }
   double get_m12_coeff_MassWB() const { return MassWB_m12_coeff; }
   double get_Azero_coeff_MassG() const { return MassG_Azero_coeff; }
   double get_m12_coeff_MassG() const { return MassG_m12_coeff; }

   void calculate_coefficients(double);

private:
   struct EWSB_args {
      CMSSM_semianalytic<Two_scale>* model;
      unsigned ewsb_loop_order;
   };

   // input parameters
   CMSSM_semianalytic_input_parameters<Two_scale> input;

   std::size_t number_of_ewsb_iterations;
   unsigned ewsb_loop_order;
   double precision;              ///< RG running precision
   double ewsb_iteration_precision;
   static const std::size_t number_of_fit_points = 4;

   // solution for ewsb_loop_order > 1
   static const std::size_t number_of_ewsb_equations = 2;

   Eigen::Array<double,number_of_tadpole_equations,1> ewsb_solution;

   void ewsb_equations(double[number_of_ewsb_equations]) const;
   static int ewsb_equations(const gsl_vector*, void*, gsl_vector*);

   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(unsigned);
   int solve_ewsb_iteratively_with(EWSB_solver*, const double[number_of_ewsb_equations]);
   int check_ewsb_solution(double);
   int ewsb_initial_guess(double[number_of_ewsb_equations]);
   int ewsb_step(double[number_of_ewsb_equations]);
   static int ewsb_step(const gsl_vector*, void*, gsl_vector*);
   static int tadpole_equations(const gsl_vector*, void*, gsl_vector*);

   void set_soft_parameters_at_input_scale(double,double,double,double);
   void set_soft_parameters_at_current_scale(double,double,double,double);

   Eigen::Matrix<double,3,3> TYd_Azero_coeff;
   Eigen::Matrix<double,3,3> TYd_m12_coeff;
   Eigen::Matrix<double,3,3> TYe_Azero_coeff;
   Eigen::Matrix<double,3,3> TYe_m12_coeff;
   Eigen::Matrix<double,3,3> TYu_Azero_coeff;
   Eigen::Matrix<double,3,3> TYu_m12_coeff;
   double BMu_BMu0_coeff;
   double BMu_Azero_coeff;
   double BMu_m12_coeff;
   Eigen::Matrix<double,3,3> mq2_m02_coeff;
   Eigen::Matrix<double,3,3> mq2_m122_coeff;
   Eigen::Matrix<double,3,3> mq2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mq2_Azero2_coeff;
   Eigen::Matrix<double,3,3> ml2_m02_coeff;
   Eigen::Matrix<double,3,3> ml2_m122_coeff;
   Eigen::Matrix<double,3,3> ml2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> ml2_Azero2_coeff;
   double mHd2_m02_coeff;
   double mHd2_m122_coeff;
   double mHd2_Azerom12_coeff;
   double mHd2_Azero2_coeff;
   double mHu2_m02_coeff;
   double mHu2_m122_coeff;
   double mHu2_Azerom12_coeff;
   double mHu2_Azero2_coeff;
   Eigen::Matrix<double,3,3> md2_m02_coeff;
   Eigen::Matrix<double,3,3> md2_m122_coeff;
   Eigen::Matrix<double,3,3> md2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> md2_Azero2_coeff;
   Eigen::Matrix<double,3,3> mu2_m02_coeff;
   Eigen::Matrix<double,3,3> mu2_m122_coeff;
   Eigen::Matrix<double,3,3> mu2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> mu2_Azero2_coeff;
   Eigen::Matrix<double,3,3> me2_m02_coeff;
   Eigen::Matrix<double,3,3> me2_m122_coeff;
   Eigen::Matrix<double,3,3> me2_Azerom12_coeff;
   Eigen::Matrix<double,3,3> me2_Azero2_coeff;
   double MassB_Azero_coeff;
   double MassB_m12_coeff;
   double MassWB_Azero_coeff;
   double MassWB_m12_coeff;
   double MassG_Azero_coeff;
   double MassG_m12_coeff;

   bool has_previous_ewsb_solution;
};

std::ostream& operator<<(std::ostream&, const CMSSM_semianalytic<Two_scale>&);

} // namespace flexiblesusy

#endif
