// ====================================================================
// Routines used in scan class
// ====================================================================

#ifndef CSE6SSM_SCAN_UTILITIES_H
#define CSE6SSM_SCAN_UTILITIES_H

#include "CSE6SSM_two_scale_model_slha.hpp"
#include "CSE6SSM_semi_two_scale_model_slha.hpp"
#include "CSE6SSM_info.hpp"
#include "CSE6SSM_utilities.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <cassert>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

#define PHYSICAL(p) model.get_physical().p
#define PHYSICALSLHA(p) model.get_physical_slha().p
#define DRBARSLHA(p) model.get_drbar_slha().p
#define MODELPARAMETER(p) model.get_##p()
#define COEFFICIENT(p,q) model.get_##q##_coeff_##p()


namespace flexiblesusy {

void write_CSE6SSM_inputs(const CSE6SSM_input_parameters<Two_scale>&,
                          std::ostream &, std::size_t);
void write_CSE6SSM_inputs_list(std::ostream &, std::size_t);
void write_CSE6SSM_inputs(
   const CSE6SSM_semianalytic_input_parameters<Two_scale>&,
   std::ostream &, std::size_t
   );
void write_CSE6SSM_semianalytic_inputs_list(std::ostream &, std::size_t);
std::valarray<double> to_valarray(double);

template <class Scalar, int M, int N>
std::valarray<double> to_valarray(const Eigen::Array<Scalar, M, N>& v)
{
   return std::valarray<double>(v.data(), v.size());
}

template <int M, int N>
std::valarray<double> to_valarray(const Eigen::Matrix<double, M, N>& v)
{
   return std::valarray<double>(v.data(), v.size());
}

// this needs to extract real and imaginary parts, then
// interleave them in the valarray
template <int M, int N>
std::valarray<double> to_valarray(
   const Eigen::Matrix<std::complex<double>, M, N>& v
   )
{
   std::size_t num_entries = v.size();
   Eigen::Matrix<double,M,N> re_v = v.real();
   std::valarray<double> real_parts(re_v.data(), num_entries);
   Eigen::Matrix<double,M,N> im_v = v.imag();
   std::valarray<double>imag_parts(im_v.data(), num_entries);

   // combine into a single valarray
   std::valarray<double> result(2 * v.size());
   result[std::slice(0, num_entries, 1)] = real_parts;
   result[std::slice(num_entries, num_entries, 1)] = imag_parts;
   return result;
}

class CSE6SSM_pole_mass_writer {
public:
   CSE6SSM_pole_mass_writer();
   ~CSE6SSM_pole_mass_writer() {}

   template <class T>
   void extract_pole_masses(const CSE6SSM<T>&);
   void write_pole_masses_comment_line(std::ostream &) const;
   void write_pole_masses_line(std::ostream &) const;

private:
   struct TPoleMass {
      std::string name;
      std::valarray<double> masses;
      TPoleMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };
   typedef std::vector<TPoleMass> TPoleMasses;
   TPoleMasses pole_masses;

   CSE6SSM_input_parameters<Two_scale> pole_masses_inputs;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> pole_masses_problems;

   double pole_masses_scale;
   unsigned width;
};

class CSE6SSM_semianalytic_pole_mass_writer {
public:
   CSE6SSM_semianalytic_pole_mass_writer();
   ~CSE6SSM_semianalytic_pole_mass_writer() {}

   template <class T>
   void extract_pole_masses(const CSE6SSM_semianalytic<T>&);
   void write_pole_masses_comment_line(std::ostream &) const;
   void write_pole_masses_line(std::ostream &) const;

private:
   struct TPoleMass {
      std::string name;
      std::valarray<double> masses;
      TPoleMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };
   typedef std::vector<TPoleMass> TPoleMasses;
   TPoleMasses pole_masses;

   CSE6SSM_semianalytic_input_parameters<Two_scale> pole_masses_inputs;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> pole_masses_problems;

   double pole_masses_scale;
   double m0Sqr;
   unsigned width;
};

class CSE6SSM_drbar_values_writer {
public:
   CSE6SSM_drbar_values_writer();
   ~CSE6SSM_drbar_values_writer() {}

   template <class T>
   void extract_drbar_masses(const CSE6SSM<T>&);
   template <class T>
   void extract_drbar_susy_pars(const CSE6SSM<T>&);
   template <class T>
   void extract_drbar_soft_pars(const CSE6SSM<T>&);
   template <class T>
   void extract_drbar_mixings(const CSE6SSM<T>&);

   void write_drbar_masses_comment_line(std::ostream &) const;
   void write_drbar_susy_pars_comment_line(std::ostream &) const;
   void write_drbar_soft_pars_comment_line(std::ostream &) const;
   void write_drbar_mixings_comment_line(std::ostream &) const;

   void write_drbar_masses_line(std::ostream &) const;
   void write_drbar_susy_pars_line(std::ostream &) const;
   void write_drbar_soft_pars_line(std::ostream &) const;
   void write_drbar_mixings_line(std::ostream &) const;

private:
   struct TMass {
      std::string name;
      std::valarray<double> masses;
      TMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TParameter {
      std::string name;
      std::size_t mass_dimension;
      std::size_t rows;
      std::size_t cols;
      std::valarray<double> values;
      TParameter(const std::string& name_,
                 std::size_t mass_dimension_,
                 std::size_t rows_, std::size_t cols_,
                 const std::valarray<double>& values_)
         : name(name_)
         , mass_dimension(mass_dimension_)
         , rows(rows_)
         , cols(cols_)
         , values(values_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TMixing {
      std::string name;
      std::size_t dimension;
      bool is_real;
      std::valarray<double> mixings;
      TMixing(const std::string& name_, std::size_t dimension_,
              bool is_real_,
              const std::valarray<double>& mixings_)
         : name(name_)
         , dimension(dimension_)
         , is_real(is_real_)
         , mixings(mixings_)
         {}
   };

   typedef std::vector<TMass> TRunningMasses;
   typedef std::vector<TParameter> TSusyPars;
   typedef std::vector<TParameter> TSoftPars;
   typedef std::vector<TMixing> TMixings;

   TRunningMasses drbar_masses;
   TSusyPars drbar_susy_pars;
   TSoftPars drbar_soft_pars;
   TMixings drbar_mixings;

   CSE6SSM_input_parameters<Two_scale> drbar_masses_inputs;
   CSE6SSM_input_parameters<Two_scale> drbar_susy_pars_inputs;
   CSE6SSM_input_parameters<Two_scale> drbar_soft_pars_inputs;
   CSE6SSM_input_parameters<Two_scale> drbar_mixings_inputs;

   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_masses_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_susy_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_soft_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_mixings_problems;

   double drbar_masses_scale;
   double drbar_susy_pars_scale;
   double drbar_soft_pars_scale;
   double drbar_mixings_scale;
   unsigned width;
};

class CSE6SSM_semianalytic_drbar_values_writer {
public:
   CSE6SSM_semianalytic_drbar_values_writer();
   ~CSE6SSM_semianalytic_drbar_values_writer() {}

   template <class T>
   void extract_drbar_masses(const CSE6SSM_semianalytic<T>&);
   template <class T>
   void extract_drbar_susy_pars(const CSE6SSM_semianalytic<T>&);
   template <class T>
   void extract_drbar_soft_pars(const CSE6SSM_semianalytic<T>&);
   template <class T>
   void extract_drbar_mixings(const CSE6SSM_semianalytic<T>&);

   void write_drbar_masses_comment_line(std::ostream &) const;
   void write_drbar_susy_pars_comment_line(std::ostream &) const;
   void write_drbar_soft_pars_comment_line(std::ostream &) const;
   void write_drbar_mixings_comment_line(std::ostream &) const;

   void write_drbar_masses_line(std::ostream &) const;
   void write_drbar_susy_pars_line(std::ostream &) const;
   void write_drbar_soft_pars_line(std::ostream &) const;
   void write_drbar_mixings_line(std::ostream &) const;

private:
   struct TMass {
      std::string name;
      std::valarray<double> masses;
      TMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TParameter {
      std::string name;
      std::size_t mass_dimension;
      std::size_t rows;
      std::size_t cols;
      std::valarray<double> values;
      TParameter(const std::string& name_,
                 std::size_t mass_dimension_,
                 std::size_t rows_, std::size_t cols_,
                 const std::valarray<double>& values_)
         : name(name_)
         , mass_dimension(mass_dimension_)
         , rows(rows_)
         , cols(cols_)
         , values(values_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TMixing {
      std::string name;
      std::size_t dimension;
      bool is_real;
      std::valarray<double> mixings;
      TMixing(const std::string& name_, std::size_t dimension_,
              bool is_real_,
              const std::valarray<double>& mixings_)
         : name(name_)
         , dimension(dimension_)
         , is_real(is_real_)
         , mixings(mixings_)
         {}
   };

   typedef std::vector<TMass> TRunningMasses;
   typedef std::vector<TParameter> TSusyPars;
   typedef std::vector<TParameter> TSoftPars;
   typedef std::vector<TMixing> TMixings;

   TRunningMasses drbar_masses;
   TSusyPars drbar_susy_pars;
   TSoftPars drbar_soft_pars;
   TMixings drbar_mixings;

   CSE6SSM_semianalytic_input_parameters<Two_scale> drbar_masses_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> drbar_susy_pars_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> drbar_soft_pars_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> drbar_mixings_inputs;

   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_masses_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_susy_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_soft_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> drbar_mixings_problems;

   double drbar_masses_scale;
   double drbar_susy_pars_scale;
   double drbar_soft_pars_scale;
   double drbar_mixings_scale;
   unsigned width;
};

/// As above, but in SLHA convention
class CSE6SSM_slha_values_writer {
public:
   CSE6SSM_slha_values_writer();
   ~CSE6SSM_slha_values_writer() {}

   void set_high_scale(double high_scale_) { high_scale = high_scale_; }
   void set_susy_scale(double susy_scale_) { susy_scale = susy_scale_; }
   void set_low_scale(double low_scale_) { low_scale = low_scale_; }

   template <class T>
   void extract_slha_pole_masses(const CSE6SSM_slha<T>&);
   template <class T>
   void extract_slha_running_masses(const CSE6SSM_slha<T>&);
   template <class T>
   void extract_slha_susy_pars(const CSE6SSM_slha<T>&);
   template <class T>
   void extract_slha_soft_pars(const CSE6SSM_slha<T>&);
   template <class T>
   void extract_slha_pole_mixings(const CSE6SSM_slha<T>&);
   template <class T>
   void extract_slha_running_mixings(const CSE6SSM_slha<T>&);

   void write_slha_pole_masses_comment_line(std::ostream &) const;
   void write_slha_running_masses_comment_line(std::ostream &) const;
   void write_slha_susy_pars_comment_line(std::ostream &) const;
   void write_slha_soft_pars_comment_line(std::ostream &) const;
   void write_slha_pole_mixings_comment_line(std::ostream &) const;
   void write_slha_running_mixings_comment_line(std::ostream &) const;

   void write_slha_pole_masses_line(std::ostream &) const;
   void write_slha_running_masses_line(std::ostream &) const;
   void write_slha_susy_pars_line(std::ostream &) const;
   void write_slha_soft_pars_line(std::ostream &) const;
   void write_slha_pole_mixings_line(std::ostream &) const;
   void write_slha_running_mixings_line(std::ostream &) const;

private:
   struct TMass {
      std::string name;
      std::valarray<double> masses;
      TMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TParameter {
      std::string name;
      std::size_t mass_dimension;
      std::size_t rows;
      std::size_t cols;
      std::valarray<double> values;
      TParameter(const std::string& name_,
                 std::size_t mass_dimension_,
                 std::size_t rows_, std::size_t cols_,
                 const std::valarray<double>& values_)
         : name(name_)
         , mass_dimension(mass_dimension_)
         , rows(rows_)
         , cols(cols_)
         , values(values_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TMixing {
      std::string name;
      std::size_t dimension;
      bool is_real;
      std::valarray<double> mixings;
      TMixing(const std::string& name_, std::size_t dimension_,
              bool is_real_,
              const std::valarray<double>& mixings_)
         : name(name_)
         , dimension(dimension_)
         , is_real(is_real_)
         , mixings(mixings_)
         {}
   };

   typedef std::vector<TMass> TPoleMasses;
   typedef std::vector<TMass> TRunningMasses;
   typedef std::vector<TParameter> TSusyPars;
   typedef std::vector<TParameter> TSoftPars;
   typedef std::vector<TMixing> TMixings;

   TPoleMasses slha_pole_masses;
   TRunningMasses slha_running_masses;
   TSusyPars slha_susy_pars;
   TSoftPars slha_soft_pars;
   TMixings slha_pole_mixings;
   TMixings slha_running_mixings;

   CSE6SSM_input_parameters<Two_scale> slha_pole_masses_inputs;
   CSE6SSM_input_parameters<Two_scale> slha_running_masses_inputs;
   CSE6SSM_input_parameters<Two_scale> slha_susy_pars_inputs;
   CSE6SSM_input_parameters<Two_scale> slha_soft_pars_inputs;
   CSE6SSM_input_parameters<Two_scale> slha_pole_mixings_inputs;
   CSE6SSM_input_parameters<Two_scale> slha_running_mixings_inputs;

   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_pole_masses_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_running_masses_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_susy_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_soft_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_pole_mixings_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_running_mixings_problems;

   double high_scale;
   double susy_scale;
   double low_scale;
   unsigned width;
};

class CSE6SSM_semianalytic_slha_values_writer {
public:
   CSE6SSM_semianalytic_slha_values_writer();
   ~CSE6SSM_semianalytic_slha_values_writer() {}

   void set_high_scale(double high_scale_) { high_scale = high_scale_; }
   void set_susy_scale(double susy_scale_) { susy_scale = susy_scale_; }
   void set_low_scale(double low_scale_) { low_scale = low_scale_; }

   template <class T>
   void extract_slha_pole_masses(const CSE6SSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_running_masses(const CSE6SSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_susy_pars(const CSE6SSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_soft_pars(const CSE6SSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_pole_mixings(const CSE6SSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_running_mixings(const CSE6SSM_semianalytic_slha<T>&);

   void write_slha_pole_masses_comment_line(std::ostream &) const;
   void write_slha_running_masses_comment_line(std::ostream &) const;
   void write_slha_susy_pars_comment_line(std::ostream &) const;
   void write_slha_soft_pars_comment_line(std::ostream &) const;
   void write_slha_pole_mixings_comment_line(std::ostream &) const;
   void write_slha_running_mixings_comment_line(std::ostream &) const;

   void write_slha_pole_masses_line(std::ostream &) const;
   void write_slha_running_masses_line(std::ostream &) const;
   void write_slha_susy_pars_line(std::ostream &) const;
   void write_slha_soft_pars_line(std::ostream &) const;
   void write_slha_pole_mixings_line(std::ostream &) const;
   void write_slha_running_mixings_line(std::ostream &) const;

private:
   struct TMass {
      std::string name;
      std::valarray<double> masses;
      TMass(const std::string& name_, const std::valarray<double>& masses_)
         : name(name_)
         , masses(masses_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TParameter {
      std::string name;
      std::size_t mass_dimension;
      std::size_t rows;
      std::size_t cols;
      std::valarray<double> values;
      TParameter(const std::string& name_,
                 std::size_t mass_dimension_,
                 std::size_t rows_, std::size_t cols_,
                 const std::valarray<double>& values_)
         : name(name_)
         , mass_dimension(mass_dimension_)
         , rows(rows_)
         , cols(cols_)
         , values(values_)
         {}
   };

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TMixing {
      std::string name;
      std::size_t dimension;
      bool is_real;
      std::valarray<double> mixings;
      TMixing(const std::string& name_, std::size_t dimension_,
              bool is_real_,
              const std::valarray<double>& mixings_)
         : name(name_)
         , dimension(dimension_)
         , is_real(is_real_)
         , mixings(mixings_)
         {}
   };

   typedef std::vector<TMass> TPoleMasses;
   typedef std::vector<TMass> TRunningMasses;
   typedef std::vector<TParameter> TSusyPars;
   typedef std::vector<TParameter> TSoftPars;
   typedef std::vector<TMixing> TMixings;

   TPoleMasses slha_pole_masses;
   TRunningMasses slha_running_masses;
   TSusyPars slha_susy_pars;
   TSoftPars slha_soft_pars;
   TMixings slha_pole_mixings;
   TMixings slha_running_mixings;

   CSE6SSM_semianalytic_input_parameters<Two_scale> slha_pole_masses_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> slha_running_masses_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> slha_susy_pars_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> slha_soft_pars_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> slha_pole_mixings_inputs;
   CSE6SSM_semianalytic_input_parameters<Two_scale> slha_running_mixings_inputs;

   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_pole_masses_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_running_masses_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_susy_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_soft_pars_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_pole_mixings_problems;
   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> slha_running_mixings_problems;

   double high_scale;
   double susy_scale;
   double low_scale;
   double slha_pole_masses_m0Sqr;
   double slha_running_masses_m0Sqr;
   double slha_susy_pars_m0Sqr;
   double slha_soft_pars_m0Sqr;
   double slha_pole_mixings_m0Sqr;
   double slha_running_mixings_m0Sqr;
   unsigned width;
};

class CSE6SSM_semianalytic_coefficients_writer {
public:
   CSE6SSM_semianalytic_coefficients_writer();
   ~CSE6SSM_semianalytic_coefficients_writer() {}

   void set_high_scale(double high_scale_) { high_scale = high_scale_; }
   void set_susy_scale(double susy_scale_) { susy_scale = susy_scale_; }
   void set_low_scale(double low_scale_) { low_scale = low_scale_; }

   template <class T>
   void extract_coefficients(const CSE6SSM_semianalytic<T>&);

   void write_coefficients_comment_line(std::ostream &) const;

   void write_coefficients_line(std::ostream &) const;

private:

   // note: stored in column-major format, i.e.
   // elements are (0,0), (1,0), (2,0),...
   struct TCoefficient {
      std::string name;
      std::size_t mass_dimension;
      std::size_t rows;
      std::size_t cols;
      std::valarray<double> values;
      TCoefficient(const std::string& name_,
                   std::size_t mass_dimension_,
                   std::size_t rows_, std::size_t cols_,
                   const std::valarray<double>& values_)
         : name(name_)
         , mass_dimension(mass_dimension_)
         , rows(rows_)
         , cols(cols_)
         , values(values_)
         {
            assert(rows_ * cols_ == values_.size() &&
                   "TCoefficient: number of elements must match "
                   "number of values");
         }
   };

   typedef std::vector<TCoefficient> TCoefficients;

   TCoefficients coefficients;

   CSE6SSM_semianalytic_input_parameters<Two_scale> coefficients_inputs;

   Problems<CSE6SSM_info::NUMBER_OF_PARTICLES> coefficients_problems;

   double high_scale;
   double susy_scale;
   double low_scale;
   double coefficients_m0Sqr;
   unsigned width;
};

template <class T>
void CSE6SSM_pole_mass_writer::extract_pole_masses(const CSE6SSM<T>& model)
{
   pole_masses.clear();
   pole_masses_scale = model.get_scale();
   pole_masses_inputs = model.get_input();
   pole_masses_problems = model.get_problems();

   pole_masses.push_back(TPoleMass("MGlu", to_valarray(PHYSICAL(MGlu))));
   pole_masses.push_back(TPoleMass("MChaP", to_valarray(PHYSICAL(MChaP))));
   pole_masses.push_back(TPoleMass("MVZp", to_valarray(PHYSICAL(MVZp))));
   pole_masses.push_back(TPoleMass("MSd", to_valarray(PHYSICAL(MSd))));
   pole_masses.push_back(TPoleMass("MSv", to_valarray(PHYSICAL(MSv))));
   pole_masses.push_back(TPoleMass("MSu", to_valarray(PHYSICAL(MSu))));
   pole_masses.push_back(TPoleMass("MSe", to_valarray(PHYSICAL(MSe))));
   pole_masses.push_back(TPoleMass("MSDX", to_valarray(PHYSICAL(MSDX))));
   pole_masses.push_back(TPoleMass("Mhh", to_valarray(PHYSICAL(Mhh))));
   pole_masses.push_back(TPoleMass("MAh", to_valarray(PHYSICAL(MAh))));
   pole_masses.push_back(TPoleMass("MHpm", to_valarray(PHYSICAL(MHpm))));
   pole_masses.push_back(TPoleMass("MChi", to_valarray(PHYSICAL(MChi))));
   pole_masses.push_back(TPoleMass("MCha", to_valarray(PHYSICAL(MCha))));
   pole_masses.push_back(TPoleMass("MFDX", to_valarray(PHYSICAL(MFDX))));
   pole_masses.push_back(TPoleMass("MSHI0", to_valarray(PHYSICAL(MSHI0))));
   pole_masses.push_back(TPoleMass("MSHIPM", to_valarray(PHYSICAL(MSHIPM))));
   pole_masses.push_back(TPoleMass("MChaI", to_valarray(PHYSICAL(MChaI))));
   pole_masses.push_back(TPoleMass("MChiI", to_valarray(PHYSICAL(MChiI))));
   pole_masses.push_back(TPoleMass("MSHp0", to_valarray(PHYSICAL(MSHp0))));
   pole_masses.push_back(TPoleMass("MSHpp", to_valarray(PHYSICAL(MSHpp))));
   pole_masses.push_back(TPoleMass("MChiP", to_valarray(PHYSICAL(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      pole_masses.push_back(TPoleMass("MFd", to_valarray(PHYSICAL(MFd))));
      pole_masses.push_back(TPoleMass("MFe", to_valarray(PHYSICAL(MFe))));
      pole_masses.push_back(TPoleMass("MFu", to_valarray(PHYSICAL(MFu))));
      pole_masses.push_back(TPoleMass("MFv", to_valarray(PHYSICAL(MFv))));
      pole_masses.push_back(TPoleMass("MVG", to_valarray(PHYSICAL(MVG))));
      pole_masses.push_back(TPoleMass("MVP", to_valarray(PHYSICAL(MVP))));
      pole_masses.push_back(TPoleMass("MVWm", to_valarray(PHYSICAL(MVWm))));
      pole_masses.push_back(TPoleMass("MVZ", to_valarray(PHYSICAL(MVZ))));

   }
}

template <class T>
void CSE6SSM_semianalytic_pole_mass_writer::extract_pole_masses(
   const CSE6SSM_semianalytic<T>& model
   )
{
   pole_masses.clear();
   pole_masses_scale = model.get_scale();
   pole_masses_inputs = model.get_input();
   pole_masses_problems = model.get_problems();
   m0Sqr = model.get_ewsb_output_parameter(0);

   pole_masses.push_back(TPoleMass("MGlu", to_valarray(PHYSICAL(MGlu))));
   pole_masses.push_back(TPoleMass("MChaP", to_valarray(PHYSICAL(MChaP))));
   pole_masses.push_back(TPoleMass("MVZp", to_valarray(PHYSICAL(MVZp))));
   pole_masses.push_back(TPoleMass("MSd", to_valarray(PHYSICAL(MSd))));
   pole_masses.push_back(TPoleMass("MSv", to_valarray(PHYSICAL(MSv))));
   pole_masses.push_back(TPoleMass("MSu", to_valarray(PHYSICAL(MSu))));
   pole_masses.push_back(TPoleMass("MSe", to_valarray(PHYSICAL(MSe))));
   pole_masses.push_back(TPoleMass("MSDX", to_valarray(PHYSICAL(MSDX))));
   pole_masses.push_back(TPoleMass("Mhh", to_valarray(PHYSICAL(Mhh))));
   pole_masses.push_back(TPoleMass("MAh", to_valarray(PHYSICAL(MAh))));
   pole_masses.push_back(TPoleMass("MHpm", to_valarray(PHYSICAL(MHpm))));
   pole_masses.push_back(TPoleMass("MChi", to_valarray(PHYSICAL(MChi))));
   pole_masses.push_back(TPoleMass("MCha", to_valarray(PHYSICAL(MCha))));
   pole_masses.push_back(TPoleMass("MFDX", to_valarray(PHYSICAL(MFDX))));
   pole_masses.push_back(TPoleMass("MSHI0", to_valarray(PHYSICAL(MSHI0))));
   pole_masses.push_back(TPoleMass("MSHIPM", to_valarray(PHYSICAL(MSHIPM))));
   pole_masses.push_back(TPoleMass("MChaI", to_valarray(PHYSICAL(MChaI))));
   pole_masses.push_back(TPoleMass("MChiI", to_valarray(PHYSICAL(MChiI))));
   pole_masses.push_back(TPoleMass("MSHp0", to_valarray(PHYSICAL(MSHp0))));
   pole_masses.push_back(TPoleMass("MSHpp", to_valarray(PHYSICAL(MSHpp))));
   pole_masses.push_back(TPoleMass("MChiP", to_valarray(PHYSICAL(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      pole_masses.push_back(TPoleMass("MFd", to_valarray(PHYSICAL(MFd))));
      pole_masses.push_back(TPoleMass("MFe", to_valarray(PHYSICAL(MFe))));
      pole_masses.push_back(TPoleMass("MFu", to_valarray(PHYSICAL(MFu))));
      pole_masses.push_back(TPoleMass("MFv", to_valarray(PHYSICAL(MFv))));
      pole_masses.push_back(TPoleMass("MVG", to_valarray(PHYSICAL(MVG))));
      pole_masses.push_back(TPoleMass("MVP", to_valarray(PHYSICAL(MVP))));
      pole_masses.push_back(TPoleMass("MVWm", to_valarray(PHYSICAL(MVWm))));
      pole_masses.push_back(TPoleMass("MVZ", to_valarray(PHYSICAL(MVZ))));

   }
}

template <class T>
void CSE6SSM_drbar_values_writer::extract_drbar_masses(const CSE6SSM<T>& model)
{
   drbar_masses.clear();
   drbar_masses_scale = model.get_scale();
   drbar_masses_inputs = model.get_input();
   drbar_masses_problems = model.get_problems();

   drbar_masses.push_back(TMass("MVG", to_valarray(MODELPARAMETER(MVG))));
   drbar_masses.push_back(TMass("MGlu", to_valarray(MODELPARAMETER(MGlu))));
   drbar_masses.push_back(TMass("MFv", to_valarray(MODELPARAMETER(MFv))));
   drbar_masses.push_back(TMass("MChaP", to_valarray(MODELPARAMETER(MChaP))));
   drbar_masses.push_back(TMass("MVP", to_valarray(MODELPARAMETER(MVP))));
   drbar_masses.push_back(TMass("MVZ", to_valarray(MODELPARAMETER(MVZ))));
   drbar_masses.push_back(TMass("MVZp", to_valarray(MODELPARAMETER(MVZp))));
   drbar_masses.push_back(TMass("MSd", to_valarray(MODELPARAMETER(MSd))));
   drbar_masses.push_back(TMass("MSv", to_valarray(MODELPARAMETER(MSv))));
   drbar_masses.push_back(TMass("MSu", to_valarray(MODELPARAMETER(MSu))));
   drbar_masses.push_back(TMass("MSe", to_valarray(MODELPARAMETER(MSe))));
   drbar_masses.push_back(TMass("MSDX", to_valarray(MODELPARAMETER(MSDX))));
   drbar_masses.push_back(TMass("Mhh", to_valarray(MODELPARAMETER(Mhh))));
   drbar_masses.push_back(TMass("MAh", to_valarray(MODELPARAMETER(MAh))));
   drbar_masses.push_back(TMass("MHpm", to_valarray(MODELPARAMETER(MHpm))));
   drbar_masses.push_back(TMass("MChi", to_valarray(MODELPARAMETER(MChi))));
   drbar_masses.push_back(TMass("MCha", to_valarray(MODELPARAMETER(MCha))));
   drbar_masses.push_back(TMass("MFe", to_valarray(MODELPARAMETER(MFe))));
   drbar_masses.push_back(TMass("MFd", to_valarray(MODELPARAMETER(MFd))));
   drbar_masses.push_back(TMass("MFu", to_valarray(MODELPARAMETER(MFu))));
   drbar_masses.push_back(TMass("MFDX", to_valarray(MODELPARAMETER(MFDX))));
   drbar_masses.push_back(TMass("MSHI0", to_valarray(MODELPARAMETER(MSHI0))));
   drbar_masses.push_back(TMass("MSHIPM", to_valarray(MODELPARAMETER(MSHIPM))));
   drbar_masses.push_back(TMass("MChaI", to_valarray(MODELPARAMETER(MChaI))));
   drbar_masses.push_back(TMass("MChiI", to_valarray(MODELPARAMETER(MChiI))));
   drbar_masses.push_back(TMass("MSHp0", to_valarray(MODELPARAMETER(MSHp0))));
   drbar_masses.push_back(TMass("MSHpp", to_valarray(MODELPARAMETER(MSHpp))));
   drbar_masses.push_back(TMass("MChiP", to_valarray(MODELPARAMETER(MChiP))));
   drbar_masses.push_back(TMass("MVWm", to_valarray(MODELPARAMETER(MVWm))));
}

template <class T>
void CSE6SSM_semianalytic_drbar_values_writer::extract_drbar_masses(
   const CSE6SSM_semianalytic<T>& model
   )
{
   drbar_masses.clear();
   drbar_masses_scale = model.get_scale();
   drbar_masses_inputs = model.get_input();
   drbar_masses_problems = model.get_problems();

   drbar_masses.push_back(TMass("MVG", to_valarray(MODELPARAMETER(MVG))));
   drbar_masses.push_back(TMass("MGlu", to_valarray(MODELPARAMETER(MGlu))));
   drbar_masses.push_back(TMass("MFv", to_valarray(MODELPARAMETER(MFv))));
   drbar_masses.push_back(TMass("MChaP", to_valarray(MODELPARAMETER(MChaP))));
   drbar_masses.push_back(TMass("MVP", to_valarray(MODELPARAMETER(MVP))));
   drbar_masses.push_back(TMass("MVZ", to_valarray(MODELPARAMETER(MVZ))));
   drbar_masses.push_back(TMass("MVZp", to_valarray(MODELPARAMETER(MVZp))));
   drbar_masses.push_back(TMass("MSd", to_valarray(MODELPARAMETER(MSd))));
   drbar_masses.push_back(TMass("MSv", to_valarray(MODELPARAMETER(MSv))));
   drbar_masses.push_back(TMass("MSu", to_valarray(MODELPARAMETER(MSu))));
   drbar_masses.push_back(TMass("MSe", to_valarray(MODELPARAMETER(MSe))));
   drbar_masses.push_back(TMass("MSDX", to_valarray(MODELPARAMETER(MSDX))));
   drbar_masses.push_back(TMass("Mhh", to_valarray(MODELPARAMETER(Mhh))));
   drbar_masses.push_back(TMass("MAh", to_valarray(MODELPARAMETER(MAh))));
   drbar_masses.push_back(TMass("MHpm", to_valarray(MODELPARAMETER(MHpm))));
   drbar_masses.push_back(TMass("MChi", to_valarray(MODELPARAMETER(MChi))));
   drbar_masses.push_back(TMass("MCha", to_valarray(MODELPARAMETER(MCha))));
   drbar_masses.push_back(TMass("MFe", to_valarray(MODELPARAMETER(MFe))));
   drbar_masses.push_back(TMass("MFd", to_valarray(MODELPARAMETER(MFd))));
   drbar_masses.push_back(TMass("MFu", to_valarray(MODELPARAMETER(MFu))));
   drbar_masses.push_back(TMass("MFDX", to_valarray(MODELPARAMETER(MFDX))));
   drbar_masses.push_back(TMass("MSHI0", to_valarray(MODELPARAMETER(MSHI0))));
   drbar_masses.push_back(TMass("MSHIPM", to_valarray(MODELPARAMETER(MSHIPM))));
   drbar_masses.push_back(TMass("MChaI", to_valarray(MODELPARAMETER(MChaI))));
   drbar_masses.push_back(TMass("MChiI", to_valarray(MODELPARAMETER(MChiI))));
   drbar_masses.push_back(TMass("MSHp0", to_valarray(MODELPARAMETER(MSHp0))));
   drbar_masses.push_back(TMass("MSHpp", to_valarray(MODELPARAMETER(MSHpp))));
   drbar_masses.push_back(TMass("MChiP", to_valarray(MODELPARAMETER(MChiP))));
   drbar_masses.push_back(TMass("MVWm", to_valarray(MODELPARAMETER(MVWm))));
}

template <class T>
void CSE6SSM_drbar_values_writer::extract_drbar_susy_pars(
   const CSE6SSM<T>& model
   )
{
   drbar_susy_pars.clear();
   drbar_susy_pars_scale = model.get_scale();
   drbar_susy_pars_inputs = model.get_input();
   drbar_susy_pars_problems = model.get_problems();

   drbar_susy_pars.push_back(TParameter("Yd", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Yd))));
   drbar_susy_pars.push_back(TParameter("hE", 0, 3, 2,
                                        to_valarray(MODELPARAMETER(hE))));
   drbar_susy_pars.push_back(TParameter("Ye", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Ye))));
   drbar_susy_pars.push_back(TParameter("SigmaL", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(SigmaL))));
   drbar_susy_pars.push_back(TParameter("KappaPr", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(KappaPr))));
   drbar_susy_pars.push_back(TParameter("Sigmax", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(Sigmax))));
   drbar_susy_pars.push_back(TParameter("gD", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(gD))));
   drbar_susy_pars.push_back(TParameter("Kappa", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Kappa))));
   drbar_susy_pars.push_back(TParameter("Lambda12", 0, 2, 2,
                                        to_valarray(MODELPARAMETER(Lambda12))));
   drbar_susy_pars.push_back(TParameter("Lambdax", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(Lambdax))));
   drbar_susy_pars.push_back(TParameter("fu", 0, 3, 2,
                                        to_valarray(MODELPARAMETER(fu))));
   drbar_susy_pars.push_back(TParameter("fd", 0, 3, 2,
                                        to_valarray(MODELPARAMETER(fd))));
   drbar_susy_pars.push_back(TParameter("Yu", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Yu))));
   drbar_susy_pars.push_back(TParameter("MuPr", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MuPr))));
   drbar_susy_pars.push_back(TParameter("MuPhi", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MuPhi))));
   drbar_susy_pars.push_back(TParameter("XiF", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(XiF))));
   drbar_susy_pars.push_back(TParameter("g1", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g1))));
   drbar_susy_pars.push_back(TParameter("g2", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g2))));
   drbar_susy_pars.push_back(TParameter("g3", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g3))));
   drbar_susy_pars.push_back(TParameter("g1p", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g1p))));
   drbar_susy_pars.push_back(TParameter("vd", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vd))));
   drbar_susy_pars.push_back(TParameter("vu", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vu))));
   drbar_susy_pars.push_back(TParameter("vs", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vs))));
   drbar_susy_pars.push_back(TParameter("vsb", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vsb))));
   drbar_susy_pars.push_back(TParameter("vphi", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vphi))));
}

template <class T>
void CSE6SSM_semianalytic_drbar_values_writer::extract_drbar_susy_pars(
   const CSE6SSM_semianalytic<T>& model
   )
{
   drbar_susy_pars.clear();
   drbar_susy_pars_scale = model.get_scale();
   drbar_susy_pars_inputs = model.get_input();
   drbar_susy_pars_problems = model.get_problems();

   drbar_susy_pars.push_back(TParameter("Yd", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Yd))));
   drbar_susy_pars.push_back(TParameter("hE", 0, 3, 2,
                                        to_valarray(MODELPARAMETER(hE))));
   drbar_susy_pars.push_back(TParameter("Ye", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Ye))));
   drbar_susy_pars.push_back(TParameter("SigmaL", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(SigmaL))));
   drbar_susy_pars.push_back(TParameter("KappaPr", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(KappaPr))));
   drbar_susy_pars.push_back(TParameter("Sigmax", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(Sigmax))));
   drbar_susy_pars.push_back(TParameter("gD", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(gD))));
   drbar_susy_pars.push_back(TParameter("Kappa", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Kappa))));
   drbar_susy_pars.push_back(TParameter("Lambda12", 0, 2, 2,
                                        to_valarray(MODELPARAMETER(Lambda12))));
   drbar_susy_pars.push_back(TParameter("Lambdax", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(Lambdax))));
   drbar_susy_pars.push_back(TParameter("fu", 0, 3, 2,
                                        to_valarray(MODELPARAMETER(fu))));
   drbar_susy_pars.push_back(TParameter("fd", 0, 3, 2,
                                        to_valarray(MODELPARAMETER(fd))));
   drbar_susy_pars.push_back(TParameter("Yu", 0, 3, 3,
                                        to_valarray(MODELPARAMETER(Yu))));
   drbar_susy_pars.push_back(TParameter("MuPr", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MuPr))));
   drbar_susy_pars.push_back(TParameter("MuPhi", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MuPhi))));
   drbar_susy_pars.push_back(TParameter("XiF", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(XiF))));
   drbar_susy_pars.push_back(TParameter("g1", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g1))));
   drbar_susy_pars.push_back(TParameter("g2", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g2))));
   drbar_susy_pars.push_back(TParameter("g3", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g3))));
   drbar_susy_pars.push_back(TParameter("g1p", 0, 1, 1,
                                        to_valarray(MODELPARAMETER(g1p))));
   drbar_susy_pars.push_back(TParameter("vd", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vd))));
   drbar_susy_pars.push_back(TParameter("vu", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vu))));
   drbar_susy_pars.push_back(TParameter("vs", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vs))));
   drbar_susy_pars.push_back(TParameter("vsb", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vsb))));
   drbar_susy_pars.push_back(TParameter("vphi", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(vphi))));
}

template <class T>
void CSE6SSM_drbar_values_writer::extract_drbar_soft_pars(
   const CSE6SSM<T>& model
   )
{
   drbar_soft_pars.clear();
   drbar_soft_pars_scale = model.get_scale();
   drbar_soft_pars_inputs = model.get_input();
   drbar_soft_pars_problems = model.get_problems();

   drbar_soft_pars.push_back(TParameter("TYd", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TYd))));
   drbar_soft_pars.push_back(TParameter("ThE", 1, 3, 2,
                                        to_valarray(MODELPARAMETER(ThE))));
   drbar_soft_pars.push_back(TParameter("TYe", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TYe))));
   drbar_soft_pars.push_back(TParameter("TSigmaL", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(TSigmaL))));
   drbar_soft_pars.push_back(TParameter("TKappaPr", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(TKappaPr))));
   drbar_soft_pars.push_back(TParameter("TSigmax", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(TSigmax))));
   drbar_soft_pars.push_back(TParameter("TgD", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TgD))));
   drbar_soft_pars.push_back(TParameter("TKappa", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TKappa))));
   drbar_soft_pars.push_back(
      TParameter("TLambda12", 1, 2, 2, to_valarray(MODELPARAMETER(TLambda12))));
   drbar_soft_pars.push_back(
      TParameter("TLambdax", 1, 1, 1, to_valarray(MODELPARAMETER(TLambdax))));
   drbar_soft_pars.push_back(TParameter("Tfu", 1, 3, 2,
                                        to_valarray(MODELPARAMETER(Tfu))));
   drbar_soft_pars.push_back(TParameter("Tfd", 1, 3, 2,
                                        to_valarray(MODELPARAMETER(Tfd))));
   drbar_soft_pars.push_back(TParameter("TYu", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TYu))));
   drbar_soft_pars.push_back(TParameter("BMuPr", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(BMuPr))));
   drbar_soft_pars.push_back(TParameter("BMuPhi", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(BMuPhi))));
   drbar_soft_pars.push_back(TParameter("LXiF", 3, 1, 1,
                                        to_valarray(MODELPARAMETER(LXiF))));
   drbar_soft_pars.push_back(TParameter("mq2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mq2))));
   drbar_soft_pars.push_back(TParameter("ml2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(ml2))));
   drbar_soft_pars.push_back(TParameter("mHd2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHd2))));
   drbar_soft_pars.push_back(TParameter("mHu2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHu2))));
   drbar_soft_pars.push_back(TParameter("md2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(md2))));
   drbar_soft_pars.push_back(TParameter("mu2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mu2))));
   drbar_soft_pars.push_back(TParameter("me2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(me2))));
   drbar_soft_pars.push_back(TParameter("ms2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(ms2))));
   drbar_soft_pars.push_back(TParameter("msbar2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(msbar2))));
   drbar_soft_pars.push_back(TParameter("mH1I2", 2, 2, 2,
                                        to_valarray(MODELPARAMETER(mH1I2))));
   drbar_soft_pars.push_back(TParameter("mH2I2", 2, 2, 2,
                                        to_valarray(MODELPARAMETER(mH2I2))));
   drbar_soft_pars.push_back(TParameter("mSI2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mSI2))));
   drbar_soft_pars.push_back(TParameter("mDx2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mDx2))));
   drbar_soft_pars.push_back(TParameter("mDxbar2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mDxbar2))));
   drbar_soft_pars.push_back(TParameter("mHp2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHp2))));
   drbar_soft_pars.push_back(TParameter("mHpbar2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHpbar2))));
   drbar_soft_pars.push_back(TParameter("mphi2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mphi2))));
   drbar_soft_pars.push_back(TParameter("MassB", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassB))));
   drbar_soft_pars.push_back(TParameter("MassWB", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassWB))));
   drbar_soft_pars.push_back(TParameter("MassG", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassG))));
   drbar_soft_pars.push_back(TParameter("MassBp", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassBp))));
}

template <class T>
void CSE6SSM_semianalytic_drbar_values_writer::extract_drbar_soft_pars(
   const CSE6SSM_semianalytic<T>& model
   )
{
   drbar_soft_pars.clear();
   drbar_soft_pars_scale = model.get_scale();
   drbar_soft_pars_inputs = model.get_input();
   drbar_soft_pars_problems = model.get_problems();

   drbar_soft_pars.push_back(TParameter("TYd", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TYd))));
   drbar_soft_pars.push_back(TParameter("ThE", 1, 3, 2,
                                        to_valarray(MODELPARAMETER(ThE))));
   drbar_soft_pars.push_back(TParameter("TYe", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TYe))));
   drbar_soft_pars.push_back(TParameter("TSigmaL", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(TSigmaL))));
   drbar_soft_pars.push_back(TParameter("TKappaPr", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(TKappaPr))));
   drbar_soft_pars.push_back(TParameter("TSigmax", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(TSigmax))));
   drbar_soft_pars.push_back(TParameter("TgD", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TgD))));
   drbar_soft_pars.push_back(TParameter("TKappa", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TKappa))));
   drbar_soft_pars.push_back(
      TParameter("TLambda12", 1, 2, 2, to_valarray(MODELPARAMETER(TLambda12))));
   drbar_soft_pars.push_back(
      TParameter("TLambdax", 1, 1, 1, to_valarray(MODELPARAMETER(TLambdax))));
   drbar_soft_pars.push_back(TParameter("Tfu", 1, 3, 2,
                                        to_valarray(MODELPARAMETER(Tfu))));
   drbar_soft_pars.push_back(TParameter("Tfd", 1, 3, 2,
                                        to_valarray(MODELPARAMETER(Tfd))));
   drbar_soft_pars.push_back(TParameter("TYu", 1, 3, 3,
                                        to_valarray(MODELPARAMETER(TYu))));
   drbar_soft_pars.push_back(TParameter("BMuPr", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(BMuPr))));
   drbar_soft_pars.push_back(TParameter("BMuPhi", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(BMuPhi))));
   drbar_soft_pars.push_back(TParameter("LXiF", 3, 1, 1,
                                        to_valarray(MODELPARAMETER(LXiF))));
   drbar_soft_pars.push_back(TParameter("mq2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mq2))));
   drbar_soft_pars.push_back(TParameter("ml2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(ml2))));
   drbar_soft_pars.push_back(TParameter("mHd2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHd2))));
   drbar_soft_pars.push_back(TParameter("mHu2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHu2))));
   drbar_soft_pars.push_back(TParameter("md2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(md2))));
   drbar_soft_pars.push_back(TParameter("mu2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mu2))));
   drbar_soft_pars.push_back(TParameter("me2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(me2))));
   drbar_soft_pars.push_back(TParameter("ms2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(ms2))));
   drbar_soft_pars.push_back(TParameter("msbar2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(msbar2))));
   drbar_soft_pars.push_back(TParameter("mH1I2", 2, 2, 2,
                                        to_valarray(MODELPARAMETER(mH1I2))));
   drbar_soft_pars.push_back(TParameter("mH2I2", 2, 2, 2,
                                        to_valarray(MODELPARAMETER(mH2I2))));
   drbar_soft_pars.push_back(TParameter("mSI2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mSI2))));
   drbar_soft_pars.push_back(TParameter("mDx2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mDx2))));
   drbar_soft_pars.push_back(TParameter("mDxbar2", 2, 3, 3,
                                        to_valarray(MODELPARAMETER(mDxbar2))));
   drbar_soft_pars.push_back(TParameter("mHp2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHp2))));
   drbar_soft_pars.push_back(TParameter("mHpbar2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mHpbar2))));
   drbar_soft_pars.push_back(TParameter("mphi2", 2, 1, 1,
                                        to_valarray(MODELPARAMETER(mphi2))));
   drbar_soft_pars.push_back(TParameter("MassB", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassB))));
   drbar_soft_pars.push_back(TParameter("MassWB", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassWB))));
   drbar_soft_pars.push_back(TParameter("MassG", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassG))));
   drbar_soft_pars.push_back(TParameter("MassBp", 1, 1, 1,
                                        to_valarray(MODELPARAMETER(MassBp))));
}

template <class T>
void CSE6SSM_drbar_values_writer::extract_drbar_mixings(const CSE6SSM<T>& model)
{
   drbar_mixings.clear();
   drbar_mixings_scale = model.get_scale();
   drbar_mixings_inputs = model.get_input();
   drbar_mixings_problems = model.get_problems();

   drbar_mixings.push_back(TMixing("ZD", 6, true,
                                   to_valarray(MODELPARAMETER(ZD))));
   drbar_mixings.push_back(TMixing("ZV", 3, true,
                                   to_valarray(MODELPARAMETER(ZV))));
   drbar_mixings.push_back(TMixing("ZU", 6, true,
                                   to_valarray(MODELPARAMETER(ZU))));
   drbar_mixings.push_back(TMixing("ZE", 6, true,
                                   to_valarray(MODELPARAMETER(ZE))));
   drbar_mixings.push_back(TMixing("ZDX", 6, true,
                                   to_valarray(MODELPARAMETER(ZDX))));
   drbar_mixings.push_back(TMixing("ZH", 5, true,
                                   to_valarray(MODELPARAMETER(ZH))));
   drbar_mixings.push_back(TMixing("ZA", 5, true,
                                   to_valarray(MODELPARAMETER(ZA))));
   drbar_mixings.push_back(TMixing("ZP", 2, true,
                                   to_valarray(MODELPARAMETER(ZP))));
   drbar_mixings.push_back(TMixing("ZN", 8, false,
                                   to_valarray(MODELPARAMETER(ZN))));
   drbar_mixings.push_back(TMixing("UM", 2, false,
                                   to_valarray(MODELPARAMETER(UM))));
   drbar_mixings.push_back(TMixing("UP", 2, false,
                                   to_valarray(MODELPARAMETER(UP))));
   drbar_mixings.push_back(TMixing("ZEL", 3, false,
                                   to_valarray(MODELPARAMETER(ZEL))));
   drbar_mixings.push_back(TMixing("ZER", 3, false,
                                   to_valarray(MODELPARAMETER(ZER))));
   drbar_mixings.push_back(TMixing("ZDL", 3, false,
                                   to_valarray(MODELPARAMETER(ZDL))));
   drbar_mixings.push_back(TMixing("ZDR", 3, false,
                                   to_valarray(MODELPARAMETER(ZDR))));
   drbar_mixings.push_back(TMixing("ZUL", 3, false,
                                   to_valarray(MODELPARAMETER(ZUL))));
   drbar_mixings.push_back(TMixing("ZUR", 3, false,
                                   to_valarray(MODELPARAMETER(ZUR))));
   drbar_mixings.push_back(TMixing("ZDXL", 3, false,
                                   to_valarray(MODELPARAMETER(ZDXL))));
   drbar_mixings.push_back(TMixing("ZDXR", 3, false,
                                   to_valarray(MODELPARAMETER(ZDXR))));
   drbar_mixings.push_back(TMixing("UHI0", 7, true,
                                   to_valarray(MODELPARAMETER(UHI0))));
   drbar_mixings.push_back(TMixing("UHIPM", 4, true,
                                   to_valarray(MODELPARAMETER(UHIPM))));
   drbar_mixings.push_back(TMixing("ZMI", 2, false,
                                   to_valarray(MODELPARAMETER(ZMI))));
   drbar_mixings.push_back(TMixing("ZPI", 2, false,
                                   to_valarray(MODELPARAMETER(ZPI))));
   drbar_mixings.push_back(TMixing("ZNI", 7, false,
                                   to_valarray(MODELPARAMETER(ZNI))));
   drbar_mixings.push_back(TMixing("UHp0", 2, true,
                                   to_valarray(MODELPARAMETER(UHp0))));
   drbar_mixings.push_back(TMixing("UHpp", 2, true,
                                   to_valarray(MODELPARAMETER(UHpp))));
   drbar_mixings.push_back(TMixing("ZNp", 2, false,
                                   to_valarray(MODELPARAMETER(ZNp))));
}

template <class T>
void CSE6SSM_semianalytic_drbar_values_writer::extract_drbar_mixings(
   const CSE6SSM_semianalytic<T>& model
   )
{
   drbar_mixings.clear();
   drbar_mixings_scale = model.get_scale();
   drbar_mixings_inputs = model.get_input();
   drbar_mixings_problems = model.get_problems();

   drbar_mixings.push_back(TMixing("ZD", 6, true,
                                   to_valarray(MODELPARAMETER(ZD))));
   drbar_mixings.push_back(TMixing("ZV", 3, true,
                                   to_valarray(MODELPARAMETER(ZV))));
   drbar_mixings.push_back(TMixing("ZU", 6, true,
                                   to_valarray(MODELPARAMETER(ZU))));
   drbar_mixings.push_back(TMixing("ZE", 6, true,
                                   to_valarray(MODELPARAMETER(ZE))));
   drbar_mixings.push_back(TMixing("ZDX", 6, true,
                                   to_valarray(MODELPARAMETER(ZDX))));
   drbar_mixings.push_back(TMixing("ZH", 5, true,
                                   to_valarray(MODELPARAMETER(ZH))));
   drbar_mixings.push_back(TMixing("ZA", 5, true,
                                   to_valarray(MODELPARAMETER(ZA))));
   drbar_mixings.push_back(TMixing("ZP", 2, true,
                                   to_valarray(MODELPARAMETER(ZP))));
   drbar_mixings.push_back(TMixing("ZN", 8, false,
                                   to_valarray(MODELPARAMETER(ZN))));
   drbar_mixings.push_back(TMixing("UM", 2, false,
                                   to_valarray(MODELPARAMETER(UM))));
   drbar_mixings.push_back(TMixing("UP", 2, false,
                                   to_valarray(MODELPARAMETER(UP))));
   drbar_mixings.push_back(TMixing("ZEL", 3, false,
                                   to_valarray(MODELPARAMETER(ZEL))));
   drbar_mixings.push_back(TMixing("ZER", 3, false,
                                   to_valarray(MODELPARAMETER(ZER))));
   drbar_mixings.push_back(TMixing("ZDL", 3, false,
                                   to_valarray(MODELPARAMETER(ZDL))));
   drbar_mixings.push_back(TMixing("ZDR", 3, false,
                                   to_valarray(MODELPARAMETER(ZDR))));
   drbar_mixings.push_back(TMixing("ZUL", 3, false,
                                   to_valarray(MODELPARAMETER(ZUL))));
   drbar_mixings.push_back(TMixing("ZUR", 3, false,
                                   to_valarray(MODELPARAMETER(ZUR))));
   drbar_mixings.push_back(TMixing("ZDXL", 3, false,
                                   to_valarray(MODELPARAMETER(ZDXL))));
   drbar_mixings.push_back(TMixing("ZDXR", 3, false,
                                   to_valarray(MODELPARAMETER(ZDXR))));
   drbar_mixings.push_back(TMixing("UHI0", 7, true,
                                   to_valarray(MODELPARAMETER(UHI0))));
   drbar_mixings.push_back(TMixing("UHIPM", 4, true,
                                   to_valarray(MODELPARAMETER(UHIPM))));
   drbar_mixings.push_back(TMixing("ZMI", 2, false,
                                   to_valarray(MODELPARAMETER(ZMI))));
   drbar_mixings.push_back(TMixing("ZPI", 2, false,
                                   to_valarray(MODELPARAMETER(ZPI))));
   drbar_mixings.push_back(TMixing("ZNI", 7, false,
                                   to_valarray(MODELPARAMETER(ZNI))));
   drbar_mixings.push_back(TMixing("UHp0", 2, true,
                                   to_valarray(MODELPARAMETER(UHp0))));
   drbar_mixings.push_back(TMixing("UHpp", 2, true,
                                   to_valarray(MODELPARAMETER(UHpp))));
   drbar_mixings.push_back(TMixing("ZNp", 2, false,
                                   to_valarray(MODELPARAMETER(ZNp))));
}

template <class T>
void CSE6SSM_slha_values_writer::extract_slha_pole_masses(
   const CSE6SSM_slha<T>& model
   )
{
   slha_pole_masses.clear();
   slha_pole_masses_inputs = model.get_input();
   slha_pole_masses_problems = model.get_problems();

   slha_pole_masses.push_back(TMass("MGlu", to_valarray(PHYSICALSLHA(MGlu))));
   slha_pole_masses.push_back(TMass("MChaP", to_valarray(PHYSICALSLHA(MChaP))));
   slha_pole_masses.push_back(TMass("MVZp", to_valarray(PHYSICALSLHA(MVZp))));
   slha_pole_masses.push_back(TMass("MSd", to_valarray(PHYSICALSLHA(MSd))));
   slha_pole_masses.push_back(TMass("MSv", to_valarray(PHYSICALSLHA(MSv))));
   slha_pole_masses.push_back(TMass("MSu", to_valarray(PHYSICALSLHA(MSu))));
   slha_pole_masses.push_back(TMass("MSe", to_valarray(PHYSICALSLHA(MSe))));
   slha_pole_masses.push_back(TMass("MSDX", to_valarray(PHYSICALSLHA(MSDX))));
   slha_pole_masses.push_back(TMass("Mhh", to_valarray(PHYSICALSLHA(Mhh))));
   slha_pole_masses.push_back(TMass("MAh", to_valarray(PHYSICALSLHA(MAh))));
   slha_pole_masses.push_back(TMass("MHpm", to_valarray(PHYSICALSLHA(MHpm))));
   slha_pole_masses.push_back(TMass("MChi", to_valarray(PHYSICALSLHA(MChi))));
   slha_pole_masses.push_back(TMass("MCha", to_valarray(PHYSICALSLHA(MCha))));
   slha_pole_masses.push_back(TMass("MFDX", to_valarray(PHYSICALSLHA(MFDX))));
   slha_pole_masses.push_back(TMass("MSHI0", to_valarray(PHYSICALSLHA(MSHI0))));
   slha_pole_masses.push_back(TMass("MSHIPM",
                                    to_valarray(PHYSICALSLHA(MSHIPM))));
   slha_pole_masses.push_back(TMass("MChaI", to_valarray(PHYSICALSLHA(MChaI))));
   slha_pole_masses.push_back(TMass("MChiI", to_valarray(PHYSICALSLHA(MChiI))));
   slha_pole_masses.push_back(TMass("MSHp0", to_valarray(PHYSICALSLHA(MSHp0))));
   slha_pole_masses.push_back(TMass("MSHpp", to_valarray(PHYSICALSLHA(MSHpp))));
   slha_pole_masses.push_back(TMass("MChiP", to_valarray(PHYSICALSLHA(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_pole_masses.push_back(TMass("MFd", to_valarray(PHYSICALSLHA(MFd))));
      slha_pole_masses.push_back(TMass("MFe", to_valarray(PHYSICALSLHA(MFe))));
      slha_pole_masses.push_back(TMass("MFu", to_valarray(PHYSICALSLHA(MFu))));
      slha_pole_masses.push_back(TMass("MFv", to_valarray(PHYSICALSLHA(MFv))));
      slha_pole_masses.push_back(TMass("MVG", to_valarray(PHYSICALSLHA(MVG))));
      slha_pole_masses.push_back(TMass("MVP", to_valarray(PHYSICALSLHA(MVP))));
      slha_pole_masses.push_back(TMass("MVWm",
                                       to_valarray(PHYSICALSLHA(MVWm))));
      slha_pole_masses.push_back(TMass("MVZ", to_valarray(PHYSICALSLHA(MVZ))));

   }
}

template <class T>
void CSE6SSM_semianalytic_slha_values_writer::extract_slha_pole_masses(
   const CSE6SSM_semianalytic_slha<T>& model
   )
{
   slha_pole_masses.clear();
   slha_pole_masses_inputs = model.get_input();
   slha_pole_masses_problems = model.get_problems();
   slha_pole_masses_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_pole_masses.push_back(TMass("MGlu", to_valarray(PHYSICALSLHA(MGlu))));
   slha_pole_masses.push_back(TMass("MChaP", to_valarray(PHYSICALSLHA(MChaP))));
   slha_pole_masses.push_back(TMass("MVZp", to_valarray(PHYSICALSLHA(MVZp))));
   slha_pole_masses.push_back(TMass("MSd", to_valarray(PHYSICALSLHA(MSd))));
   slha_pole_masses.push_back(TMass("MSv", to_valarray(PHYSICALSLHA(MSv))));
   slha_pole_masses.push_back(TMass("MSu", to_valarray(PHYSICALSLHA(MSu))));
   slha_pole_masses.push_back(TMass("MSe", to_valarray(PHYSICALSLHA(MSe))));
   slha_pole_masses.push_back(TMass("MSDX", to_valarray(PHYSICALSLHA(MSDX))));
   slha_pole_masses.push_back(TMass("Mhh", to_valarray(PHYSICALSLHA(Mhh))));
   slha_pole_masses.push_back(TMass("MAh", to_valarray(PHYSICALSLHA(MAh))));
   slha_pole_masses.push_back(TMass("MHpm", to_valarray(PHYSICALSLHA(MHpm))));
   slha_pole_masses.push_back(TMass("MChi", to_valarray(PHYSICALSLHA(MChi))));
   slha_pole_masses.push_back(TMass("MCha", to_valarray(PHYSICALSLHA(MCha))));
   slha_pole_masses.push_back(TMass("MFDX", to_valarray(PHYSICALSLHA(MFDX))));
   slha_pole_masses.push_back(TMass("MSHI0", to_valarray(PHYSICALSLHA(MSHI0))));
   slha_pole_masses.push_back(TMass("MSHIPM",
                                    to_valarray(PHYSICALSLHA(MSHIPM))));
   slha_pole_masses.push_back(TMass("MChaI", to_valarray(PHYSICALSLHA(MChaI))));
   slha_pole_masses.push_back(TMass("MChiI", to_valarray(PHYSICALSLHA(MChiI))));
   slha_pole_masses.push_back(TMass("MSHp0", to_valarray(PHYSICALSLHA(MSHp0))));
   slha_pole_masses.push_back(TMass("MSHpp", to_valarray(PHYSICALSLHA(MSHpp))));
   slha_pole_masses.push_back(TMass("MChiP", to_valarray(PHYSICALSLHA(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_pole_masses.push_back(TMass("MFd", to_valarray(PHYSICALSLHA(MFd))));
      slha_pole_masses.push_back(TMass("MFe", to_valarray(PHYSICALSLHA(MFe))));
      slha_pole_masses.push_back(TMass("MFu", to_valarray(PHYSICALSLHA(MFu))));
      slha_pole_masses.push_back(TMass("MFv", to_valarray(PHYSICALSLHA(MFv))));
      slha_pole_masses.push_back(TMass("MVG", to_valarray(PHYSICALSLHA(MVG))));
      slha_pole_masses.push_back(TMass("MVP", to_valarray(PHYSICALSLHA(MVP))));
      slha_pole_masses.push_back(TMass("MVWm",
                                       to_valarray(PHYSICALSLHA(MVWm))));
      slha_pole_masses.push_back(TMass("MVZ", to_valarray(PHYSICALSLHA(MVZ))));

   }
}

template <class T>
void CSE6SSM_slha_values_writer::extract_slha_running_masses(
   const CSE6SSM_slha<T>& model
   )
{
   slha_running_masses.clear();
   slha_running_masses_inputs = model.get_input();
   slha_running_masses_problems = model.get_problems();

   slha_running_masses.push_back(TMass("DRbarMGlu",
                                       to_valarray(DRBARSLHA(MGlu))));
   slha_running_masses.push_back(TMass("DRbarMChaP",
                                       to_valarray(DRBARSLHA(MChaP))));
   slha_running_masses.push_back(TMass("DRbarMVZp",
                                       to_valarray(DRBARSLHA(MVZp))));
   slha_running_masses.push_back(TMass("DRbarMSd",
                                       to_valarray(DRBARSLHA(MSd))));
   slha_running_masses.push_back(TMass("DRbarMSv",
                                       to_valarray(DRBARSLHA(MSv))));
   slha_running_masses.push_back(TMass("DRbarMSu",
                                       to_valarray(DRBARSLHA(MSu))));
   slha_running_masses.push_back(TMass("DRbarMSe",
                                       to_valarray(DRBARSLHA(MSe))));
   slha_running_masses.push_back(TMass("DRbarMSDX",
                                       to_valarray(DRBARSLHA(MSDX))));
   slha_running_masses.push_back(TMass("DRbarMhh",
                                       to_valarray(DRBARSLHA(Mhh))));
   slha_running_masses.push_back(TMass("DRbarMAh",
                                       to_valarray(DRBARSLHA(MAh))));
   slha_running_masses.push_back(TMass("DRbarMHpm",
                                       to_valarray(DRBARSLHA(MHpm))));
   slha_running_masses.push_back(TMass("DRbarMChi",
                                       to_valarray(DRBARSLHA(MChi))));
   slha_running_masses.push_back(TMass("DRbarMCha",
                                       to_valarray(DRBARSLHA(MCha))));
   slha_running_masses.push_back(TMass("DRbarMFDX",
                                       to_valarray(DRBARSLHA(MFDX))));
   slha_running_masses.push_back(TMass("DRbarMSHI0",
                                       to_valarray(DRBARSLHA(MSHI0))));
   slha_running_masses.push_back(TMass("DRbarMSHIPM",
                                       to_valarray(DRBARSLHA(MSHIPM))));
   slha_running_masses.push_back(TMass("DRbarMChaI",
                                       to_valarray(DRBARSLHA(MChaI))));
   slha_running_masses.push_back(TMass("DRbarMChiI",
                                       to_valarray(DRBARSLHA(MChiI))));
   slha_running_masses.push_back(TMass("DRbarMSHp0",
                                       to_valarray(DRBARSLHA(MSHp0))));
   slha_running_masses.push_back(TMass("DRbarMSHpp",
                                       to_valarray(DRBARSLHA(MSHpp))));
   slha_running_masses.push_back(TMass("DRbarMChiP",
                                       to_valarray(DRBARSLHA(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_running_masses.push_back(TMass("DRbarMFd",
                                          to_valarray(DRBARSLHA(MFd))));
      slha_running_masses.push_back(TMass("DRbarMFe",
                                          to_valarray(DRBARSLHA(MFe))));
      slha_running_masses.push_back(TMass("DRbarMFu",
                                          to_valarray(DRBARSLHA(MFu))));
      slha_running_masses.push_back(TMass("DRbarMFv",
                                          to_valarray(DRBARSLHA(MFv))));
      slha_running_masses.push_back(TMass("DRbarMVG",
                                          to_valarray(DRBARSLHA(MVG))));
      slha_running_masses.push_back(TMass("DRbarMVP",
                                          to_valarray(DRBARSLHA(MVP))));
      slha_running_masses.push_back(TMass("DRbarMVWm",
                                          to_valarray(DRBARSLHA(MVWm))));
      slha_running_masses.push_back(TMass("DRbarMVZ",
                                          to_valarray(DRBARSLHA(MVZ))));

   }
}

template <class T>
void CSE6SSM_semianalytic_slha_values_writer::extract_slha_running_masses(
   const CSE6SSM_semianalytic_slha<T>& model
   )
{
   slha_running_masses.clear();
   slha_running_masses_inputs = model.get_input();
   slha_running_masses_problems = model.get_problems();
   slha_running_masses_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_running_masses.push_back(TMass("DRbarMGlu",
                                       to_valarray(DRBARSLHA(MGlu))));
   slha_running_masses.push_back(TMass("DRbarMChaP",
                                       to_valarray(DRBARSLHA(MChaP))));
   slha_running_masses.push_back(TMass("DRbarMVZp",
                                       to_valarray(DRBARSLHA(MVZp))));
   slha_running_masses.push_back(TMass("DRbarMSd",
                                       to_valarray(DRBARSLHA(MSd))));
   slha_running_masses.push_back(TMass("DRbarMSv",
                                       to_valarray(DRBARSLHA(MSv))));
   slha_running_masses.push_back(TMass("DRbarMSu",
                                       to_valarray(DRBARSLHA(MSu))));
   slha_running_masses.push_back(TMass("DRbarMSe",
                                       to_valarray(DRBARSLHA(MSe))));
   slha_running_masses.push_back(TMass("DRbarMSDX",
                                       to_valarray(DRBARSLHA(MSDX))));
   slha_running_masses.push_back(TMass("DRbarMhh",
                                       to_valarray(DRBARSLHA(Mhh))));
   slha_running_masses.push_back(TMass("DRbarMAh",
                                       to_valarray(DRBARSLHA(MAh))));
   slha_running_masses.push_back(TMass("DRbarMHpm",
                                       to_valarray(DRBARSLHA(MHpm))));
   slha_running_masses.push_back(TMass("DRbarMChi",
                                       to_valarray(DRBARSLHA(MChi))));
   slha_running_masses.push_back(TMass("DRbarMCha",
                                       to_valarray(DRBARSLHA(MCha))));
   slha_running_masses.push_back(TMass("DRbarMFDX",
                                       to_valarray(DRBARSLHA(MFDX))));
   slha_running_masses.push_back(TMass("DRbarMSHI0",
                                       to_valarray(DRBARSLHA(MSHI0))));
   slha_running_masses.push_back(TMass("DRbarMSHIPM",
                                       to_valarray(DRBARSLHA(MSHIPM))));
   slha_running_masses.push_back(TMass("DRbarMChaI",
                                       to_valarray(DRBARSLHA(MChaI))));
   slha_running_masses.push_back(TMass("DRbarMChiI",
                                       to_valarray(DRBARSLHA(MChiI))));
   slha_running_masses.push_back(TMass("DRbarMSHp0",
                                       to_valarray(DRBARSLHA(MSHp0))));
   slha_running_masses.push_back(TMass("DRbarMSHpp",
                                       to_valarray(DRBARSLHA(MSHpp))));
   slha_running_masses.push_back(TMass("DRbarMChiP",
                                       to_valarray(DRBARSLHA(MChiP))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_running_masses.push_back(TMass("DRbarMFd",
                                          to_valarray(DRBARSLHA(MFd))));
      slha_running_masses.push_back(TMass("DRbarMFe",
                                          to_valarray(DRBARSLHA(MFe))));
      slha_running_masses.push_back(TMass("DRbarMFu",
                                          to_valarray(DRBARSLHA(MFu))));
      slha_running_masses.push_back(TMass("DRbarMFv",
                                          to_valarray(DRBARSLHA(MFv))));
      slha_running_masses.push_back(TMass("DRbarMVG",
                                          to_valarray(DRBARSLHA(MVG))));
      slha_running_masses.push_back(TMass("DRbarMVP",
                                          to_valarray(DRBARSLHA(MVP))));
      slha_running_masses.push_back(TMass("DRbarMVWm",
                                          to_valarray(DRBARSLHA(MVWm))));
      slha_running_masses.push_back(TMass("DRbarMVZ",
                                          to_valarray(DRBARSLHA(MVZ))));

   }
}

template <class T>
void CSE6SSM_slha_values_writer::extract_slha_susy_pars(
   const CSE6SSM_slha<T>& model
   )
{
   slha_susy_pars.clear();
   slha_susy_pars_inputs = model.get_input();
   slha_susy_pars_problems = model.get_problems();

   // convert SLHA Yukawas to matrices
   Eigen::Matrix<double,3,3> Ye_matrix(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd_matrix(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yu_matrix(Eigen::Matrix<double,3,3>::Zero());

   for (std::size_t i = 0; i < 3; ++i) {
      Ye_matrix(i,i) = MODELPARAMETER(Ye_slha)(i);
      Yd_matrix(i,i) = MODELPARAMETER(Yd_slha)(i);
      Yu_matrix(i,i) = MODELPARAMETER(Yu_slha)(i);
   }

   slha_susy_pars.push_back(TParameter("Yd", 0, 3, 3, to_valarray(Yd_matrix)));
   slha_susy_pars.push_back(TParameter("hE", 0, 3, 2,
                                       to_valarray(MODELPARAMETER(hE))));
   slha_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(Ye_matrix)));
   slha_susy_pars.push_back(TParameter("SigmaL", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(SigmaL))));
   slha_susy_pars.push_back(TParameter("KappaPr", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(KappaPr))));
   slha_susy_pars.push_back(TParameter("Sigmax", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(Sigmax))));
   slha_susy_pars.push_back(TParameter("gD", 0, 3, 3,
                                       to_valarray(MODELPARAMETER(gD))));
   slha_susy_pars.push_back(TParameter("Kappa", 0, 3, 3,
                                       to_valarray(MODELPARAMETER(Kappa))));
   slha_susy_pars.push_back(TParameter("Lambda12", 0, 2, 2,
                                       to_valarray(MODELPARAMETER(Lambda12))));
   slha_susy_pars.push_back(TParameter("Lambdax", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(Lambdax))));
   slha_susy_pars.push_back(TParameter("fu", 0, 3, 2,
                                       to_valarray(MODELPARAMETER(fu))));
   slha_susy_pars.push_back(TParameter("fd", 0, 3, 2,
                                       to_valarray(MODELPARAMETER(fd))));
   slha_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(Yu_matrix)));
   slha_susy_pars.push_back(TParameter("MuPr", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MuPr))));
   slha_susy_pars.push_back(TParameter("MuPhi", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MuPhi))));
   slha_susy_pars.push_back(TParameter("XiF", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(XiF))));
   slha_susy_pars.push_back(
      TParameter("gY", 0, 1, 1,
                 to_valarray(MODELPARAMETER(g1) * 0.7745966692414834)));
   slha_susy_pars.push_back(TParameter("g2", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(g2))));
   slha_susy_pars.push_back(TParameter("g3", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(g3))));
   slha_susy_pars.push_back(TParameter("g1p", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(g1p))));
   slha_susy_pars.push_back(TParameter("vd", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vd))));
   slha_susy_pars.push_back(TParameter("vu", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vu))));
   slha_susy_pars.push_back(TParameter("vs", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vs))));
   slha_susy_pars.push_back(TParameter("vsb", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vsb))));
   slha_susy_pars.push_back(TParameter("vphi", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vphi))));
   slha_susy_pars.push_back(TParameter("QS", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(QS))));
   slha_susy_pars.push_back(
      TParameter("Beta", 0, 1, 1,
                 to_valarray(ArcTan((MODELPARAMETER(vu))
                                    / (MODELPARAMETER(vd))))));
}

template <class T>
void CSE6SSM_semianalytic_slha_values_writer::extract_slha_susy_pars(
   const CSE6SSM_semianalytic_slha<T>& model
   )
{
   slha_susy_pars.clear();
   slha_susy_pars_inputs = model.get_input();
   slha_susy_pars_problems = model.get_problems();
   slha_susy_pars_m0Sqr = model.get_ewsb_output_parameter(0);

   // convert SLHA Yukawas to matrices
   Eigen::Matrix<double,3,3> Ye_matrix(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd_matrix(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yu_matrix(Eigen::Matrix<double,3,3>::Zero());

   for (std::size_t i = 0; i < 3; ++i) {
      Ye_matrix(i,i) = MODELPARAMETER(Ye_slha)(i);
      Yd_matrix(i,i) = MODELPARAMETER(Yd_slha)(i);
      Yu_matrix(i,i) = MODELPARAMETER(Yu_slha)(i);
   }

   slha_susy_pars.push_back(TParameter("Yd", 0, 3, 3, to_valarray(Yd_matrix)));
   slha_susy_pars.push_back(TParameter("hE", 0, 3, 2,
                                       to_valarray(MODELPARAMETER(hE))));
   slha_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(Ye_matrix)));
   slha_susy_pars.push_back(TParameter("SigmaL", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(SigmaL))));
   slha_susy_pars.push_back(TParameter("KappaPr", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(KappaPr))));
   slha_susy_pars.push_back(TParameter("Sigmax", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(Sigmax))));
   slha_susy_pars.push_back(TParameter("gD", 0, 3, 3,
                                       to_valarray(MODELPARAMETER(gD))));
   slha_susy_pars.push_back(TParameter("Kappa", 0, 3, 3,
                                       to_valarray(MODELPARAMETER(Kappa))));
   slha_susy_pars.push_back(TParameter("Lambda12", 0, 2, 2,
                                       to_valarray(MODELPARAMETER(Lambda12))));
   slha_susy_pars.push_back(TParameter("Lambdax", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(Lambdax))));
   slha_susy_pars.push_back(TParameter("fu", 0, 3, 2,
                                       to_valarray(MODELPARAMETER(fu))));
   slha_susy_pars.push_back(TParameter("fd", 0, 3, 2,
                                       to_valarray(MODELPARAMETER(fd))));
   slha_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(Yu_matrix)));
   slha_susy_pars.push_back(TParameter("MuPr", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MuPr))));
   slha_susy_pars.push_back(TParameter("MuPhi", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MuPhi))));
   slha_susy_pars.push_back(TParameter("XiF", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(XiF))));
   slha_susy_pars.push_back(
      TParameter("gY", 0, 1, 1,
                 to_valarray(MODELPARAMETER(g1) * 0.7745966692414834)));
   slha_susy_pars.push_back(TParameter("g2", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(g2))));
   slha_susy_pars.push_back(TParameter("g3", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(g3))));
   slha_susy_pars.push_back(TParameter("g1p", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(g1p))));
   slha_susy_pars.push_back(TParameter("vd", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vd))));
   slha_susy_pars.push_back(TParameter("vu", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vu))));
   slha_susy_pars.push_back(TParameter("vs", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vs))));
   slha_susy_pars.push_back(TParameter("vsb", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vsb))));
   slha_susy_pars.push_back(TParameter("vphi", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(vphi))));
   slha_susy_pars.push_back(TParameter("QS", 0, 1, 1,
                                       to_valarray(MODELPARAMETER(QS))));
   slha_susy_pars.push_back(
      TParameter("Beta", 0, 1, 1,
                 to_valarray(ArcTan((MODELPARAMETER(vu))
                                    / (MODELPARAMETER(vd))))));
}

template <class T>
void CSE6SSM_slha_values_writer::extract_slha_soft_pars(
   const CSE6SSM_slha<T>& model
   )
{
   slha_soft_pars.clear();
   slha_soft_pars_inputs = model.get_input();
   slha_soft_pars_problems = model.get_problems();

   slha_soft_pars.push_back(TParameter("TYd", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TYd_slha))));
   slha_soft_pars.push_back(TParameter("ThE", 1, 3, 2,
                                       to_valarray(MODELPARAMETER(ThE))));
   slha_soft_pars.push_back(TParameter("TYe", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TYe_slha))));
   slha_soft_pars.push_back(TParameter("TSigmaL", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TSigmaL))));
   slha_soft_pars.push_back(TParameter("TKappaPr", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TKappaPr))));
   slha_soft_pars.push_back(TParameter("TSigmax", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TSigmax))));
   slha_soft_pars.push_back(TParameter("TgD", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TgD))));
   slha_soft_pars.push_back(TParameter("TKappa", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TKappa))));
   slha_soft_pars.push_back(TParameter("TLambda12", 1, 2, 2,
                                       to_valarray(MODELPARAMETER(TLambda12))));
   slha_soft_pars.push_back(TParameter("TLambdax", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TLambdax))));
   slha_soft_pars.push_back(TParameter("Tfu", 1, 3, 2,
                                       to_valarray(MODELPARAMETER(Tfu))));
   slha_soft_pars.push_back(TParameter("Tfd", 1, 3, 2,
                                       to_valarray(MODELPARAMETER(Tfd))));
   slha_soft_pars.push_back(TParameter("TYu", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TYu_slha))));
   slha_soft_pars.push_back(TParameter("BMuPr", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(BMuPr))));
   slha_soft_pars.push_back(TParameter("BMuPhi", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(BMuPhi))));
   slha_soft_pars.push_back(TParameter("LXiF", 3, 1, 1,
                                       to_valarray(MODELPARAMETER(LXiF))));
   slha_soft_pars.push_back(TParameter("mq2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mq2_slha))));
   slha_soft_pars.push_back(TParameter("ml2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(ml2_slha))));
   slha_soft_pars.push_back(TParameter("mHd2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHd2))));
   slha_soft_pars.push_back(TParameter("mHu2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHu2))));
   slha_soft_pars.push_back(TParameter("md2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(md2_slha))));
   slha_soft_pars.push_back(TParameter("mu2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mu2_slha))));
   slha_soft_pars.push_back(TParameter("me2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(me2_slha))));
   slha_soft_pars.push_back(TParameter("ms2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(ms2))));
   slha_soft_pars.push_back(TParameter("msbar2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(msbar2))));
   slha_soft_pars.push_back(TParameter("mH1I2", 2, 2, 2,
                                       to_valarray(MODELPARAMETER(mH1I2))));
   slha_soft_pars.push_back(TParameter("mH2I2", 2, 2, 2,
                                       to_valarray(MODELPARAMETER(mH2I2))));
   slha_soft_pars.push_back(TParameter("mSI2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mSI2))));
   slha_soft_pars.push_back(TParameter("mDx2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mDx2))));
   slha_soft_pars.push_back(TParameter("mDxbar2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mDxbar2))));
   slha_soft_pars.push_back(TParameter("mHp2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHp2))));
   slha_soft_pars.push_back(TParameter("mHpbar2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHpbar2))));
   slha_soft_pars.push_back(TParameter("mphi2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mphi2))));
   slha_soft_pars.push_back(TParameter("MassB", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassB))));
   slha_soft_pars.push_back(TParameter("MassWB", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassWB))));
   slha_soft_pars.push_back(TParameter("MassG", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassG))));
   slha_soft_pars.push_back(TParameter("MassBp", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassBp))));
   slha_soft_pars.push_back(
      TParameter("Re(PhaseGlu)", 0, 1, 1,
                 to_valarray(Re(MODELPARAMETER(PhaseGlu)))));
   slha_soft_pars.push_back(
      TParameter("Im(PhaseGlu)", 0, 1, 1,
                 to_valarray(Im(MODELPARAMETER(PhaseGlu)))));
}

template <class T>
void CSE6SSM_semianalytic_slha_values_writer::extract_slha_soft_pars(
   const CSE6SSM_semianalytic_slha<T>& model
   )
{
   slha_soft_pars.clear();
   slha_soft_pars_inputs = model.get_input();
   slha_soft_pars_problems = model.get_problems();
   slha_soft_pars_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_soft_pars.push_back(TParameter("TYd", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TYd_slha))));
   slha_soft_pars.push_back(TParameter("ThE", 1, 3, 2,
                                       to_valarray(MODELPARAMETER(ThE))));
   slha_soft_pars.push_back(TParameter("TYe", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TYe_slha))));
   slha_soft_pars.push_back(TParameter("TSigmaL", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TSigmaL))));
   slha_soft_pars.push_back(TParameter("TKappaPr", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TKappaPr))));
   slha_soft_pars.push_back(TParameter("TSigmax", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TSigmax))));
   slha_soft_pars.push_back(TParameter("TgD", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TgD))));
   slha_soft_pars.push_back(TParameter("TKappa", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TKappa))));
   slha_soft_pars.push_back(TParameter("TLambda12", 1, 2, 2,
                                       to_valarray(MODELPARAMETER(TLambda12))));
   slha_soft_pars.push_back(TParameter("TLambdax", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(TLambdax))));
   slha_soft_pars.push_back(TParameter("Tfu", 1, 3, 2,
                                       to_valarray(MODELPARAMETER(Tfu))));
   slha_soft_pars.push_back(TParameter("Tfd", 1, 3, 2,
                                       to_valarray(MODELPARAMETER(Tfd))));
   slha_soft_pars.push_back(TParameter("TYu", 1, 3, 3,
                                       to_valarray(MODELPARAMETER(TYu_slha))));
   slha_soft_pars.push_back(TParameter("BMuPr", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(BMuPr))));
   slha_soft_pars.push_back(TParameter("BMuPhi", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(BMuPhi))));
   slha_soft_pars.push_back(TParameter("LXiF", 3, 1, 1,
                                       to_valarray(MODELPARAMETER(LXiF))));
   slha_soft_pars.push_back(TParameter("mq2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mq2_slha))));
   slha_soft_pars.push_back(TParameter("ml2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(ml2_slha))));
   slha_soft_pars.push_back(TParameter("mHd2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHd2))));
   slha_soft_pars.push_back(TParameter("mHu2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHu2))));
   slha_soft_pars.push_back(TParameter("md2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(md2_slha))));
   slha_soft_pars.push_back(TParameter("mu2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mu2_slha))));
   slha_soft_pars.push_back(TParameter("me2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(me2_slha))));
   slha_soft_pars.push_back(TParameter("ms2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(ms2))));
   slha_soft_pars.push_back(TParameter("msbar2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(msbar2))));
   slha_soft_pars.push_back(TParameter("mH1I2", 2, 2, 2,
                                       to_valarray(MODELPARAMETER(mH1I2))));
   slha_soft_pars.push_back(TParameter("mH2I2", 2, 2, 2,
                                       to_valarray(MODELPARAMETER(mH2I2))));
   slha_soft_pars.push_back(TParameter("mSI2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mSI2))));
   slha_soft_pars.push_back(TParameter("mDx2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mDx2))));
   slha_soft_pars.push_back(TParameter("mDxbar2", 2, 3, 3,
                                       to_valarray(MODELPARAMETER(mDxbar2))));
   slha_soft_pars.push_back(TParameter("mHp2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHp2))));
   slha_soft_pars.push_back(TParameter("mHpbar2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mHpbar2))));
   slha_soft_pars.push_back(TParameter("mphi2", 2, 1, 1,
                                       to_valarray(MODELPARAMETER(mphi2))));
   slha_soft_pars.push_back(TParameter("MassB", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassB))));
   slha_soft_pars.push_back(TParameter("MassWB", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassWB))));
   slha_soft_pars.push_back(TParameter("MassG", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassG))));
   slha_soft_pars.push_back(TParameter("MassBp", 1, 1, 1,
                                       to_valarray(MODELPARAMETER(MassBp))));
   slha_soft_pars.push_back(
      TParameter("Re(PhaseGlu)", 0, 1, 1,
                 to_valarray(Re(MODELPARAMETER(PhaseGlu)))));
   slha_soft_pars.push_back(
      TParameter("Im(PhaseGlu)", 0, 1, 1,
                 to_valarray(Im(MODELPARAMETER(PhaseGlu)))));
}

template <class T>
void CSE6SSM_slha_values_writer::extract_slha_pole_mixings(
   const CSE6SSM_slha<T>& model
   )
{
   slha_pole_mixings.clear();
   slha_pole_mixings_inputs = model.get_input();
   slha_pole_mixings_problems = model.get_problems();

   slha_pole_mixings.push_back(TMixing("ZD", 6, true,
                                       to_valarray(PHYSICALSLHA(ZD))));
   slha_pole_mixings.push_back(TMixing("ZV", 3, true,
                                       to_valarray(PHYSICALSLHA(ZV))));
   slha_pole_mixings.push_back(TMixing("ZU", 6, true,
                                       to_valarray(PHYSICALSLHA(ZU))));
   slha_pole_mixings.push_back(TMixing("ZE", 6, true,
                                       to_valarray(PHYSICALSLHA(ZE))));
   slha_pole_mixings.push_back(TMixing("ZDX", 6, true,
                                       to_valarray(PHYSICALSLHA(ZDX))));
   slha_pole_mixings.push_back(TMixing("ZH", 5, true,
                                       to_valarray(PHYSICALSLHA(ZH))));
   slha_pole_mixings.push_back(TMixing("ZA", 5, true,
                                       to_valarray(PHYSICALSLHA(ZA))));
   slha_pole_mixings.push_back(TMixing("ZP", 2, true,
                                       to_valarray(PHYSICALSLHA(ZP))));
   slha_pole_mixings.push_back(TMixing("ZN", 8, false,
                                       to_valarray(PHYSICALSLHA(ZN))));
   slha_pole_mixings.push_back(TMixing("UM", 2, false,
                                       to_valarray(PHYSICALSLHA(UM))));
   slha_pole_mixings.push_back(TMixing("UP", 2, false,
                                       to_valarray(PHYSICALSLHA(UP))));
   slha_pole_mixings.push_back(TMixing("ZEL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZEL))));
   slha_pole_mixings.push_back(TMixing("ZER", 3, false,
                                       to_valarray(PHYSICALSLHA(ZER))));
   slha_pole_mixings.push_back(TMixing("ZDL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDL))));
   slha_pole_mixings.push_back(TMixing("ZDR", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDR))));
   slha_pole_mixings.push_back(TMixing("ZUL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZUL))));
   slha_pole_mixings.push_back(TMixing("ZUR", 3, false,
                                       to_valarray(PHYSICALSLHA(ZUR))));
   slha_pole_mixings.push_back(TMixing("ZDXL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDXL))));
   slha_pole_mixings.push_back(TMixing("ZDXR", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDXR))));
   slha_pole_mixings.push_back(TMixing("UHI0", 7, true,
                                       to_valarray(PHYSICALSLHA(UHI0))));
   slha_pole_mixings.push_back(TMixing("UHIPM", 4, true,
                                       to_valarray(PHYSICALSLHA(UHIPM))));
   slha_pole_mixings.push_back(TMixing("ZMI", 2, false,
                                       to_valarray(PHYSICALSLHA(ZMI))));
   slha_pole_mixings.push_back(TMixing("ZPI", 2, false,
                                       to_valarray(PHYSICALSLHA(ZPI))));
   slha_pole_mixings.push_back(TMixing("ZNI", 7, false,
                                       to_valarray(PHYSICALSLHA(ZNI))));
   slha_pole_mixings.push_back(TMixing("UHp0", 2, true,
                                       to_valarray(PHYSICALSLHA(UHp0))));
   slha_pole_mixings.push_back(TMixing("UHpp", 2, true,
                                       to_valarray(PHYSICALSLHA(UHpp))));
   slha_pole_mixings.push_back(TMixing("ZNp", 2, false,
                                       to_valarray(PHYSICALSLHA(ZNp))));
}

template <class T>
void CSE6SSM_semianalytic_slha_values_writer::extract_slha_pole_mixings(
   const CSE6SSM_semianalytic_slha<T>& model
   )
{
   slha_pole_mixings.clear();
   slha_pole_mixings_inputs = model.get_input();
   slha_pole_mixings_problems = model.get_problems();
   slha_pole_mixings_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_pole_mixings.push_back(TMixing("ZD", 6, true,
                                       to_valarray(PHYSICALSLHA(ZD))));
   slha_pole_mixings.push_back(TMixing("ZV", 3, true,
                                       to_valarray(PHYSICALSLHA(ZV))));
   slha_pole_mixings.push_back(TMixing("ZU", 6, true,
                                       to_valarray(PHYSICALSLHA(ZU))));
   slha_pole_mixings.push_back(TMixing("ZE", 6, true,
                                       to_valarray(PHYSICALSLHA(ZE))));
   slha_pole_mixings.push_back(TMixing("ZDX", 6, true,
                                       to_valarray(PHYSICALSLHA(ZDX))));
   slha_pole_mixings.push_back(TMixing("ZH", 5, true,
                                       to_valarray(PHYSICALSLHA(ZH))));
   slha_pole_mixings.push_back(TMixing("ZA", 5, true,
                                       to_valarray(PHYSICALSLHA(ZA))));
   slha_pole_mixings.push_back(TMixing("ZP", 2, true,
                                       to_valarray(PHYSICALSLHA(ZP))));
   slha_pole_mixings.push_back(TMixing("ZN", 8, false,
                                       to_valarray(PHYSICALSLHA(ZN))));
   slha_pole_mixings.push_back(TMixing("UM", 2, false,
                                       to_valarray(PHYSICALSLHA(UM))));
   slha_pole_mixings.push_back(TMixing("UP", 2, false,
                                       to_valarray(PHYSICALSLHA(UP))));
   slha_pole_mixings.push_back(TMixing("ZEL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZEL))));
   slha_pole_mixings.push_back(TMixing("ZER", 3, false,
                                       to_valarray(PHYSICALSLHA(ZER))));
   slha_pole_mixings.push_back(TMixing("ZDL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDL))));
   slha_pole_mixings.push_back(TMixing("ZDR", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDR))));
   slha_pole_mixings.push_back(TMixing("ZUL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZUL))));
   slha_pole_mixings.push_back(TMixing("ZUR", 3, false,
                                       to_valarray(PHYSICALSLHA(ZUR))));
   slha_pole_mixings.push_back(TMixing("ZDXL", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDXL))));
   slha_pole_mixings.push_back(TMixing("ZDXR", 3, false,
                                       to_valarray(PHYSICALSLHA(ZDXR))));
   slha_pole_mixings.push_back(TMixing("UHI0", 7, true,
                                       to_valarray(PHYSICALSLHA(UHI0))));
   slha_pole_mixings.push_back(TMixing("UHIPM", 4, true,
                                       to_valarray(PHYSICALSLHA(UHIPM))));
   slha_pole_mixings.push_back(TMixing("ZMI", 2, false,
                                       to_valarray(PHYSICALSLHA(ZMI))));
   slha_pole_mixings.push_back(TMixing("ZPI", 2, false,
                                       to_valarray(PHYSICALSLHA(ZPI))));
   slha_pole_mixings.push_back(TMixing("ZNI", 7, false,
                                       to_valarray(PHYSICALSLHA(ZNI))));
   slha_pole_mixings.push_back(TMixing("UHp0", 2, true,
                                       to_valarray(PHYSICALSLHA(UHp0))));
   slha_pole_mixings.push_back(TMixing("UHpp", 2, true,
                                       to_valarray(PHYSICALSLHA(UHpp))));
   slha_pole_mixings.push_back(TMixing("ZNp", 2, false,
                                       to_valarray(PHYSICALSLHA(ZNp))));
}

template <class T>
void CSE6SSM_slha_values_writer::extract_slha_running_mixings(
   const CSE6SSM_slha<T>& model
   )
{
   slha_running_mixings.clear();
   slha_running_mixings_inputs = model.get_input();
   slha_running_mixings_problems = model.get_problems();

   slha_running_mixings.push_back(TMixing("DRbarZD", 6, true,
                                          to_valarray(DRBARSLHA(ZD))));
   slha_running_mixings.push_back(TMixing("DRbarZV", 3, true,
                                          to_valarray(DRBARSLHA(ZV))));
   slha_running_mixings.push_back(TMixing("DRbarZU", 6, true,
                                          to_valarray(DRBARSLHA(ZU))));
   slha_running_mixings.push_back(TMixing("DRbarZE", 6, true,
                                          to_valarray(DRBARSLHA(ZE))));
   slha_running_mixings.push_back(TMixing("DRbarZDX", 6, true,
                                          to_valarray(DRBARSLHA(ZDX))));
   slha_running_mixings.push_back(TMixing("DRbarZH", 5, true,
                                          to_valarray(DRBARSLHA(ZH))));
   slha_running_mixings.push_back(TMixing("DRbarZA", 5, true,
                                          to_valarray(DRBARSLHA(ZA))));
   slha_running_mixings.push_back(TMixing("DRbarZP", 2, true,
                                          to_valarray(DRBARSLHA(ZP))));
   slha_running_mixings.push_back(TMixing("DRbarZN", 8, false,
                                          to_valarray(DRBARSLHA(ZN))));
   slha_running_mixings.push_back(TMixing("DRbarUM", 2, false,
                                          to_valarray(DRBARSLHA(UM))));
   slha_running_mixings.push_back(TMixing("DRbarUP", 2, false,
                                          to_valarray(DRBARSLHA(UP))));
   slha_running_mixings.push_back(TMixing("DRbarZEL", 3, false,
                                          to_valarray(DRBARSLHA(ZEL))));
   slha_running_mixings.push_back(TMixing("DRbarZER", 3, false,
                                          to_valarray(DRBARSLHA(ZER))));
   slha_running_mixings.push_back(TMixing("DRbarZDL", 3, false,
                                          to_valarray(DRBARSLHA(ZDL))));
   slha_running_mixings.push_back(TMixing("DRbarZDR", 3, false,
                                          to_valarray(DRBARSLHA(ZDR))));
   slha_running_mixings.push_back(TMixing("DRbarZUL", 3, false,
                                          to_valarray(DRBARSLHA(ZUL))));
   slha_running_mixings.push_back(TMixing("DRbarZUR", 3, false,
                                          to_valarray(DRBARSLHA(ZUR))));
   slha_running_mixings.push_back(TMixing("DRbarZDXL", 3, false,
                                          to_valarray(DRBARSLHA(ZDXL))));
   slha_running_mixings.push_back(TMixing("DRbarZDXR", 3, false,
                                          to_valarray(DRBARSLHA(ZDXR))));
   slha_running_mixings.push_back(TMixing("DRbarUHI0", 7, true,
                                          to_valarray(DRBARSLHA(UHI0))));
   slha_running_mixings.push_back(TMixing("DRbarUHIPM", 4, true,
                                          to_valarray(DRBARSLHA(UHIPM))));
   slha_running_mixings.push_back(TMixing("DRbarZMI", 2, false,
                                          to_valarray(DRBARSLHA(ZMI))));
   slha_running_mixings.push_back(TMixing("DRbarZPI", 2, false,
                                          to_valarray(DRBARSLHA(ZPI))));
   slha_running_mixings.push_back(TMixing("DRbarZNI", 7, false,
                                          to_valarray(DRBARSLHA(ZNI))));
   slha_running_mixings.push_back(TMixing("DRbarUHp0", 2, true,
                                          to_valarray(DRBARSLHA(UHp0))));
   slha_running_mixings.push_back(TMixing("DRbarUHpp", 2, true,
                                          to_valarray(DRBARSLHA(UHpp))));
   slha_running_mixings.push_back(TMixing("DRbarZNp", 2, false,
                                          to_valarray(DRBARSLHA(ZNp))));
}

template <class T>
void CSE6SSM_semianalytic_slha_values_writer::extract_slha_running_mixings(
   const CSE6SSM_semianalytic_slha<T>& model
   )
{
   slha_running_mixings.clear();
   slha_running_mixings_inputs = model.get_input();
   slha_running_mixings_problems = model.get_problems();
   slha_running_mixings_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_running_mixings.push_back(TMixing("DRbarZD", 6, true,
                                          to_valarray(DRBARSLHA(ZD))));
   slha_running_mixings.push_back(TMixing("DRbarZV", 3, true,
                                          to_valarray(DRBARSLHA(ZV))));
   slha_running_mixings.push_back(TMixing("DRbarZU", 6, true,
                                          to_valarray(DRBARSLHA(ZU))));
   slha_running_mixings.push_back(TMixing("DRbarZE", 6, true,
                                          to_valarray(DRBARSLHA(ZE))));
   slha_running_mixings.push_back(TMixing("DRbarZDX", 6, true,
                                          to_valarray(DRBARSLHA(ZDX))));
   slha_running_mixings.push_back(TMixing("DRbarZH", 5, true,
                                          to_valarray(DRBARSLHA(ZH))));
   slha_running_mixings.push_back(TMixing("DRbarZA", 5, true,
                                          to_valarray(DRBARSLHA(ZA))));
   slha_running_mixings.push_back(TMixing("DRbarZP", 2, true,
                                          to_valarray(DRBARSLHA(ZP))));
   slha_running_mixings.push_back(TMixing("DRbarZN", 8, false,
                                          to_valarray(DRBARSLHA(ZN))));
   slha_running_mixings.push_back(TMixing("DRbarUM", 2, false,
                                          to_valarray(DRBARSLHA(UM))));
   slha_running_mixings.push_back(TMixing("DRbarUP", 2, false,
                                          to_valarray(DRBARSLHA(UP))));
   slha_running_mixings.push_back(TMixing("DRbarZEL", 3, false,
                                          to_valarray(DRBARSLHA(ZEL))));
   slha_running_mixings.push_back(TMixing("DRbarZER", 3, false,
                                          to_valarray(DRBARSLHA(ZER))));
   slha_running_mixings.push_back(TMixing("DRbarZDL", 3, false,
                                          to_valarray(DRBARSLHA(ZDL))));
   slha_running_mixings.push_back(TMixing("DRbarZDR", 3, false,
                                          to_valarray(DRBARSLHA(ZDR))));
   slha_running_mixings.push_back(TMixing("DRbarZUL", 3, false,
                                          to_valarray(DRBARSLHA(ZUL))));
   slha_running_mixings.push_back(TMixing("DRbarZUR", 3, false,
                                          to_valarray(DRBARSLHA(ZUR))));
   slha_running_mixings.push_back(TMixing("DRbarZDXL", 3, false,
                                          to_valarray(DRBARSLHA(ZDXL))));
   slha_running_mixings.push_back(TMixing("DRbarZDXR", 3, false,
                                          to_valarray(DRBARSLHA(ZDXR))));
   slha_running_mixings.push_back(TMixing("DRbarUHI0", 7, true,
                                          to_valarray(DRBARSLHA(UHI0))));
   slha_running_mixings.push_back(TMixing("DRbarUHIPM", 4, true,
                                          to_valarray(DRBARSLHA(UHIPM))));
   slha_running_mixings.push_back(TMixing("DRbarZMI", 2, false,
                                          to_valarray(DRBARSLHA(ZMI))));
   slha_running_mixings.push_back(TMixing("DRbarZPI", 2, false,
                                          to_valarray(DRBARSLHA(ZPI))));
   slha_running_mixings.push_back(TMixing("DRbarZNI", 7, false,
                                          to_valarray(DRBARSLHA(ZNI))));
   slha_running_mixings.push_back(TMixing("DRbarUHp0", 2, true,
                                          to_valarray(DRBARSLHA(UHp0))));
   slha_running_mixings.push_back(TMixing("DRbarUHpp", 2, true,
                                          to_valarray(DRBARSLHA(UHpp))));
   slha_running_mixings.push_back(TMixing("DRbarZNp", 2, false,
                                          to_valarray(DRBARSLHA(ZNp))));
}

template <class T>
void CSE6SSM_semianalytic_coefficients_writer::extract_coefficients(
   const CSE6SSM_semianalytic<T>& model
   )
{
   coefficients.clear();
   coefficients_inputs = model.get_input();
   coefficients_problems = model.get_problems();
   coefficients_m0Sqr = model.get_ewsb_output_parameter(0);

   coefficients.push_back(
      TCoefficient("TYd_Azero_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TYd, Azero))));
   coefficients.push_back(
      TCoefficient("TYd_m12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TYd, m12))));
   coefficients.push_back(
      TCoefficient("ThE_Azero_coeff", 0, 3, 2,
                   to_valarray(COEFFICIENT(ThE, Azero))));
   coefficients.push_back(
      TCoefficient("ThE_m12_coeff", 0, 3, 2,
                   to_valarray(COEFFICIENT(ThE, m12))));
   coefficients.push_back(
      TCoefficient("TYe_Azero_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TYe, Azero))));
   coefficients.push_back(
      TCoefficient("TYe_m12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TYe, m12))));
   coefficients.push_back(
      TCoefficient("TSigmaL_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TSigmaL, Azero))));
   coefficients.push_back(
      TCoefficient("TSigmaL_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TSigmaL, m12))));
   coefficients.push_back(
      TCoefficient("TKappaPr_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TKappaPr, Azero))));
   coefficients.push_back(
      TCoefficient("TKappaPr_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TKappaPr, m12))));
   coefficients.push_back(
      TCoefficient("TSigmax_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TSigmax, Azero))));
   coefficients.push_back(
      TCoefficient("TSigmax_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TSigmax, m12))));
   coefficients.push_back(
      TCoefficient("TgD_Azero_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TgD, Azero))));
   coefficients.push_back(
      TCoefficient("TgD_m12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TgD, m12))));
   coefficients.push_back(
      TCoefficient("TKappa_Azero_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TKappa, Azero))));
   coefficients.push_back(
      TCoefficient("TKappa_m12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TKappa, m12))));
   coefficients.push_back(
      TCoefficient("TLambda12_Azero_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(TLambda12, Azero))));
   coefficients.push_back(
      TCoefficient("TLambda12_m12_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(TLambda12, m12))));
   coefficients.push_back(
      TCoefficient("TLambdax_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TLambdax, Azero))));
   coefficients.push_back(
      TCoefficient("TLambdax_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(TLambdax, m12))));
   coefficients.push_back(
      TCoefficient("Tfu_Azero_coeff", 0, 3, 2,
                   to_valarray(COEFFICIENT(Tfu, Azero))));
   coefficients.push_back(
      TCoefficient("Tfu_m12_coeff", 0, 3, 2,
                   to_valarray(COEFFICIENT(Tfu, m12))));
   coefficients.push_back(
      TCoefficient("Tfd_Azero_coeff", 0, 3, 2,
                   to_valarray(COEFFICIENT(Tfd, Azero))));
   coefficients.push_back(
      TCoefficient("Tfd_m12_coeff", 0, 3, 2,
                   to_valarray(COEFFICIENT(Tfd, m12))));
   coefficients.push_back(
      TCoefficient("TYu_Azero_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TYu, Azero))));
   coefficients.push_back(
      TCoefficient("TYu_m12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(TYu, m12))));
   coefficients.push_back(
      TCoefficient("BMuPr_BMuPr_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(BMuPr, BMuPr))));
   coefficients.push_back(
      TCoefficient("BMuPr_BMuPhi_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(BMuPr, BMuPhi))));
   coefficients.push_back(
      TCoefficient("BMuPr_Azero_coeff", 1, 1, 1,
                   to_valarray(COEFFICIENT(BMuPr, Azero))));
   coefficients.push_back(
      TCoefficient("BMuPr_m12_coeff", 1, 1, 1,
                   to_valarray(COEFFICIENT(BMuPr, m12))));
   coefficients.push_back(
      TCoefficient("BMuPhi_BMuPr_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(BMuPhi, BMuPr))));
   coefficients.push_back(
      TCoefficient("BMuPhi_BMuPhi_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(BMuPhi, BMuPhi))));
   coefficients.push_back(
      TCoefficient("BMuPhi_Azero_coeff", 1, 1, 1,
                   to_valarray(COEFFICIENT(BMuPhi, Azero))));
   coefficients.push_back(
      TCoefficient("BMuPhi_m12_coeff", 1, 1, 1,
                   to_valarray(COEFFICIENT(BMuPhi, m12))));
   coefficients.push_back(
      TCoefficient("mq2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mq2, m02))));
   coefficients.push_back(
      TCoefficient("mq2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mq2, m122))));
   coefficients.push_back(
      TCoefficient("mq2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mq2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mq2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mq2, Azero2))));
   coefficients.push_back(
      TCoefficient("ml2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(ml2, m02))));
   coefficients.push_back(
      TCoefficient("ml2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(ml2, m122))));
   coefficients.push_back(
      TCoefficient("ml2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(ml2, Azerom12))));
   coefficients.push_back(
      TCoefficient("ml2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(ml2, Azero2))));
   coefficients.push_back(
      TCoefficient("mHd2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHd2, m02))));
   coefficients.push_back(
      TCoefficient("mHd2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHd2, m122))));
   coefficients.push_back(
      TCoefficient("mHd2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHd2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mHd2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHd2, Azero2))));
   coefficients.push_back(
      TCoefficient("mHu2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHu2, m02))));
   coefficients.push_back(
      TCoefficient("mHu2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHu2, m122))));
   coefficients.push_back(
      TCoefficient("mHu2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHu2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mHu2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHu2, Azero2))));
   coefficients.push_back(
      TCoefficient("md2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(md2, m02))));
   coefficients.push_back(
      TCoefficient("md2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(md2, m122))));
   coefficients.push_back(
      TCoefficient("md2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(md2, Azerom12))));
   coefficients.push_back(
      TCoefficient("md2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(md2, Azero2))));
   coefficients.push_back(
      TCoefficient("mu2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mu2, m02))));
   coefficients.push_back(
      TCoefficient("mu2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mu2, m122))));
   coefficients.push_back(
      TCoefficient("mu2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mu2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mu2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mu2, Azero2))));
   coefficients.push_back(
      TCoefficient("me2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(me2, m02))));
   coefficients.push_back(
      TCoefficient("me2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(me2, m122))));
   coefficients.push_back(
      TCoefficient("me2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(me2, Azerom12))));
   coefficients.push_back(
      TCoefficient("me2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(me2, Azero2))));
   coefficients.push_back(
      TCoefficient("ms2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(ms2, m02))));
   coefficients.push_back(
      TCoefficient("ms2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(ms2, m122))));
   coefficients.push_back(
      TCoefficient("ms2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(ms2, Azerom12))));
   coefficients.push_back(
      TCoefficient("ms2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(ms2, Azero2))));
   coefficients.push_back(
      TCoefficient("msbar2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(msbar2, m02))));
   coefficients.push_back(
      TCoefficient("msbar2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(msbar2, m122))));
   coefficients.push_back(
      TCoefficient("msbar2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(msbar2, Azerom12))));
   coefficients.push_back(
      TCoefficient("msbar2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(msbar2, Azero2))));
   coefficients.push_back(
      TCoefficient("mH1I2_m02_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH1I2, m02))));
   coefficients.push_back(
      TCoefficient("mH1I2_m122_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH1I2, m122))));
   coefficients.push_back(
      TCoefficient("mH1I2_Azerom12_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH1I2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mH1I2_Azero2_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH1I2, Azero2))));
   coefficients.push_back(
      TCoefficient("mH2I2_m02_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH2I2, m02))));
   coefficients.push_back(
      TCoefficient("mH2I2_m122_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH2I2, m122))));
   coefficients.push_back(
      TCoefficient("mH2I2_Azerom12_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH2I2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mH2I2_Azero2_coeff", 0, 2, 2,
                   to_valarray(COEFFICIENT(mH2I2, Azero2))));
   coefficients.push_back(
      TCoefficient("mSI2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mSI2, m02))));
   coefficients.push_back(
      TCoefficient("mSI2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mSI2, m122))));
   coefficients.push_back(
      TCoefficient("mSI2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mSI2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mSI2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mSI2, Azero2))));
   coefficients.push_back(
      TCoefficient("mDx2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDx2, m02))));
   coefficients.push_back(
      TCoefficient("mDx2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDx2, m122))));
   coefficients.push_back(
      TCoefficient("mDx2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDx2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mDx2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDx2, Azero2))));
   coefficients.push_back(
      TCoefficient("mDxbar2_m02_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDxbar2, m02))));
   coefficients.push_back(
      TCoefficient("mDxbar2_m122_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDxbar2, m122))));
   coefficients.push_back(
      TCoefficient("mDxbar2_Azerom12_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDxbar2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mDxbar2_Azero2_coeff", 0, 3, 3,
                   to_valarray(COEFFICIENT(mDxbar2, Azero2))));
   coefficients.push_back(
      TCoefficient("mHp2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHp2, m02))));
   coefficients.push_back(
      TCoefficient("mHp2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHp2, m122))));
   coefficients.push_back(
      TCoefficient("mHp2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHp2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mHp2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHp2, Azero2))));
   coefficients.push_back(
      TCoefficient("mHpbar2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHpbar2, m02))));
   coefficients.push_back(
      TCoefficient("mHpbar2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHpbar2, m122))));
   coefficients.push_back(
      TCoefficient("mHpbar2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHpbar2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mHpbar2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mHpbar2, Azero2))));
   coefficients.push_back(
      TCoefficient("mphi2_m02_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mphi2, m02))));
   coefficients.push_back(
      TCoefficient("mphi2_m122_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mphi2, m122))));
   coefficients.push_back(
      TCoefficient("mphi2_Azerom12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mphi2, Azerom12))));
   coefficients.push_back(
      TCoefficient("mphi2_Azero2_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(mphi2, Azero2))));
   coefficients.push_back(
      TCoefficient("MassB_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassB, Azero))));
   coefficients.push_back(
      TCoefficient("MassB_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassB, m12))));
   coefficients.push_back(
      TCoefficient("MassWB_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassWB, Azero))));
   coefficients.push_back(
      TCoefficient("MassWB_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassWB, m12))));
   coefficients.push_back(
      TCoefficient("MassG_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassG, Azero))));
   coefficients.push_back(
      TCoefficient("MassG_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassG, m12))));
   coefficients.push_back(
      TCoefficient("MassBp_Azero_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassBp, Azero))));
   coefficients.push_back(
      TCoefficient("MassBp_m12_coeff", 0, 1, 1,
                   to_valarray(COEFFICIENT(MassBp, m12))));
}

} // namespace flexiblesusy

#endif
