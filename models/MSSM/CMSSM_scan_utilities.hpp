// ====================================================================
// Routines used in scan class
// ====================================================================

#ifndef CMSSM_SCAN_UTILITIES_H
#define CMSSM_SCAN_UTILITIES_H

#include "CMSSM_two_scale_model_slha.hpp"
#include "CMSSM_semi_two_scale_model_slha.hpp"
#include "MSSM_info.hpp"
#include "MSSM_utilities.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

#define PHYSICAL(p) model.get_physical().p
#define PHYSICALSLHA(p) model.get_physical_slha().p
#define DRBARSLHA(p) model.get_drbar_slha().p
#define MODELPARAMETER(p) model.get_##p()


namespace flexiblesusy {

void write_CMSSM_inputs(const CMSSM_input_parameters<Two_scale>&, std::ostream &, std::size_t);
void write_CMSSM_inputs_list(std::ostream &, std::size_t);
void write_CMSSM_inputs(const CMSSM_semianalytic_input_parameters<Two_scale>&, std::ostream &, std::size_t);
void write_CMSSM_semianalytic_inputs_list(std::ostream &, std::size_t);
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
std::valarray<double> to_valarray(const Eigen::Matrix<std::complex<double>, M, N> & v)
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

class CMSSM_pole_mass_writer {
public:
   CMSSM_pole_mass_writer();
   ~CMSSM_pole_mass_writer() {}

   template <class T>
   void extract_pole_masses(const CMSSM<T>&);
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

   CMSSM_input_parameters<Two_scale> pole_masses_inputs;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> pole_masses_problems;

   double pole_masses_scale;
   unsigned width;
};

class CMSSM_semianalytic_pole_mass_writer {
public:
   CMSSM_semianalytic_pole_mass_writer();
   ~CMSSM_semianalytic_pole_mass_writer() {}

   template <class T>
   void extract_pole_masses(const CMSSM_semianalytic<T>&);
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

   CMSSM_semianalytic_input_parameters<Two_scale> pole_masses_inputs;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> pole_masses_problems;

   double pole_masses_scale;
   double m0Sqr;
   unsigned width;
};

class CMSSM_drbar_values_writer {
public:
   CMSSM_drbar_values_writer();
   ~CMSSM_drbar_values_writer() {}

   template <class T>
   void extract_drbar_masses(const CMSSM<T>&);
   template <class T>
   void extract_drbar_susy_pars(const CMSSM<T>&);
   template <class T>
   void extract_drbar_soft_pars(const CMSSM<T>&);
   template <class T>
   void extract_drbar_mixings(const CMSSM<T>&);

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
      TMixing(const std::string name_, std::size_t dimension_,
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

   CMSSM_input_parameters<Two_scale> drbar_masses_inputs;
   CMSSM_input_parameters<Two_scale> drbar_susy_pars_inputs;
   CMSSM_input_parameters<Two_scale> drbar_soft_pars_inputs;
   CMSSM_input_parameters<Two_scale> drbar_mixings_inputs;

   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_masses_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_susy_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_soft_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_mixings_problems;

   double drbar_masses_scale;
   double drbar_susy_pars_scale;
   double drbar_soft_pars_scale;
   double drbar_mixings_scale;
   unsigned width;
};

class CMSSM_semianalytic_drbar_values_writer {
public:
   CMSSM_semianalytic_drbar_values_writer();
   ~CMSSM_semianalytic_drbar_values_writer() {}

   template <class T>
   void extract_drbar_masses(const CMSSM_semianalytic<T>&);
   template <class T>
   void extract_drbar_susy_pars(const CMSSM_semianalytic<T>&);
   template <class T>
   void extract_drbar_soft_pars(const CMSSM_semianalytic<T>&);
   template <class T>
   void extract_drbar_mixings(const CMSSM_semianalytic<T>&);

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
      TMixing(const std::string name_, std::size_t dimension_,
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

   CMSSM_semianalytic_input_parameters<Two_scale> drbar_masses_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> drbar_susy_pars_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> drbar_soft_pars_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> drbar_mixings_inputs;

   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_masses_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_susy_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_soft_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> drbar_mixings_problems;

   double drbar_masses_scale;
   double drbar_susy_pars_scale;
   double drbar_soft_pars_scale;
   double drbar_mixings_scale;
   unsigned width;
};

/// As above, but in SLHA convention
class CMSSM_slha_values_writer {
public:
   CMSSM_slha_values_writer();
   ~CMSSM_slha_values_writer() {}

   void set_high_scale(double high_scale_) { high_scale = high_scale_; }
   void set_susy_scale(double susy_scale_) { susy_scale = susy_scale_; }
   void set_low_scale(double low_scale_) { low_scale = low_scale_; }

   template <class T>
   void extract_slha_pole_masses(const CMSSM_slha<T>&);
   template <class T>
   void extract_slha_running_masses(const CMSSM_slha<T>&);
   template <class T>
   void extract_slha_susy_pars(const CMSSM_slha<T>&);
   template <class T>
   void extract_slha_soft_pars(const CMSSM_slha<T>&);
   template <class T>
   void extract_slha_pole_mixings(const CMSSM_slha<T>&);
   template <class T>
   void extract_slha_running_mixings(const CMSSM_slha<T>&);

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
      TMixing(const std::string name_, std::size_t dimension_,
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

   CMSSM_input_parameters<Two_scale> slha_pole_masses_inputs;
   CMSSM_input_parameters<Two_scale> slha_running_masses_inputs;
   CMSSM_input_parameters<Two_scale> slha_susy_pars_inputs;
   CMSSM_input_parameters<Two_scale> slha_soft_pars_inputs;
   CMSSM_input_parameters<Two_scale> slha_pole_mixings_inputs;
   CMSSM_input_parameters<Two_scale> slha_running_mixings_inputs;

   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_pole_masses_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_running_masses_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_susy_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_soft_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_pole_mixings_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_running_mixings_problems;

   double high_scale;
   double susy_scale;
   double low_scale;
   unsigned width;
};

class CMSSM_semianalytic_slha_values_writer {
public:
   CMSSM_semianalytic_slha_values_writer();
   ~CMSSM_semianalytic_slha_values_writer() {}

   void set_high_scale(double high_scale_) { high_scale = high_scale_; }
   void set_susy_scale(double susy_scale_) { susy_scale = susy_scale_; }
   void set_low_scale(double low_scale_) { low_scale = low_scale_; }

   template <class T>
   void extract_slha_pole_masses(const CMSSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_running_masses(const CMSSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_susy_pars(const CMSSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_soft_pars(const CMSSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_pole_mixings(const CMSSM_semianalytic_slha<T>&);
   template <class T>
   void extract_slha_running_mixings(const CMSSM_semianalytic_slha<T>&);

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
      TMixing(const std::string name_, std::size_t dimension_,
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

   CMSSM_semianalytic_input_parameters<Two_scale> slha_pole_masses_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> slha_running_masses_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> slha_susy_pars_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> slha_soft_pars_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> slha_pole_mixings_inputs;
   CMSSM_semianalytic_input_parameters<Two_scale> slha_running_mixings_inputs;

   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_pole_masses_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_running_masses_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_susy_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_soft_pars_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_pole_mixings_problems;
   Problems<MSSM_info::NUMBER_OF_PARTICLES> slha_running_mixings_problems;

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

template <class T>
void CMSSM_pole_mass_writer::extract_pole_masses(const CMSSM<T>& model)
{
   pole_masses.clear();
   pole_masses_scale = model.get_scale();
   pole_masses_inputs = model.get_input();
   pole_masses_problems = model.get_problems();

   pole_masses.push_back(TPoleMass("MGlu", to_valarray(PHYSICAL(MGlu))));
   pole_masses.push_back(TPoleMass("MSd", to_valarray(PHYSICAL(MSd))));
   pole_masses.push_back(TPoleMass("MSv", to_valarray(PHYSICAL(MSv))));
   pole_masses.push_back(TPoleMass("MSu", to_valarray(PHYSICAL(MSu))));
   pole_masses.push_back(TPoleMass("MSe", to_valarray(PHYSICAL(MSe))));
   pole_masses.push_back(TPoleMass("Mhh", to_valarray(PHYSICAL(Mhh))));
   pole_masses.push_back(TPoleMass("MAh", to_valarray(PHYSICAL(MAh))));
   pole_masses.push_back(TPoleMass("MHpm", to_valarray(PHYSICAL(MHpm))));
   pole_masses.push_back(TPoleMass("MChi", to_valarray(PHYSICAL(MChi))));
   pole_masses.push_back(TPoleMass("MCha", to_valarray(PHYSICAL(MCha))));

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
void CMSSM_semianalytic_pole_mass_writer::extract_pole_masses(const CMSSM_semianalytic<T>& model)
{
   pole_masses.clear();
   pole_masses_scale = model.get_scale();
   pole_masses_inputs = model.get_input();
   pole_masses_problems = model.get_problems();
   m0Sqr = model.get_ewsb_output_parameter(0);

   pole_masses.push_back(TPoleMass("MGlu", to_valarray(PHYSICAL(MGlu))));
   pole_masses.push_back(TPoleMass("MSd", to_valarray(PHYSICAL(MSd))));
   pole_masses.push_back(TPoleMass("MSv", to_valarray(PHYSICAL(MSv))));
   pole_masses.push_back(TPoleMass("MSu", to_valarray(PHYSICAL(MSu))));
   pole_masses.push_back(TPoleMass("MSe", to_valarray(PHYSICAL(MSe))));
   pole_masses.push_back(TPoleMass("Mhh", to_valarray(PHYSICAL(Mhh))));
   pole_masses.push_back(TPoleMass("MAh", to_valarray(PHYSICAL(MAh))));
   pole_masses.push_back(TPoleMass("MHpm", to_valarray(PHYSICAL(MHpm))));
   pole_masses.push_back(TPoleMass("MChi", to_valarray(PHYSICAL(MChi))));
   pole_masses.push_back(TPoleMass("MCha", to_valarray(PHYSICAL(MCha))));

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
void CMSSM_drbar_values_writer::extract_drbar_masses(const CMSSM<T>& model)
{
   drbar_masses.clear();
   drbar_masses_scale = model.get_scale();
   drbar_masses_inputs = model.get_input();
   drbar_masses_problems = model.get_problems();

   drbar_masses.push_back(TMass("MVG", to_valarray(MODELPARAMETER(MVG))));
   drbar_masses.push_back(TMass("MGlu", to_valarray(MODELPARAMETER(MGlu))));
   drbar_masses.push_back(TMass("MFv", to_valarray(MODELPARAMETER(MFv))));
   drbar_masses.push_back(TMass("MVP", to_valarray(MODELPARAMETER(MVP))));
   drbar_masses.push_back(TMass("MVZ", to_valarray(MODELPARAMETER(MVZ))));
   drbar_masses.push_back(TMass("MSd", to_valarray(MODELPARAMETER(MSd))));
   drbar_masses.push_back(TMass("MSv", to_valarray(MODELPARAMETER(MSv))));
   drbar_masses.push_back(TMass("MSu", to_valarray(MODELPARAMETER(MSu))));
   drbar_masses.push_back(TMass("MSe", to_valarray(MODELPARAMETER(MSe))));
   drbar_masses.push_back(TMass("Mhh", to_valarray(MODELPARAMETER(Mhh))));
   drbar_masses.push_back(TMass("MAh", to_valarray(MODELPARAMETER(MAh))));
   drbar_masses.push_back(TMass("MHpm", to_valarray(MODELPARAMETER(MHpm))));
   drbar_masses.push_back(TMass("MChi", to_valarray(MODELPARAMETER(MChi))));
   drbar_masses.push_back(TMass("MCha", to_valarray(MODELPARAMETER(MCha))));
   drbar_masses.push_back(TMass("MFe", to_valarray(MODELPARAMETER(MFe))));
   drbar_masses.push_back(TMass("MFd", to_valarray(MODELPARAMETER(MFd))));
   drbar_masses.push_back(TMass("MFu", to_valarray(MODELPARAMETER(MFu))));
   drbar_masses.push_back(TMass("MVWm", to_valarray(MODELPARAMETER(MVWm))));
}

template <class T>
void CMSSM_semianalytic_drbar_values_writer::extract_drbar_masses(const CMSSM_semianalytic<T>& model)
{
   drbar_masses.clear();
   drbar_masses_scale = model.get_scale();
   drbar_masses_inputs = model.get_input();
   drbar_masses_problems = model.get_problems();

   drbar_masses.push_back(TMass("MVG", to_valarray(MODELPARAMETER(MVG))));
   drbar_masses.push_back(TMass("MGlu", to_valarray(MODELPARAMETER(MGlu))));
   drbar_masses.push_back(TMass("MFv", to_valarray(MODELPARAMETER(MFv))));
   drbar_masses.push_back(TMass("MVP", to_valarray(MODELPARAMETER(MVP))));
   drbar_masses.push_back(TMass("MVZ", to_valarray(MODELPARAMETER(MVZ))));
   drbar_masses.push_back(TMass("MSd", to_valarray(MODELPARAMETER(MSd))));
   drbar_masses.push_back(TMass("MSv", to_valarray(MODELPARAMETER(MSv))));
   drbar_masses.push_back(TMass("MSu", to_valarray(MODELPARAMETER(MSu))));
   drbar_masses.push_back(TMass("MSe", to_valarray(MODELPARAMETER(MSe))));
   drbar_masses.push_back(TMass("Mhh", to_valarray(MODELPARAMETER(Mhh))));
   drbar_masses.push_back(TMass("MAh", to_valarray(MODELPARAMETER(MAh))));
   drbar_masses.push_back(TMass("MHpm", to_valarray(MODELPARAMETER(MHpm))));
   drbar_masses.push_back(TMass("MChi", to_valarray(MODELPARAMETER(MChi))));
   drbar_masses.push_back(TMass("MCha", to_valarray(MODELPARAMETER(MCha))));
   drbar_masses.push_back(TMass("MFe", to_valarray(MODELPARAMETER(MFe))));
   drbar_masses.push_back(TMass("MFd", to_valarray(MODELPARAMETER(MFd))));
   drbar_masses.push_back(TMass("MFu", to_valarray(MODELPARAMETER(MFu))));
   drbar_masses.push_back(TMass("MVWm", to_valarray(MODELPARAMETER(MVWm))));
}

template <class T>
void CMSSM_drbar_values_writer::extract_drbar_susy_pars(const CMSSM<T>& model)
{
   drbar_susy_pars.clear();
   drbar_susy_pars_scale = model.get_scale();
   drbar_susy_pars_inputs = model.get_input();
   drbar_susy_pars_problems = model.get_problems();

   drbar_susy_pars.push_back(TParameter("Yd", 0, 3, 3, to_valarray(MODELPARAMETER(Yd))));
   drbar_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(MODELPARAMETER(Ye))));
   drbar_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(MODELPARAMETER(Yu))));
   drbar_susy_pars.push_back(TParameter("Mu", 1, 1, 1, to_valarray(MODELPARAMETER(Mu))));
   drbar_susy_pars.push_back(TParameter("g1", 0, 1, 1, to_valarray(MODELPARAMETER(g1))));
   drbar_susy_pars.push_back(TParameter("g2", 0, 1, 1, to_valarray(MODELPARAMETER(g2))));
   drbar_susy_pars.push_back(TParameter("g3", 0, 1, 1, to_valarray(MODELPARAMETER(g3))));
   drbar_susy_pars.push_back(TParameter("vd", 1, 1, 1, to_valarray(MODELPARAMETER(vd))));
   drbar_susy_pars.push_back(TParameter("vu", 1, 1, 1, to_valarray(MODELPARAMETER(vu))));
}

template <class T>
void CMSSM_semianalytic_drbar_values_writer::extract_drbar_susy_pars(const CMSSM_semianalytic<T>& model)
{
   drbar_susy_pars.clear();
   drbar_susy_pars_scale = model.get_scale();
   drbar_susy_pars_inputs = model.get_input();
   drbar_susy_pars_problems = model.get_problems();

   drbar_susy_pars.push_back(TParameter("Yd", 0, 3, 3, to_valarray(MODELPARAMETER(Yd))));
   drbar_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(MODELPARAMETER(Ye))));
   drbar_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(MODELPARAMETER(Yu))));
   drbar_susy_pars.push_back(TParameter("Mu", 1, 1, 1, to_valarray(MODELPARAMETER(Mu))));
   drbar_susy_pars.push_back(TParameter("g1", 0, 1, 1, to_valarray(MODELPARAMETER(g1))));
   drbar_susy_pars.push_back(TParameter("g2", 0, 1, 1, to_valarray(MODELPARAMETER(g2))));
   drbar_susy_pars.push_back(TParameter("g3", 0, 1, 1, to_valarray(MODELPARAMETER(g3))));
   drbar_susy_pars.push_back(TParameter("vd", 1, 1, 1, to_valarray(MODELPARAMETER(vd))));
   drbar_susy_pars.push_back(TParameter("vu", 1, 1, 1, to_valarray(MODELPARAMETER(vu))));
}

template <class T>
void CMSSM_drbar_values_writer::extract_drbar_soft_pars(const CMSSM<T>& model)
{
   drbar_soft_pars.clear();
   drbar_soft_pars_scale = model.get_scale();
   drbar_soft_pars_inputs = model.get_input();
   drbar_soft_pars_problems = model.get_problems();

   drbar_soft_pars.push_back(TParameter("TYd", 1, 3, 3, to_valarray(MODELPARAMETER(TYd))));
   drbar_soft_pars.push_back(TParameter("TYe", 1, 3, 3, to_valarray(MODELPARAMETER(TYe))));
   drbar_soft_pars.push_back(TParameter("TYu", 1, 3, 3, to_valarray(MODELPARAMETER(TYu))));
   drbar_soft_pars.push_back(TParameter("BMu", 2, 1, 1, to_valarray(MODELPARAMETER(BMu))));
   drbar_soft_pars.push_back(TParameter("mq2", 2, 3, 3, to_valarray(MODELPARAMETER(mq2))));
   drbar_soft_pars.push_back(TParameter("ml2", 2, 3, 3, to_valarray(MODELPARAMETER(ml2))));
   drbar_soft_pars.push_back(TParameter("mHd2", 2, 1, 1, to_valarray(MODELPARAMETER(mHd2))));
   drbar_soft_pars.push_back(TParameter("mHu2", 2, 1, 1, to_valarray(MODELPARAMETER(mHu2))));
   drbar_soft_pars.push_back(TParameter("md2", 2, 3, 3, to_valarray(MODELPARAMETER(md2))));
   drbar_soft_pars.push_back(TParameter("mu2", 2, 3, 3, to_valarray(MODELPARAMETER(mu2))));
   drbar_soft_pars.push_back(TParameter("me2", 2, 3, 3, to_valarray(MODELPARAMETER(me2))));
   drbar_soft_pars.push_back(TParameter("MassB", 1, 1, 1, to_valarray(MODELPARAMETER(MassB))));
   drbar_soft_pars.push_back(TParameter("MassWB", 1, 1, 1, to_valarray(MODELPARAMETER(MassWB))));
   drbar_soft_pars.push_back(TParameter("MassG", 1, 1, 1, to_valarray(MODELPARAMETER(MassG))));
}

template <class T>
void CMSSM_semianalytic_drbar_values_writer::extract_drbar_soft_pars(const CMSSM_semianalytic<T>& model)
{
   drbar_soft_pars.clear();
   drbar_soft_pars_scale = model.get_scale();
   drbar_soft_pars_inputs = model.get_input();
   drbar_soft_pars_problems = model.get_problems();

   drbar_soft_pars.push_back(TParameter("TYd", 1, 3, 3, to_valarray(MODELPARAMETER(TYd))));
   drbar_soft_pars.push_back(TParameter("TYe", 1, 3, 3, to_valarray(MODELPARAMETER(TYe))));
   drbar_soft_pars.push_back(TParameter("TYu", 1, 3, 3, to_valarray(MODELPARAMETER(TYu))));
   drbar_soft_pars.push_back(TParameter("BMu", 2, 1, 1, to_valarray(MODELPARAMETER(BMu))));
   drbar_soft_pars.push_back(TParameter("mq2", 2, 3, 3, to_valarray(MODELPARAMETER(mq2))));
   drbar_soft_pars.push_back(TParameter("ml2", 2, 3, 3, to_valarray(MODELPARAMETER(ml2))));
   drbar_soft_pars.push_back(TParameter("mHd2", 2, 1, 1, to_valarray(MODELPARAMETER(mHd2))));
   drbar_soft_pars.push_back(TParameter("mHu2", 2, 1, 1, to_valarray(MODELPARAMETER(mHu2))));
   drbar_soft_pars.push_back(TParameter("md2", 2, 3, 3, to_valarray(MODELPARAMETER(md2))));
   drbar_soft_pars.push_back(TParameter("mu2", 2, 3, 3, to_valarray(MODELPARAMETER(mu2))));
   drbar_soft_pars.push_back(TParameter("me2", 2, 3, 3, to_valarray(MODELPARAMETER(me2))));
   drbar_soft_pars.push_back(TParameter("MassB", 1, 1, 1, to_valarray(MODELPARAMETER(MassB))));
   drbar_soft_pars.push_back(TParameter("MassWB", 1, 1, 1, to_valarray(MODELPARAMETER(MassWB))));
   drbar_soft_pars.push_back(TParameter("MassG", 1, 1, 1, to_valarray(MODELPARAMETER(MassG))));
}

template <class T>
void CMSSM_drbar_values_writer::extract_drbar_mixings(const CMSSM<T>& model)
{
   drbar_mixings.clear();
   drbar_mixings_scale = model.get_scale();
   drbar_mixings_inputs = model.get_input();
   drbar_mixings_problems = model.get_problems();

   drbar_mixings.push_back(TMixing("ZD", 6, true, to_valarray(MODELPARAMETER(ZD))));
   drbar_mixings.push_back(TMixing("ZV", 3, true, to_valarray(MODELPARAMETER(ZV))));
   drbar_mixings.push_back(TMixing("ZU", 6, true, to_valarray(MODELPARAMETER(ZU))));
   drbar_mixings.push_back(TMixing("ZE", 6, true, to_valarray(MODELPARAMETER(ZE))));
   drbar_mixings.push_back(TMixing("ZH", 2, true, to_valarray(MODELPARAMETER(ZH))));
   drbar_mixings.push_back(TMixing("ZA", 2, true, to_valarray(MODELPARAMETER(ZA))));
   drbar_mixings.push_back(TMixing("ZP", 2, true, to_valarray(MODELPARAMETER(ZP))));
   drbar_mixings.push_back(TMixing("ZN", 4, false, to_valarray(MODELPARAMETER(ZN))));
   drbar_mixings.push_back(TMixing("UM", 2, false, to_valarray(MODELPARAMETER(UM))));
   drbar_mixings.push_back(TMixing("UP", 2, false, to_valarray(MODELPARAMETER(UP))));
   drbar_mixings.push_back(TMixing("ZEL", 3, false, to_valarray(MODELPARAMETER(ZEL))));
   drbar_mixings.push_back(TMixing("ZER", 3, false, to_valarray(MODELPARAMETER(ZER))));
   drbar_mixings.push_back(TMixing("ZDL", 3, false, to_valarray(MODELPARAMETER(ZDL))));
   drbar_mixings.push_back(TMixing("ZDR", 3, false, to_valarray(MODELPARAMETER(ZDR))));
   drbar_mixings.push_back(TMixing("ZUL", 3, false, to_valarray(MODELPARAMETER(ZUL))));
   drbar_mixings.push_back(TMixing("ZUR", 3, false, to_valarray(MODELPARAMETER(ZUR))));
}

template <class T>
void CMSSM_semianalytic_drbar_values_writer::extract_drbar_mixings(const CMSSM_semianalytic<T>& model)
{
   drbar_mixings.clear();
   drbar_mixings_scale = model.get_scale();
   drbar_mixings_inputs = model.get_input();
   drbar_mixings_problems = model.get_problems();

   drbar_mixings.push_back(TMixing("ZD", 6, true, to_valarray(MODELPARAMETER(ZD))));
   drbar_mixings.push_back(TMixing("ZV", 3, true, to_valarray(MODELPARAMETER(ZV))));
   drbar_mixings.push_back(TMixing("ZU", 6, true, to_valarray(MODELPARAMETER(ZU))));
   drbar_mixings.push_back(TMixing("ZE", 6, true, to_valarray(MODELPARAMETER(ZE))));
   drbar_mixings.push_back(TMixing("ZH", 2, true, to_valarray(MODELPARAMETER(ZH))));
   drbar_mixings.push_back(TMixing("ZA", 2, true, to_valarray(MODELPARAMETER(ZA))));
   drbar_mixings.push_back(TMixing("ZP", 2, true, to_valarray(MODELPARAMETER(ZP))));
   drbar_mixings.push_back(TMixing("ZN", 4, false, to_valarray(MODELPARAMETER(ZN))));
   drbar_mixings.push_back(TMixing("UM", 2, false, to_valarray(MODELPARAMETER(UM))));
   drbar_mixings.push_back(TMixing("UP", 2, false, to_valarray(MODELPARAMETER(UP))));
   drbar_mixings.push_back(TMixing("ZEL", 3, false, to_valarray(MODELPARAMETER(ZEL))));
   drbar_mixings.push_back(TMixing("ZER", 3, false, to_valarray(MODELPARAMETER(ZER))));
   drbar_mixings.push_back(TMixing("ZDL", 3, false, to_valarray(MODELPARAMETER(ZDL))));
   drbar_mixings.push_back(TMixing("ZDR", 3, false, to_valarray(MODELPARAMETER(ZDR))));
   drbar_mixings.push_back(TMixing("ZUL", 3, false, to_valarray(MODELPARAMETER(ZUL))));
   drbar_mixings.push_back(TMixing("ZUR", 3, false, to_valarray(MODELPARAMETER(ZUR))));
}

template <class T>
void CMSSM_slha_values_writer::extract_slha_pole_masses(const CMSSM_slha<T>& model)
{
   slha_pole_masses.clear();
   slha_pole_masses_inputs = model.get_input();
   slha_pole_masses_problems = model.get_problems();

   slha_pole_masses.push_back(TMass("MGlu", to_valarray(PHYSICALSLHA(MGlu))));
   slha_pole_masses.push_back(TMass("MSd", to_valarray(PHYSICALSLHA(MSd))));
   slha_pole_masses.push_back(TMass("MSv", to_valarray(PHYSICALSLHA(MSv))));
   slha_pole_masses.push_back(TMass("MSu", to_valarray(PHYSICALSLHA(MSu))));
   slha_pole_masses.push_back(TMass("MSe", to_valarray(PHYSICALSLHA(MSe))));
   slha_pole_masses.push_back(TMass("Mhh", to_valarray(PHYSICALSLHA(Mhh))));
   slha_pole_masses.push_back(TMass("MAh", to_valarray(PHYSICALSLHA(MAh))));
   slha_pole_masses.push_back(TMass("MHpm", to_valarray(PHYSICALSLHA(MHpm))));
   slha_pole_masses.push_back(TMass("MChi", to_valarray(PHYSICALSLHA(MChi))));
   slha_pole_masses.push_back(TMass("MCha", to_valarray(PHYSICALSLHA(MCha))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_pole_masses.push_back(TMass("MFd", to_valarray(PHYSICALSLHA(MFd))));
      slha_pole_masses.push_back(TMass("MFe", to_valarray(PHYSICALSLHA(MFe))));
      slha_pole_masses.push_back(TMass("MFu", to_valarray(PHYSICALSLHA(MFu))));
      slha_pole_masses.push_back(TMass("MFv", to_valarray(PHYSICALSLHA(MFv))));
      slha_pole_masses.push_back(TMass("MVG", to_valarray(PHYSICALSLHA(MVG))));
      slha_pole_masses.push_back(TMass("MVP", to_valarray(PHYSICALSLHA(MVP))));
      slha_pole_masses.push_back(TMass("MVWm", to_valarray(PHYSICALSLHA(MVWm))));
      slha_pole_masses.push_back(TMass("MVZ", to_valarray(PHYSICALSLHA(MVZ))));

   }
}

template <class T>
void CMSSM_semianalytic_slha_values_writer::extract_slha_pole_masses(const CMSSM_semianalytic_slha<T>& model)
{
   slha_pole_masses.clear();
   slha_pole_masses_inputs = model.get_input();
   slha_pole_masses_problems = model.get_problems();
   slha_pole_masses_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_pole_masses.push_back(TMass("MGlu", to_valarray(PHYSICALSLHA(MGlu))));
   slha_pole_masses.push_back(TMass("MSd", to_valarray(PHYSICALSLHA(MSd))));
   slha_pole_masses.push_back(TMass("MSv", to_valarray(PHYSICALSLHA(MSv))));
   slha_pole_masses.push_back(TMass("MSu", to_valarray(PHYSICALSLHA(MSu))));
   slha_pole_masses.push_back(TMass("MSe", to_valarray(PHYSICALSLHA(MSe))));
   slha_pole_masses.push_back(TMass("Mhh", to_valarray(PHYSICALSLHA(Mhh))));
   slha_pole_masses.push_back(TMass("MAh", to_valarray(PHYSICALSLHA(MAh))));
   slha_pole_masses.push_back(TMass("MHpm", to_valarray(PHYSICALSLHA(MHpm))));
   slha_pole_masses.push_back(TMass("MChi", to_valarray(PHYSICALSLHA(MChi))));
   slha_pole_masses.push_back(TMass("MCha", to_valarray(PHYSICALSLHA(MCha))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_pole_masses.push_back(TMass("MFd", to_valarray(PHYSICALSLHA(MFd))));
      slha_pole_masses.push_back(TMass("MFe", to_valarray(PHYSICALSLHA(MFe))));
      slha_pole_masses.push_back(TMass("MFu", to_valarray(PHYSICALSLHA(MFu))));
      slha_pole_masses.push_back(TMass("MFv", to_valarray(PHYSICALSLHA(MFv))));
      slha_pole_masses.push_back(TMass("MVG", to_valarray(PHYSICALSLHA(MVG))));
      slha_pole_masses.push_back(TMass("MVP", to_valarray(PHYSICALSLHA(MVP))));
      slha_pole_masses.push_back(TMass("MVWm", to_valarray(PHYSICALSLHA(MVWm))));
      slha_pole_masses.push_back(TMass("MVZ", to_valarray(PHYSICALSLHA(MVZ))));

   }
}

template <class T>
void CMSSM_slha_values_writer::extract_slha_running_masses(const CMSSM_slha<T>& model)
{
   slha_running_masses.clear();
   slha_running_masses_inputs = model.get_input();
   slha_running_masses_problems = model.get_problems();

   slha_running_masses.push_back(TMass("DRbarMGlu", to_valarray(DRBARSLHA(MGlu))));
   slha_running_masses.push_back(TMass("DRbarMSd", to_valarray(DRBARSLHA(MSd))));
   slha_running_masses.push_back(TMass("DRbarMSv", to_valarray(DRBARSLHA(MSv))));
   slha_running_masses.push_back(TMass("DRbarMSu", to_valarray(DRBARSLHA(MSu))));
   slha_running_masses.push_back(TMass("DRbarMSe", to_valarray(DRBARSLHA(MSe))));
   slha_running_masses.push_back(TMass("DRbarMhh", to_valarray(DRBARSLHA(Mhh))));
   slha_running_masses.push_back(TMass("DRbarMAh", to_valarray(DRBARSLHA(MAh))));
   slha_running_masses.push_back(TMass("DRbarMHpm", to_valarray(DRBARSLHA(MHpm))));
   slha_running_masses.push_back(TMass("DRbarMChi", to_valarray(DRBARSLHA(MChi))));
   slha_running_masses.push_back(TMass("DRbarMCha", to_valarray(DRBARSLHA(MCha))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_running_masses.push_back(TMass("DRbarMFd", to_valarray(DRBARSLHA(MFd))));
      slha_running_masses.push_back(TMass("DRbarMFe", to_valarray(DRBARSLHA(MFe))));
      slha_running_masses.push_back(TMass("DRbarMFu", to_valarray(DRBARSLHA(MFu))));
      slha_running_masses.push_back(TMass("DRbarMFv", to_valarray(DRBARSLHA(MFv))));
      slha_running_masses.push_back(TMass("DRbarMVG", to_valarray(DRBARSLHA(MVG))));
      slha_running_masses.push_back(TMass("DRbarMVP", to_valarray(DRBARSLHA(MVP))));
      slha_running_masses.push_back(TMass("DRbarMVWm", to_valarray(DRBARSLHA(MVWm))));
      slha_running_masses.push_back(TMass("DRbarMVZ", to_valarray(DRBARSLHA(MVZ))));

   }
}

template <class T>
void CMSSM_semianalytic_slha_values_writer::extract_slha_running_masses(const CMSSM_semianalytic_slha<T>& model)
{
   slha_running_masses.clear();
   slha_running_masses_inputs = model.get_input();
   slha_running_masses_problems = model.get_problems();
   slha_running_masses_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_running_masses.push_back(TMass("DRbarMGlu", to_valarray(DRBARSLHA(MGlu))));
   slha_running_masses.push_back(TMass("DRbarMSd", to_valarray(DRBARSLHA(MSd))));
   slha_running_masses.push_back(TMass("DRbarMSv", to_valarray(DRBARSLHA(MSv))));
   slha_running_masses.push_back(TMass("DRbarMSu", to_valarray(DRBARSLHA(MSu))));
   slha_running_masses.push_back(TMass("DRbarMSe", to_valarray(DRBARSLHA(MSe))));
   slha_running_masses.push_back(TMass("DRbarMhh", to_valarray(DRBARSLHA(Mhh))));
   slha_running_masses.push_back(TMass("DRbarMAh", to_valarray(DRBARSLHA(MAh))));
   slha_running_masses.push_back(TMass("DRbarMHpm", to_valarray(DRBARSLHA(MHpm))));
   slha_running_masses.push_back(TMass("DRbarMChi", to_valarray(DRBARSLHA(MChi))));
   slha_running_masses.push_back(TMass("DRbarMCha", to_valarray(DRBARSLHA(MCha))));

   if (model.do_calculate_sm_pole_masses()) {
      slha_running_masses.push_back(TMass("DRbarMFd", to_valarray(DRBARSLHA(MFd))));
      slha_running_masses.push_back(TMass("DRbarMFe", to_valarray(DRBARSLHA(MFe))));
      slha_running_masses.push_back(TMass("DRbarMFu", to_valarray(DRBARSLHA(MFu))));
      slha_running_masses.push_back(TMass("DRbarMFv", to_valarray(DRBARSLHA(MFv))));
      slha_running_masses.push_back(TMass("DRbarMVG", to_valarray(DRBARSLHA(MVG))));
      slha_running_masses.push_back(TMass("DRbarMVP", to_valarray(DRBARSLHA(MVP))));
      slha_running_masses.push_back(TMass("DRbarMVWm", to_valarray(DRBARSLHA(MVWm))));
      slha_running_masses.push_back(TMass("DRbarMVZ", to_valarray(DRBARSLHA(MVZ))));

   }
}

template <class T>
void CMSSM_slha_values_writer::extract_slha_susy_pars(const CMSSM_slha<T>& model)
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
   slha_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(Ye_matrix)));
   slha_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(Yu_matrix)));
   slha_susy_pars.push_back(TParameter("Mu", 1, 1, 1, to_valarray(MODELPARAMETER(Mu))));
   slha_susy_pars.push_back(TParameter("gY", 0, 1, 1, to_valarray(MODELPARAMETER(g1) * 0.7745966692414834)));
   slha_susy_pars.push_back(TParameter("g2", 0, 1, 1, to_valarray(MODELPARAMETER(g2))));
   slha_susy_pars.push_back(TParameter("g3", 0, 1, 1, to_valarray(MODELPARAMETER(g3))));
   slha_susy_pars.push_back(TParameter("vd", 1, 1, 1, to_valarray(MODELPARAMETER(vd))));
   slha_susy_pars.push_back(TParameter("vu", 1, 1, 1, to_valarray(MODELPARAMETER(vu))));
   slha_susy_pars.push_back(TParameter("Beta", 0, 1, 1,
                                       to_valarray(ArcTan((MODELPARAMETER(vu)) / (MODELPARAMETER(vd))))));
}

template <class T>
void CMSSM_semianalytic_slha_values_writer::extract_slha_susy_pars(const CMSSM_semianalytic_slha<T>& model)
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
   slha_susy_pars.push_back(TParameter("Ye", 0, 3, 3, to_valarray(Ye_matrix)));
   slha_susy_pars.push_back(TParameter("Yu", 0, 3, 3, to_valarray(Yu_matrix)));
   slha_susy_pars.push_back(TParameter("Mu", 1, 1, 1, to_valarray(MODELPARAMETER(Mu))));
   slha_susy_pars.push_back(TParameter("gY", 0, 1, 1, to_valarray(MODELPARAMETER(g1) * 0.7745966692414834)));
   slha_susy_pars.push_back(TParameter("g2", 0, 1, 1, to_valarray(MODELPARAMETER(g2))));
   slha_susy_pars.push_back(TParameter("g3", 0, 1, 1, to_valarray(MODELPARAMETER(g3))));
   slha_susy_pars.push_back(TParameter("vd", 1, 1, 1, to_valarray(MODELPARAMETER(vd))));
   slha_susy_pars.push_back(TParameter("vu", 1, 1, 1, to_valarray(MODELPARAMETER(vu))));
   slha_susy_pars.push_back(TParameter("Beta", 0, 1, 1,
                                       to_valarray(ArcTan((MODELPARAMETER(vu)) / (MODELPARAMETER(vd))))));
}

template <class T>
void CMSSM_slha_values_writer::extract_slha_soft_pars(const CMSSM_slha<T>& model)
{
   slha_soft_pars.clear();
   slha_soft_pars_inputs = model.get_input();
   slha_soft_pars_problems = model.get_problems();

   slha_soft_pars.push_back(TParameter("TYd", 1, 3, 3, to_valarray(MODELPARAMETER(TYd_slha))));
   slha_soft_pars.push_back(TParameter("TYe", 1, 3, 3, to_valarray(MODELPARAMETER(TYe_slha))));
   slha_soft_pars.push_back(TParameter("TYu", 1, 3, 3, to_valarray(MODELPARAMETER(TYu_slha))));
   slha_soft_pars.push_back(TParameter("BMu", 2, 1, 1, to_valarray(MODELPARAMETER(BMu))));
   slha_soft_pars.push_back(TParameter("mq2", 2, 3, 3, to_valarray(MODELPARAMETER(mq2_slha))));
   slha_soft_pars.push_back(TParameter("ml2", 2, 3, 3, to_valarray(MODELPARAMETER(ml2_slha))));
   slha_soft_pars.push_back(TParameter("mHd2", 2, 1, 1, to_valarray(MODELPARAMETER(mHd2))));
   slha_soft_pars.push_back(TParameter("mHu2", 2, 1, 1, to_valarray(MODELPARAMETER(mHu2))));
   slha_soft_pars.push_back(TParameter("md2", 2, 3, 3, to_valarray(MODELPARAMETER(md2_slha))));
   slha_soft_pars.push_back(TParameter("mu2", 2, 3, 3, to_valarray(MODELPARAMETER(mu2_slha))));
   slha_soft_pars.push_back(TParameter("me2", 2, 3, 3, to_valarray(MODELPARAMETER(me2_slha))));
   slha_soft_pars.push_back(TParameter("MassB", 1, 1, 1, to_valarray(MODELPARAMETER(MassB))));
   slha_soft_pars.push_back(TParameter("MassWB", 1, 1, 1, to_valarray(MODELPARAMETER(MassWB))));
   slha_soft_pars.push_back(TParameter("MassG", 1, 1, 1, to_valarray(MODELPARAMETER(MassG))));
   slha_soft_pars.push_back(TParameter("Re(PhaseGlu)", 0, 1, 1, to_valarray(Re(MODELPARAMETER(PhaseGlu)))));
   slha_soft_pars.push_back(TParameter("Im(PhaseGlu)", 0, 1, 1, to_valarray(Im(MODELPARAMETER(PhaseGlu)))));
}

template <class T>
void CMSSM_semianalytic_slha_values_writer::extract_slha_soft_pars(const CMSSM_semianalytic_slha<T>& model)
{
   slha_soft_pars.clear();
   slha_soft_pars_inputs = model.get_input();
   slha_soft_pars_problems = model.get_problems();
   slha_soft_pars_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_soft_pars.push_back(TParameter("TYd", 1, 3, 3, to_valarray(MODELPARAMETER(TYd_slha))));
   slha_soft_pars.push_back(TParameter("TYe", 1, 3, 3, to_valarray(MODELPARAMETER(TYe_slha))));
   slha_soft_pars.push_back(TParameter("TYu", 1, 3, 3, to_valarray(MODELPARAMETER(TYu_slha))));
   slha_soft_pars.push_back(TParameter("BMu", 2, 1, 1, to_valarray(MODELPARAMETER(BMu))));
   slha_soft_pars.push_back(TParameter("mq2", 2, 3, 3, to_valarray(MODELPARAMETER(mq2_slha))));
   slha_soft_pars.push_back(TParameter("ml2", 2, 3, 3, to_valarray(MODELPARAMETER(ml2_slha))));
   slha_soft_pars.push_back(TParameter("mHd2", 2, 1, 1, to_valarray(MODELPARAMETER(mHd2))));
   slha_soft_pars.push_back(TParameter("mHu2", 2, 1, 1, to_valarray(MODELPARAMETER(mHu2))));
   slha_soft_pars.push_back(TParameter("md2", 2, 3, 3, to_valarray(MODELPARAMETER(md2_slha))));
   slha_soft_pars.push_back(TParameter("mu2", 2, 3, 3, to_valarray(MODELPARAMETER(mu2_slha))));
   slha_soft_pars.push_back(TParameter("me2", 2, 3, 3, to_valarray(MODELPARAMETER(me2_slha))));
   slha_soft_pars.push_back(TParameter("MassB", 1, 1, 1, to_valarray(MODELPARAMETER(MassB))));
   slha_soft_pars.push_back(TParameter("MassWB", 1, 1, 1, to_valarray(MODELPARAMETER(MassWB))));
   slha_soft_pars.push_back(TParameter("MassG", 1, 1, 1, to_valarray(MODELPARAMETER(MassG))));
   slha_soft_pars.push_back(TParameter("Re(PhaseGlu)", 0, 1, 1, to_valarray(Re(MODELPARAMETER(PhaseGlu)))));
   slha_soft_pars.push_back(TParameter("Im(PhaseGlu)", 0, 1, 1, to_valarray(Im(MODELPARAMETER(PhaseGlu)))));
}

template <class T>
void CMSSM_slha_values_writer::extract_slha_pole_mixings(const CMSSM_slha<T>& model)
{
   slha_pole_mixings.clear();
   slha_pole_mixings_inputs = model.get_input();
   slha_pole_mixings_problems = model.get_problems();

   slha_pole_mixings.push_back(TMixing("ZD", 6, true, to_valarray(PHYSICALSLHA(ZD))));
   slha_pole_mixings.push_back(TMixing("ZV", 3, true, to_valarray(PHYSICALSLHA(ZV))));
   slha_pole_mixings.push_back(TMixing("ZU", 6, true, to_valarray(PHYSICALSLHA(ZU))));
   slha_pole_mixings.push_back(TMixing("ZE", 6, true, to_valarray(PHYSICALSLHA(ZE))));
   slha_pole_mixings.push_back(TMixing("ZH", 2, true, to_valarray(PHYSICALSLHA(ZH))));
   slha_pole_mixings.push_back(TMixing("ZA", 2, true, to_valarray(PHYSICALSLHA(ZA))));
   slha_pole_mixings.push_back(TMixing("ZP", 2, true, to_valarray(PHYSICALSLHA(ZP))));
   slha_pole_mixings.push_back(TMixing("ZN", 4, false, to_valarray(PHYSICALSLHA(ZN))));
   slha_pole_mixings.push_back(TMixing("UM", 2, false, to_valarray(PHYSICALSLHA(UM))));
   slha_pole_mixings.push_back(TMixing("UP", 2, false, to_valarray(PHYSICALSLHA(UP))));
   slha_pole_mixings.push_back(TMixing("ZEL", 3, false, to_valarray(PHYSICALSLHA(ZEL))));
   slha_pole_mixings.push_back(TMixing("ZER", 3, false, to_valarray(PHYSICALSLHA(ZER))));
   slha_pole_mixings.push_back(TMixing("ZDL", 3, false, to_valarray(PHYSICALSLHA(ZDL))));
   slha_pole_mixings.push_back(TMixing("ZDR", 3, false, to_valarray(PHYSICALSLHA(ZDR))));
   slha_pole_mixings.push_back(TMixing("ZUL", 3, false, to_valarray(PHYSICALSLHA(ZUL))));
   slha_pole_mixings.push_back(TMixing("ZUR", 3, false, to_valarray(PHYSICALSLHA(ZUR))));
}

template <class T>
void CMSSM_semianalytic_slha_values_writer::extract_slha_pole_mixings(const CMSSM_semianalytic_slha<T>& model)
{
   slha_pole_mixings.clear();
   slha_pole_mixings_inputs = model.get_input();
   slha_pole_mixings_problems = model.get_problems();
   slha_pole_mixings_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_pole_mixings.push_back(TMixing("ZD", 6, true, to_valarray(PHYSICALSLHA(ZD))));
   slha_pole_mixings.push_back(TMixing("ZV", 3, true, to_valarray(PHYSICALSLHA(ZV))));
   slha_pole_mixings.push_back(TMixing("ZU", 6, true, to_valarray(PHYSICALSLHA(ZU))));
   slha_pole_mixings.push_back(TMixing("ZE", 6, true, to_valarray(PHYSICALSLHA(ZE))));
   slha_pole_mixings.push_back(TMixing("ZH", 2, true, to_valarray(PHYSICALSLHA(ZH))));
   slha_pole_mixings.push_back(TMixing("ZA", 2, true, to_valarray(PHYSICALSLHA(ZA))));
   slha_pole_mixings.push_back(TMixing("ZP", 2, true, to_valarray(PHYSICALSLHA(ZP))));
   slha_pole_mixings.push_back(TMixing("ZN", 4, false, to_valarray(PHYSICALSLHA(ZN))));
   slha_pole_mixings.push_back(TMixing("UM", 2, false, to_valarray(PHYSICALSLHA(UM))));
   slha_pole_mixings.push_back(TMixing("UP", 2, false, to_valarray(PHYSICALSLHA(UP))));
   slha_pole_mixings.push_back(TMixing("ZEL", 3, false, to_valarray(PHYSICALSLHA(ZEL))));
   slha_pole_mixings.push_back(TMixing("ZER", 3, false, to_valarray(PHYSICALSLHA(ZER))));
   slha_pole_mixings.push_back(TMixing("ZDL", 3, false, to_valarray(PHYSICALSLHA(ZDL))));
   slha_pole_mixings.push_back(TMixing("ZDR", 3, false, to_valarray(PHYSICALSLHA(ZDR))));
   slha_pole_mixings.push_back(TMixing("ZUL", 3, false, to_valarray(PHYSICALSLHA(ZUL))));
   slha_pole_mixings.push_back(TMixing("ZUR", 3, false, to_valarray(PHYSICALSLHA(ZUR))));
}

template <class T>
void CMSSM_slha_values_writer::extract_slha_running_mixings(const CMSSM_slha<T>& model)
{
   slha_running_mixings.clear();
   slha_running_mixings_inputs = model.get_input();
   slha_running_mixings_problems = model.get_problems();

   slha_running_mixings.push_back(TMixing("DRbarZD", 6, true, to_valarray(DRBARSLHA(ZD))));
   slha_running_mixings.push_back(TMixing("DRbarZV", 3, true, to_valarray(DRBARSLHA(ZV))));
   slha_running_mixings.push_back(TMixing("DRbarZU", 6, true, to_valarray(DRBARSLHA(ZU))));
   slha_running_mixings.push_back(TMixing("DRbarZE", 6, true, to_valarray(DRBARSLHA(ZE))));
   slha_running_mixings.push_back(TMixing("DRbarZH", 2, true, to_valarray(DRBARSLHA(ZH))));
   slha_running_mixings.push_back(TMixing("DRbarZA", 2, true, to_valarray(DRBARSLHA(ZA))));
   slha_running_mixings.push_back(TMixing("DRbarZP", 2, true, to_valarray(DRBARSLHA(ZP))));
   slha_running_mixings.push_back(TMixing("DRbarZN", 4, false, to_valarray(DRBARSLHA(ZN))));
   slha_running_mixings.push_back(TMixing("DRbarUM", 2, false, to_valarray(DRBARSLHA(UM))));
   slha_running_mixings.push_back(TMixing("DRbarUP", 2, false, to_valarray(DRBARSLHA(UP))));
   slha_running_mixings.push_back(TMixing("DRbarZEL", 3, false, to_valarray(DRBARSLHA(ZEL))));
   slha_running_mixings.push_back(TMixing("DRbarZER", 3, false, to_valarray(DRBARSLHA(ZER))));
   slha_running_mixings.push_back(TMixing("DRbarZDL", 3, false, to_valarray(DRBARSLHA(ZDL))));
   slha_running_mixings.push_back(TMixing("DRbarZDR", 3, false, to_valarray(DRBARSLHA(ZDR))));
   slha_running_mixings.push_back(TMixing("DRbarZUL", 3, false, to_valarray(DRBARSLHA(ZUL))));
   slha_running_mixings.push_back(TMixing("DRbarZUR", 3, false, to_valarray(DRBARSLHA(ZUR))));
}

template <class T>
void CMSSM_semianalytic_slha_values_writer::extract_slha_running_mixings(const CMSSM_semianalytic_slha<T>& model)
{
   slha_running_mixings.clear();
   slha_running_mixings_inputs = model.get_input();
   slha_running_mixings_problems = model.get_problems();
   slha_running_mixings_m0Sqr = model.get_ewsb_output_parameter(0);

   slha_running_mixings.push_back(TMixing("DRbarZD", 6, true, to_valarray(DRBARSLHA(ZD))));
   slha_running_mixings.push_back(TMixing("DRbarZV", 3, true, to_valarray(DRBARSLHA(ZV))));
   slha_running_mixings.push_back(TMixing("DRbarZU", 6, true, to_valarray(DRBARSLHA(ZU))));
   slha_running_mixings.push_back(TMixing("DRbarZE", 6, true, to_valarray(DRBARSLHA(ZE))));
   slha_running_mixings.push_back(TMixing("DRbarZH", 2, true, to_valarray(DRBARSLHA(ZH))));
   slha_running_mixings.push_back(TMixing("DRbarZA", 2, true, to_valarray(DRBARSLHA(ZA))));
   slha_running_mixings.push_back(TMixing("DRbarZP", 2, true, to_valarray(DRBARSLHA(ZP))));
   slha_running_mixings.push_back(TMixing("DRbarZN", 4, false, to_valarray(DRBARSLHA(ZN))));
   slha_running_mixings.push_back(TMixing("DRbarUM", 2, false, to_valarray(DRBARSLHA(UM))));
   slha_running_mixings.push_back(TMixing("DRbarUP", 2, false, to_valarray(DRBARSLHA(UP))));
   slha_running_mixings.push_back(TMixing("DRbarZEL", 3, false, to_valarray(DRBARSLHA(ZEL))));
   slha_running_mixings.push_back(TMixing("DRbarZER", 3, false, to_valarray(DRBARSLHA(ZER))));
   slha_running_mixings.push_back(TMixing("DRbarZDL", 3, false, to_valarray(DRBARSLHA(ZDL))));
   slha_running_mixings.push_back(TMixing("DRbarZDR", 3, false, to_valarray(DRBARSLHA(ZDR))));
   slha_running_mixings.push_back(TMixing("DRbarZUL", 3, false, to_valarray(DRBARSLHA(ZUL))));
   slha_running_mixings.push_back(TMixing("DRbarZUR", 3, false, to_valarray(DRBARSLHA(ZUR))));
}

} // namespace flexiblesusy

#endif
