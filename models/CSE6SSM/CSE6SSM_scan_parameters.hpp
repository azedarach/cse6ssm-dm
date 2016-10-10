// ====================================================================
// Helper class containing details of scans
// ====================================================================

#ifndef CSE6SSM_SCAN_PARAMETERS_H
#define CSE6SSM_SCAN_PARAMETERS_H

#include "CSE6SSM_two_scale_input_parameters.hpp"

#include "scan.hpp"

#include <chrono>

namespace flexiblesusy {

class CSE6SSM_scan_parameters {
public:
   CSE6SSM_scan_parameters();

   CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_, std::size_t m0_npts_,
                           double m12_lower_, double m12_upper_, std::size_t m12_npts_,
                           double TanBeta_lower_, double TanBeta_upper_, std::size_t TanBeta_npts_,
                           double SignLambdax_lower_, double SignLambdax_upper_, std::size_t SignLambdax_npts_,
                           double Azero_lower_, double Azero_upper_, std::size_t Azero_npts_,
                           double Kappa_lower_, double Kappa_upper_, std::size_t Kappa_npts_,
                           double Lambda12_lower_, double Lambda12_upper_, std::size_t Lambda12_npts_);

   CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_,
                           double m12_lower_, double m12_upper_,
                           double TanBeta_lower_, double TanBeta_upper_, 
                           double SignLambdax_lower_, double SignLambdax_upper_,
                           double Azero_lower_, double Azero_upper_, 
                           double Kappa_lower_, double Kappa_upper_,
                           double Lambda12_lower_, double Lambda12_upper_, std::size_t total_npts_);

   CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_, std::size_t m0_npts_,
                           double m12_lower_, double m12_upper_, std::size_t m12_npts_,
                           double TanBeta_lower_, double TanBeta_upper_, std::size_t TanBeta_npts_,
                           double SignLambdax_lower_, double SignLambdax_upper_, std::size_t SignLambdax_npts_,
                           double Azero_lower_, double Azero_upper_, std::size_t Azero_npts_,
                           double Kappa_lower_, double Kappa_upper_, std::size_t Kappa_npts_,
                           double Lambda12_lower_, double Lambda12_upper_, std::size_t Lambda12_npts_,
                           double output_scale_);

   CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_,
                           double m12_lower_, double m12_upper_,
                           double TanBeta_lower_, double TanBeta_upper_, 
                           double SignLambdax_lower_, double SignLambdax_upper_,
                           double Azero_lower_, double Azero_upper_, 
                           double Kappa_lower_, double Kappa_upper_, 
                           double Lambda12_lower_, double Lambda12_upper_, std::size_t total_npts_,
                           double output_scale_);

   ~CSE6SSM_scan_parameters() {}

   double get_m0_lower() const { return m0_lower; }
   double get_m0_upper() const { return m0_upper; }
   std::size_t get_m0_npts() const { return m0_npts; }
   double get_m0_incr() const { return m0_incr; }
   double get_m12_lower() const { return m12_lower; }
   double get_m12_upper() const { return m12_upper; }
   std::size_t get_m12_npts() const { return m12_npts; }
   double get_m12_incr() const { return m12_incr; }
   double get_TanBeta_lower() const { return TanBeta_lower; }
   double get_TanBeta_upper() const { return TanBeta_upper; }
   std::size_t get_TanBeta_npts() const { return TanBeta_npts; }
   double get_TanBeta_incr() const { return TanBeta_incr; }
   int get_SignLambdax_lower() const { return SignLambdax_lower; }
   int get_SignLambdax_upper() const { return SignLambdax_upper; }
   std::size_t get_SignLambdax_npts() const { return SignLambdax_npts; }
   int get_SignLambdax_incr() const { return SignLambdax_incr; }
   double get_Azero_lower() const { return Azero_lower; }
   double get_Azero_upper() const { return Azero_upper; }
   std::size_t get_Azero_npts() const { return Azero_npts; }
   double get_Azero_incr() const { return Azero_incr; }
   double get_Kappa_lower() const { return Kappa_lower; }
   double get_Kappa_upper() const { return Kappa_upper; }
   std::size_t get_Kappa_npts() const { return Kappa_npts; }
   double get_Kappa_incr() const { return Kappa_incr; }
   double get_Lambda12_lower() const { return Lambda12_lower; }
   double get_Lambda12_upper() const { return Lambda12_upper; }
   std::size_t get_Lambda12_npts() const { return Lambda12_npts; }
   double get_Lambda12_incr() const { return Lambda12_incr; }
   double get_total_npts() const { return total_npts; }
   bool get_is_grid_scan() const { return is_grid_scan; }
   double get_output_scale() const { return output_scale; }

   double get_random_m0(std::minstd_rand &);
   double get_random_m12(std::minstd_rand &);
   double get_random_TanBeta(std::minstd_rand &);
   double get_random_SignLambdax(std::minstd_rand &);
   double get_random_Azero(std::minstd_rand &);
   double get_random_Kappa(std::minstd_rand &);
   double get_random_Lambda12(std::minstd_rand &);

private:
   double m0_lower;
   double m0_upper;
   std::size_t m0_npts;
   double m0_incr;
   double m12_lower;
   double m12_upper;
   std::size_t m12_npts;
   double m12_incr;
   double TanBeta_lower;
   double TanBeta_upper;
   std::size_t TanBeta_npts;
   double TanBeta_incr;
   int SignLambdax_lower;
   int SignLambdax_upper;
   std::size_t SignLambdax_npts;
   int SignLambdax_incr;
   double Azero_lower;
   double Azero_upper;
   std::size_t Azero_npts;
   double Azero_incr;
   double Kappa_lower;
   double Kappa_upper;
   std::size_t Kappa_npts;
   double Kappa_incr;
   double Lambda12_lower;
   double Lambda12_upper;
   std::size_t Lambda12_npts;
   double Lambda12_incr;
   std::size_t total_npts;
   bool is_grid_scan;
   double output_scale;

   std::uniform_real_distribution<double> distribution;

};

} // namespace flexiblesusy

#endif
