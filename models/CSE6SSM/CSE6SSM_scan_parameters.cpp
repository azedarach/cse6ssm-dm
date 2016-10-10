#include "CSE6SSM_scan_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

   CSE6SSM_scan_parameters::CSE6SSM_scan_parameters()
      : m0_lower(0.), m0_upper(3000.), m0_npts(10)
      , m12_lower(0.), m12_upper(3000.), m12_npts(10)
      , TanBeta_lower(10.), TanBeta_upper(10.), TanBeta_npts(1)
      , SignLambdax_lower(-1), SignLambdax_upper(1), SignLambdax_npts(2)
      , Azero_lower(-1000.), Azero_upper(1000.), Azero_npts(5)
      , Kappa_lower(0.2), Kappa_upper(0.2), Kappa_npts(1)
      , Lambda12_lower(0.5), Lambda12_upper(0.5), Lambda12_npts(1)
      , is_grid_scan(true), output_scale(-1.)
   {
      m0_incr = (m0_npts > 1 ? (m0_upper - m0_lower) / (m0_npts - 1.) : 0.);
      m12_incr = (m12_npts > 1 ? (m12_upper - m12_lower) / (m12_npts - 1.) : 0.);
      TanBeta_incr = (TanBeta_npts > 1 ? (TanBeta_upper - TanBeta_lower) / (TanBeta_npts - 1.) : 0.);
      SignLambdax_incr = (SignLambdax_npts > 1 ? (SignLambdax_upper - SignLambdax_lower) : 0);
      Azero_incr = (Azero_npts > 1 ? (Azero_upper - Azero_lower) / (Azero_npts - 1.) : 0.);
      Kappa_incr = (Kappa_npts > 1 ? (Kappa_upper - Kappa_lower) / (Kappa_npts - 1.) : 0.);
      Lambda12_incr = (Lambda12_npts > 1 ? (Lambda12_upper - Lambda12_lower) / (Lambda12_npts - 1.) : 0.);

      total_npts = m0_npts * m12_npts * TanBeta_npts * SignLambdax_npts * Azero_npts * Kappa_npts * Lambda12_npts;

      distribution = std::uniform_real_distribution<double>(0.,1.);
   }

   CSE6SSM_scan_parameters::CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_, std::size_t m0_npts_,
                                                    double m12_lower_, double m12_upper_, std::size_t m12_npts_,
                                                    double TanBeta_lower_, double TanBeta_upper_, std::size_t TanBeta_npts_,
                                                    double SignLambdax_lower_, double SignLambdax_upper_, 
                                                    std::size_t SignLambdax_npts_,
                                                    double Azero_lower_, double Azero_upper_, std::size_t Azero_npts_,
                                                    double Kappa_lower_, double Kappa_upper_, std::size_t Kappa_npts_,
                                                    double Lambda12_lower_, double Lambda12_upper_, std::size_t Lambda12_npts_)
   : m0_lower(m0_lower_), m0_upper(m0_upper_), m0_npts(m0_npts_)
   , m12_lower(m12_lower_), m12_upper(m12_upper_), m12_npts(m12_npts_)
   , TanBeta_lower(TanBeta_lower_), TanBeta_upper(TanBeta_upper_), TanBeta_npts(TanBeta_npts_)
   , Azero_lower(Azero_lower_), Azero_upper(Azero_upper_), Azero_npts(Azero_npts_) 
   , Kappa_lower(Kappa_lower_), Kappa_upper(Kappa_upper_), Kappa_npts(Kappa_npts_)
   , Lambda12_lower(Lambda12_lower_), Lambda12_upper(Lambda12_upper_), Lambda12_npts(Lambda12_npts_)
   , is_grid_scan(true), output_scale(-1.) 
   {
      SignLambdax_lower = Sign(SignLambdax_lower_);
      SignLambdax_upper = Sign(SignLambdax_upper_);
      
      if (SignLambdax_lower == SignLambdax_upper) {
         SignLambdax_npts = 1;
      } else {
         SignLambdax_npts = (SignLambdax_npts_ > 1 ? 2 : 1);
      }
      
      m0_incr = (m0_npts > 1 ? (m0_upper - m0_lower) / (m0_npts - 1.) : 0.);
      m12_incr = (m12_npts > 1 ? (m12_upper - m12_lower) / (m12_npts - 1.) : 0.);
      TanBeta_incr = (TanBeta_npts > 1 ? (TanBeta_upper - TanBeta_lower) / (TanBeta_npts - 1.) : 0.);
      SignLambdax_incr = (SignLambdax_npts > 1 ? (SignLambdax_upper - SignLambdax_lower) : 0);
      Azero_incr = (Azero_npts > 1 ? (Azero_upper - Azero_lower) / (Azero_npts - 1.) : 0.);
      Kappa_incr = (Kappa_npts > 1 ? (Kappa_upper - Kappa_lower) / (Kappa_npts - 1.) : 0.);
      Lambda12_incr = (Lambda12_npts > 1 ? (Lambda12_upper - Lambda12_lower) / (Lambda12_npts - 1.) : 0.);      

      total_npts = m0_npts * m12_npts * TanBeta_npts * SignLambdax_npts * Azero_npts * Kappa_npts * Lambda12_npts;

      distribution = std::uniform_real_distribution<double>(0.,1.);
   }
   
   CSE6SSM_scan_parameters::CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_,
                                                    double m12_lower_, double m12_upper_,
                                                    double TanBeta_lower_, double TanBeta_upper_, 
                                                    double SignLambdax_lower_, double SignLambdax_upper_,
                                                    double Azero_lower_, double Azero_upper_, 
                                                    double Kappa_lower_, double Kappa_upper_,
                                                    double Lambda12_lower_, double Lambda12_upper_, 
                                                    std::size_t total_npts_)
      : m0_lower(m0_lower_), m0_upper(m0_upper_), m0_incr(0.)
      , m12_lower(m12_lower_), m12_upper(m12_upper_), m12_incr(0.)
      , TanBeta_lower(TanBeta_lower_), TanBeta_upper(TanBeta_upper_), TanBeta_incr(0.)
      , Azero_lower(Azero_lower_), Azero_upper(Azero_upper_), Azero_incr(0.)
      , Kappa_lower(Kappa_lower_), Kappa_upper(Kappa_upper_), Kappa_incr(0.)
      , Lambda12_lower(Lambda12_lower_), Lambda12_upper(Lambda12_upper_), Lambda12_incr(0.)
      , total_npts(total_npts_)
      , is_grid_scan(false), output_scale(-1.)
   {
      SignLambdax_lower = Sign(SignLambdax_lower_);
      SignLambdax_upper = Sign(SignLambdax_upper_);
      SignLambdax_incr = 0;
      
      // in this case the number of points assigned to each 
      // is irrelevant
      m0_npts = 1;
      m12_npts = 1;
      TanBeta_npts = 1;
      SignLambdax_npts = 1;
      Azero_npts = total_npts;
      Kappa_npts = 1;
      Lambda12_npts = 1;

      distribution = std::uniform_real_distribution<double>(0.,1.);
   }
   
   CSE6SSM_scan_parameters::CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_, std::size_t m0_npts_,
                                                    double m12_lower_, double m12_upper_, std::size_t m12_npts_,
                                                    double TanBeta_lower_, double TanBeta_upper_, std::size_t TanBeta_npts_,
                                                    double SignLambdax_lower_, double SignLambdax_upper_, 
                                                    std::size_t SignLambdax_npts_,
                                                    double Azero_lower_, double Azero_upper_, std::size_t Azero_npts_,
                                                    double Kappa_lower_, double Kappa_upper_, std::size_t Kappa_npts_,
                                                    double Lambda12_lower_, double Lambda12_upper_, std::size_t Lambda12_npts_,
                                                    double output_scale_)
   : m0_lower(m0_lower_), m0_upper(m0_upper_), m0_npts(m0_npts_)
   , m12_lower(m12_lower_), m12_upper(m12_upper_), m12_npts(m12_npts_)
   , TanBeta_lower(TanBeta_lower_), TanBeta_upper(TanBeta_upper_), TanBeta_npts(TanBeta_npts_)
   , Azero_lower(Azero_lower_), Azero_upper(Azero_upper_), Azero_npts(Azero_npts_)
   , Kappa_lower(Kappa_lower_), Kappa_upper(Kappa_upper_), Kappa_npts(Kappa_npts_)
   , Lambda12_lower(Lambda12_lower_), Lambda12_upper(Lambda12_upper_), Lambda12_npts(Lambda12_npts_)
   ,is_grid_scan(true), output_scale(output_scale_) 
   {
      SignLambdax_lower = Sign(SignLambdax_lower_);
      SignLambdax_upper = Sign(SignLambdax_upper_);
      
      if (SignLambdax_lower == SignLambdax_upper) {
         SignLambdax_npts = 1;
      } else {
         SignLambdax_npts = (SignLambdax_npts_ > 1 ? 2 : 1);
      }
      
      m0_incr = (m0_npts > 1 ? (m0_upper - m0_lower) / (m0_npts - 1.) : 0.);
      m12_incr = (m12_npts > 1 ? (m12_upper - m12_lower) / (m12_npts - 1.) : 0.);
      TanBeta_incr = (TanBeta_npts > 1 ? (TanBeta_upper - TanBeta_lower) / (TanBeta_npts - 1.) : 0.);
      SignLambdax_incr = (SignLambdax_npts > 1 ? (SignLambdax_upper - SignLambdax_lower) : 0);
      Azero_incr = (Azero_npts > 1 ? (Azero_upper - Azero_lower) / (Azero_npts - 1.) : 0.);
      Kappa_incr = (Kappa_npts > 1 ? (Kappa_upper - Kappa_lower) / (Kappa_npts - 1.) : 0.);
      Lambda12_incr = (Lambda12_npts > 1 ? (Lambda12_upper - Lambda12_lower) / (Lambda12_npts - 1.) : 0.);

      total_npts = m0_npts * m12_npts * TanBeta_npts * SignLambdax_npts * Azero_npts * Kappa_npts * Lambda12_npts;

      distribution = std::uniform_real_distribution<double>(0.,1.);
   }
   
   CSE6SSM_scan_parameters::CSE6SSM_scan_parameters(double m0_lower_, double m0_upper_,
                                                    double m12_lower_, double m12_upper_,
                                                    double TanBeta_lower_, double TanBeta_upper_, 
                                                    double SignLambdax_lower_, double SignLambdax_upper_,
                                                    double Azero_lower_, double Azero_upper_, 
                                                    double Kappa_lower_, double Kappa_upper_,
                                                    double Lambda12_lower_, double Lambda12_upper_,
                                                    std::size_t total_npts_, double output_scale_)
      : m0_lower(m0_lower_), m0_upper(m0_upper_), m0_incr(0.)
      , m12_lower(m12_lower_), m12_upper(m12_upper_), m12_incr(0.)
      , TanBeta_lower(TanBeta_lower_), TanBeta_upper(TanBeta_upper_), TanBeta_incr(0.)
      , Azero_lower(Azero_lower_), Azero_upper(Azero_upper_), Azero_incr(0.)
      , Kappa_lower(Kappa_lower_), Kappa_upper(Kappa_upper_), Kappa_incr(0.)
      , Lambda12_lower(Lambda12_lower_), Lambda12_upper(Lambda12_upper_), Lambda12_incr(0.)
      , total_npts(total_npts_)
      , is_grid_scan(false), output_scale(output_scale_)
   {
      SignLambdax_lower = Sign(SignLambdax_lower_);
      SignLambdax_upper = Sign(SignLambdax_upper_);
      SignLambdax_incr = 0;
      
      // in this case the number of points assigned to each 
      // is irrelevant
      m0_npts = 1;
      m12_npts = 1;
      TanBeta_npts = 1;
      SignLambdax_npts = 1;
      Azero_npts = total_npts;
      Kappa_npts = 1;
      Lambda12_npts = 1;

      distribution = std::uniform_real_distribution<double>(0.,1.);
   }

   double CSE6SSM_scan_parameters::get_random_m0(std::minstd_rand & generator)
   {
      return m0_lower + (m0_upper - m0_lower) * distribution(generator);
   }

   double CSE6SSM_scan_parameters::get_random_m12(std::minstd_rand & generator)
   {
      return m12_lower + (m12_upper - m12_lower) * distribution(generator);
   }

   double CSE6SSM_scan_parameters::get_random_TanBeta(std::minstd_rand & generator)
   {
      return TanBeta_lower + (TanBeta_upper - TanBeta_lower) * distribution(generator);
   }

   double CSE6SSM_scan_parameters::get_random_SignLambdax(std::minstd_rand & generator)
   {
      return Sign(SignLambdax_lower + (SignLambdax_upper - SignLambdax_lower) * distribution(generator));
   }

   double CSE6SSM_scan_parameters::get_random_Azero(std::minstd_rand & generator)
   {
      return Azero_lower + (Azero_upper - Azero_lower) * distribution(generator);
   }

   double CSE6SSM_scan_parameters::get_random_Kappa(std::minstd_rand & generator)
   {
      return Kappa_lower + (Kappa_upper - Kappa_lower) * distribution(generator);
   }

   double CSE6SSM_scan_parameters::get_random_Lambda12(std::minstd_rand & generator)
   {
      return Lambda12_lower + (Lambda12_upper - Lambda12_lower) * distribution(generator);
   }

} // namespace flexiblesusy
