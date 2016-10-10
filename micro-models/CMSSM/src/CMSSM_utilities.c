/**
 * @file CMSSM_utilities.c
 * @brief contains implementations of model-specific functions
 */

#include "CMSSM_utilities.h"

#include "micromegas.h"
#include "micromegas_aux.h"

#include <stdio.h>

void print_DM_candidate_info(char* candidate)
{
   if (candidate) {
      int spin_times_two;
      int charge_times_three;
      int colour_rep;

      qNumbers(candidate, &spin_times_two, &charge_times_three, &colour_rep);

      if (spin_times_two == 0)
         printf("Dark matter candidate is '%s', spin = 0, mass = %.2E\n",
                candidate, Mcdm1);
      else
         printf("Dark matter candidate is '%s', spin = %d/2, mass = %.2E\n",
                candidate, spin_times_two, Mcdm1);

      if (charge_times_three != 0)
         printf("\tNote: candidate has electric charge %d/3\n",
                charge_times_three);

      if (colour_rep != 1)
         printf("\tNote: candidate has colour charge\n");
   }
}

int has_two_DM_candidates(void)
{
   return CDM1 && CDM2;
}

/* Note no freeze-out approximation is available for 2 DM candidates */
double calculate_relic_density(int fast, double thresh, int freeze_out)
{
   double omega = darkOmega2(fast, thresh);
   return omega;
}

double calculate_single_DM_relic_density(int fast, double thresh, double* xf,
                                         int freeze_out)
{
   double omega = 0.;
   if (freeze_out) {
      omega = darkOmegaFO(xf, fast, thresh);
   } else {
      omega = darkOmega(xf, fast, thresh);
   }
   return omega;
}

double calculate_annihilation_spectra(int key, double* gamma, double* eplus,
                                      double* pbar, double* nu_e,
                                      double* nu_mu, double* nu_tau,
                                      int* error)
{
   double sigma_v = calcSpectrum(key, gamma, eplus, pbar, nu_e, nu_mu,
                                 nu_tau, error);
   return sigma_v;
}

void get_photon_flux_table(double fi, double dfi, double sigma_v,
                           double* spectrum_gamma, double* flux_gamma)
{
   gammaFluxTab(fi, dfi, sigma_v, spectrum_gamma,  flux_gamma);
}

void get_positron_flux_table(double E_min, double sigma_v,
                             double* spectrum_eplus, double* flux_eplus)
{
   posiFluxTab(E_min, sigma_v, spectrum_eplus, flux_eplus);
}

void get_antiproton_flux_table(double E_min, double sigma_v,
                               double* spectrum_pbar, double* flux_pbar)
{
   pbarFluxTab(E_min, sigma_v, spectrum_pbar, flux_pbar);
}

int calculate_nucleon_amplitudes(char* candidate,
                                 double (*loop_order)(double,double,
                                                      double,double),
                                 double *pA0, double *pA5, double *nA0,
                                 double *nA5)
{
   if (!candidate) {
      return 1;
   }

   int error = nucleonAmplitudes(candidate, loop_order, pA0, pA5, nA0, nA5);

   return error;
}

static double get_cross_section_normalisation(double mdm, double mnuc)
{
   /* note conversion to pb */
   double norm = 4.0 * 3.8937966E8 * pow(mnuc * mdm / (mnuc + mdm), 2.) / M_PI;
   return norm;
}

double get_SI_cross_section(double mdm, double mp, double amp)
{
   double norm = get_cross_section_normalisation(mdm, mp);
   return norm * amp * amp;
}

double get_SD_cross_section(double mdm, double mp, double amp)
{
   double norm = get_cross_section_normalisation(mdm, mp);
   return 3.0 * norm * amp * amp;
}
