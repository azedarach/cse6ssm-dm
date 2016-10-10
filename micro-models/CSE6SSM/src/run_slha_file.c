/**
 * @file run_slha_file.c
 * @brief run an SLHA file generated by FlexibleSUSY in micrOMEGAs
 */

#include "micromegas.h"
#include "micromegas_aux.h"

#include "CSE6SSM_utilities.h"
#include "read_slha.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

void parse_options(int*, char***);
void print_usage(void);
void set_standard_model_inputs(void);

/* Flag set by '--verbose'. */
static int verbose_flag = 0;

/* Flag set by '--enable-relic-density'. */
static int relic_density_flag = 1;

/* Flag set by '--enable-direct-detection'. */
static int direct_det_flag = 0;

/* Flag set by '--enable-indirect-detection'. */
static int indirect_det_flag = 0;

/* Optional name of output file */
static char output_filename[1024] = "run_slha_file.out";

int main(int argc, char** argv)
{
   parse_options(&argc, (char***)&argv);

   if (optind == argc) {
      fprintf(stderr, "Usage: ./run_slha_file.x [OPTION]... FILE\n"
              "Try './run_slha_file.x --help' for more information.\n");
      return EXIT_FAILURE;
   }

   set_standard_model_inputs();

   /* Read SLHA file */
   int error_code = 0;
   error_code = read_slha_file(argv[optind]);

   if (error_code) {
      fprintf(stderr, "Error: could not read SLHA file!\n");
      return error_code;
   }

   error_code = sortOddParticles(0);
   if (error_code) {
      fprintf(stderr, "Error: could not identify dark matter candidates!\n");
      return error_code;
   }

   if (verbose_flag) {
      printf("###### Model information ######\n");
      if (CDM1)
         print_DM_candidate_info(CDM1);
      if (CDM2)
         print_DM_candidate_info(CDM2);
   }

   /* Observables to calculate */
   double relic_density = 0.;
   double E_test = 0.;
   double sigma_v = 0.;
   double photon_flux = 0.;
   double positron_flux = 0.;
   double antiproton_flux = 0.;
   double proton_SI_amp = 0.;
   double proton_SI_anti_amp = 0.;
   double proton_SD_amp = 0.;
   double proton_SD_anti_amp = 0.;
   double neutron_SI_amp = 0.;
   double neutron_SI_anti_amp = 0.;
   double neutron_SD_amp = 0.;
   double neutron_SD_anti_amp = 0.;
   double proton_SI_xsec = 0.;
   double proton_SI_anti_xsec = 0.;
   double proton_SD_xsec = 0.;
   double proton_SD_anti_xsec = 0.;
   double neutron_SI_xsec = 0.;
   double neutron_SI_anti_xsec = 0.;
   double neutron_SD_xsec = 0.;
   double neutron_SD_anti_xsec = 0.;

   if (relic_density_flag) {
      if (verbose_flag) {
         printf("###### Relic density calculation ######\n");
      }

      /* For the flags VZdecay and VWdecay, a value of
         - 0 indicates 3-body final states are not included
         - 1 indicates 3-body final states are included in annihilations
         - 2 indicates 3-body final states are also included in coannihilations
       */
      VZdecay = 2;
      VWdecay = 2;

      int fast = 1;
      double Beps = 1.0e-8;
      int freeze_out = 0;

      if (verbose_flag) {
         if (fast) {
            printf("Performing fast calculation\n");
         } else {
            printf("Performing high precision calculation\n");
         }

         printf("Only including processes with B > %.2E\n", Beps);
      }

      if (has_two_DM_candidates()) {
         relic_density = calculate_relic_density(fast, Beps, freeze_out);
      } else {
         double Xf;
         relic_density = calculate_single_DM_relic_density(fast, Beps, &Xf,
                                                           freeze_out);

         if (verbose_flag) {
            double cut = 0.01;
            printf("Obtained freeze-out temp Xf = %.2E\n", Xf);
            double sum = printChannels(Xf, cut, Beps, 1, stdout);

            printf("Result of printChannels = %.8E\n", sum);
            printf("Detailed channels:\n");

            double tN1tN1_to_u3U3 = oneChannel(Xf, Beps, "~N1", "~N1",
                                               "u3", "U3");
            double tN1tN1_to_WmWp = oneChannel(Xf, Beps, "~N1", "~N1",
                                               "Wm", "Wp");
            double tN1tN1_to_ZZ = oneChannel(Xf, Beps, "~N1", "~N1",
                                             "Z", "Z");
            double tN1tN1_to_h1Z = oneChannel(Xf, Beps, "~N1", "~N1",
                                              "h1", "Z");
            double tN1tC1_to_d3U3 = oneChannel(Xf, Beps, "~N1", "~C1",
                                               "d3", "U3");
            double tN1tC1_to_d2U2 = oneChannel(Xf, Beps, "~N1", "~C1",
                                               "d2", "U2");
            double tN1tC1_to_d1U1 = oneChannel(Xf, Beps, "~N1", "~C1",
                                               "d1", "U1");
            double tN1tC1_to_h1Wm = oneChannel(Xf, Beps, "~N1", "~C1",
                                               "h1", "Wm");
            double tN1tC1_to_Nu1e1 = oneChannel(Xf, Beps, "~N1", "~C1",
                                                "Nu1", "e1");
            double tN1tC1_to_Nu2e2 = oneChannel(Xf, Beps, "~N1", "~C1",
                                                "Nu2", "e2");
            double tN1tC1_to_Nu3e3 = oneChannel(Xf, Beps, "~N1", "~C1",
                                                "Nu3", "e3");
            double tN1tC1_to_ZWm = oneChannel(Xf, Beps, "~N1", "~C1",
                                              "Z", "Wm");
            double tN2tC1_to_d1U1 = oneChannel(Xf, Beps, "~N2", "~C1",
                                               "d1", "U1");
            double tN2tC1_to_d2U2 = oneChannel(Xf, Beps, "~N2", "~C1",
                                               "d2", "U2");
            double tN2tC1_to_d3U3 = oneChannel(Xf, Beps, "~N2", "~C1",
                                               "d3", "U3");
            double tN1tN2_to_d1D1 = oneChannel(Xf, Beps, "~N1", "~N2",
                                               "d1", "D1");
            double tN1tN2_to_d2D2 = oneChannel(Xf, Beps, "~N1", "~N2",
                                               "d2", "D2");
            double tN1tN2_to_d3D3 = oneChannel(Xf, Beps, "~N1", "~N2",
                                               "d3", "D3");
            double tC1tc1_to_WmWp = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "Wm", "Wp");
            double tN1tN2_to_u1U1 = oneChannel(Xf, Beps, "~N1", "~N2",
                                               "u1", "U1");
            double tN1tN2_to_u2U2 = oneChannel(Xf, Beps, "~N1", "~N2",
                                               "u2", "U2");
            double tN1tN2_to_u3U3 = oneChannel(Xf, Beps, "~N1", "~N2",
                                               "u3", "U3");
            double tC1tc1_to_u1U1 = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "u1", "U1");
            double tC1tc1_to_u2U2 = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "u2", "U2");
            double tC1tc1_to_u3U3 = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "u3", "U3");
            double tN2tC1_to_Nu1e1 = oneChannel(Xf, Beps, "~N2", "~C1",
                                                "Nu1", "e1");
            double tN2tC1_to_Nu2e2 = oneChannel(Xf, Beps, "~N2", "~C1",
                                                "Nu2", "e2");
            double tN2tC1_to_Nu3e3 = oneChannel(Xf, Beps, "~N2", "~C1",
                                                "Nu3", "e3");
            double tN1tC1_to_AWm = oneChannel(Xf, Beps, "~N1", "~C1",
                                              "A", "Wm");
            double tC1tc1_to_d1D1 = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "d1", "D1");
            double tC1tc1_to_d2D2 = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "d2", "D2");
            double tC1tc1_to_d3D3 = oneChannel(Xf, Beps, "~C1", "~c1",
                                               "d3", "D3");
            double tN1tN1_to_h1h1 = oneChannel(Xf, Beps, "~N1", "~N1",
                                               "h1", "h1");

            double total_contributions = tN1tN1_to_u3U3 + tN1tN1_to_WmWp
               + tN1tN1_to_ZZ + tN1tN1_to_h1Z + tN1tC1_to_d3U3
               + tN1tC1_to_d2U2 + tN1tC1_to_d1U1 + tN1tC1_to_h1Wm
               + tN1tC1_to_Nu1e1 + tN1tC1_to_Nu2e2 + tN1tC1_to_Nu3e3
               + tN1tC1_to_ZWm + tN2tC1_to_d1U1 + tN2tC1_to_d2U2
               + tN2tC1_to_d3U3 + tN1tN2_to_d1D1 + tN1tN2_to_d2D2
               + tN1tN2_to_d3D3 + tC1tc1_to_WmWp + tN1tN2_to_u1U1
               + tN1tN2_to_u2U2 + tN1tN2_to_u3U3 + tC1tc1_to_u1U1
               + tC1tc1_to_u2U2 + tC1tc1_to_u3U3 + tN2tC1_to_Nu1e1
               + tN2tC1_to_Nu2e2 + tN2tC1_to_Nu3e3 + tN1tC1_to_AWm
               + tC1tc1_to_d1D1 + tC1tc1_to_d2D2 + tC1tc1_to_d3D3
               + tN1tN1_to_h1h1;

            printf("~N1 ~N1 -> u3 U3 = %.8E\n", 100.0 * tN1tN1_to_u3U3);
            printf("~N1 ~N1 -> Wm Wp = %.8E\n", 100.0 * tN1tN1_to_WmWp);
            printf("~N1 ~N1 -> Z Z = %.8E\n", 100.0 * tN1tN1_to_ZZ);
            printf("~N1 ~N1 -> h1 Z = %.8E\n", 100.0 * tN1tN1_to_h1Z);
            printf("~N1 ~C1 -> d3 U3 = %.8E\n", 100.0 * tN1tC1_to_d3U3);
            printf("~N1 ~C1 -> d2 U2 = %.8E\n", 100.0 * tN1tC1_to_d2U2);
            printf("~N1 ~C1 -> d1 U1 = %.8E\n", 100.0 * tN1tC1_to_d1U1);
            printf("~N1 ~C1 -> h1 Wm = %.8E\n", 100.0 * tN1tC1_to_h1Wm);
            printf("~N1 ~C1 -> Nu1 e1 = %.8E\n", 100.0 * tN1tC1_to_Nu1e1);
            printf("~N1 ~C1 -> Nu2 e2 = %.8E\n", 100.0 * tN1tC1_to_Nu2e2);
            printf("~N1 ~C1 -> Nu3 e3 = %.8E\n", 100.0 * tN1tC1_to_Nu3e3);
            printf("~N1 ~C1 -> Z Wm = %.8E\n", 100.0 * tN1tC1_to_ZWm);
            printf("~N2 ~C1 -> d1 U1 = %.8E\n", 100.0 * tN2tC1_to_d1U1);
            printf("~N2 ~C1 -> d2 U2 = %.8E\n", 100.0 * tN2tC1_to_d2U2);
            printf("~N2 ~C1 -> d3 U3 = %.8E\n", 100.0 * tN2tC1_to_d3U3);
            printf("~N1 ~N2 -> d1 D1 = %.8E\n", 100.0 * tN1tN2_to_d1D1);
            printf("~N1 ~N2 -> d2 D2 = %.8E\n", 100.0 * tN1tN2_to_d2D2);
            printf("~N1 ~N2 -> d3 D3 = %.8E\n", 100.0 * tN1tN2_to_d3D3);
            printf("~C1 ~c1 -> Wm Wp = %.8E\n", 100.0 * tC1tc1_to_WmWp);
            printf("~N1 ~N2 -> u1 U1 = %.8E\n", 100.0 * tN1tN2_to_u1U1);
            printf("~N1 ~N2 -> u2 U2 = %.8E\n", 100.0 * tN1tN2_to_u2U2);
            printf("~N1 ~N2 -> u3 U3 = %.8E\n", 100.0 * tN1tN2_to_u3U3);
            printf("~C1 ~c1 -> u1 U1 = %.8E\n", 100.0 * tC1tc1_to_u1U1);
            printf("~C1 ~c1 -> u2 U2 = %.8E\n", 100.0 * tC1tc1_to_u2U2);
            printf("~C1 ~c1 -> u3 U3 = %.8E\n", 100.0 * tC1tc1_to_u3U3);
            printf("~N2 ~C1 -> Nu1 e1 = %.8E\n", 100.0 * tN2tC1_to_Nu1e1);
            printf("~N2 ~C1 -> Nu2 e2 = %.8E\n", 100.0 * tN2tC1_to_Nu2e2);
            printf("~N2 ~C1 -> Nu3 e3 = %.8E\n", 100.0 * tN2tC1_to_Nu3e3);
            printf("~N1 ~C1 -> A Wm = %.8E\n", 100.0 * tN1tC1_to_AWm);
            printf("~C1 ~c1 -> d1 D1 = %.8E\n", 100.0 * tC1tc1_to_d1D1);
            printf("~C1 ~c1 -> d2 D2 = %.8E\n", 100.0 * tC1tc1_to_d2D2);
            printf("~C1 ~c1 -> d3 D3 = %.8E\n", 100.0 * tC1tc1_to_d3D3);
            printf("~N1 ~N1 -> h1 h1 = %.8E\n", 100.0 * tN1tN1_to_h1h1);

            printf("Cumulative channels:\n");
            printf("~N1 ~C1 -> di Ui = %.8E\n",
                   100.0 * (tN1tC1_to_d1U1 + tN1tC1_to_d2U2 + tN1tC1_to_d3U3)); 
            printf("~N1 ~C1 -> Nui ei = %.8E\n",
                   100.0 * (tN1tC1_to_Nu1e1 + tN1tC1_to_Nu2e2
                            + tN1tC1_to_Nu3e3)); 
            printf("~N2 ~C1 -> di Ui = %.8E\n",
                   100.0 * (tN2tC1_to_d1U1 + tN2tC1_to_d2U2 + tN2tC1_to_d3U3)); 
            printf("~N2 ~C1 -> Nui ei = %.8E\n",
                   100.0 * (tN2tC1_to_Nu1e1 + tN2tC1_to_Nu2e2
                            + tN2tC1_to_Nu3e3));
            printf("~N1 ~N2 -> di Di = %.8E\n",
                   100.0 * (tN1tN2_to_d1D1 + tN1tN2_to_d2D2 + tN1tN2_to_d3D3)); 
            printf("~N1 ~N2 -> ui Ui = %.8E\n",
                   100.0 * (tN1tN2_to_u1U1 + tN1tN2_to_u2U2 + tN1tN2_to_u3U3)); 
            printf("~C1 ~c1 -> ui Ui = %.8E\n",
                   100.0 * (tC1tc1_to_u1U1 + tC1tc1_to_u2U2 + tC1tc1_to_u3U3)); 
            printf("~C1 ~c1 -> di Di = %.8E\n",
                   100.0 * (tC1tc1_to_d1D1 + tC1tc1_to_d2D2 + tC1tc1_to_d3D3));
            printf("Total = %.8E\n", total_contributions);

            printf("Compared to freeze-out approximation:\n");
            double Xf_fo;
            double relic_density_fo =
               calculate_single_DM_relic_density(fast, Beps, &Xf_fo, 1);
            printf("Omega (freeze-out approx.) = %.8E\n", relic_density_fo);
            printf("Xf (freeze-out approx.) = %.2E\n", Xf_fo);
         }
      }

      if (verbose_flag) {
         printf("Relic density = %.8E\n", relic_density);
      }
   }

   if (indirect_det_flag) {
      if (verbose_flag) {
         printf("###### Indirect detection calculation ######\n");
      }

      int key = 3;
      if (verbose_flag) {
         key += 4;
      }
      double spectrum_gamma[NZ];
      double spectrum_eplus[NZ];
      double spectrum_pbar[NZ];
      double* spectrum_nu_e = NULL;
      double* spectrum_nu_mu = NULL;
      double* spectrum_nu_tau = NULL;

      sigma_v = calculate_annihilation_spectra(key, spectrum_gamma,
                                               spectrum_eplus,
                                               spectrum_pbar,
                                               spectrum_nu_e,
                                               spectrum_nu_mu,
                                               spectrum_nu_tau,
                                               &error_code);

      if (error_code) {
         fprintf(stderr, "Error: could not calculate annihilation spectra! "
                 "(error code %d)\n", error_code);
         return error_code;
      }

      double fi = 0.1;
      double dfi = 0.05;
      double flux_gamma[NZ];
      double flux_eplus[NZ];
      double flux_pbar[NZ];
      double E_min = 1.0;
      E_test = 0.5 * Mcdm;

      if (verbose_flag) {
         printf("Calculating photon flux for angle of sight f = %.2E rad\n"
                "and spherical region described by cone with angle "
                "%.2E rad\n", fi, dfi);
      }

      get_photon_flux_table(fi, dfi, sigma_v, spectrum_gamma, flux_gamma);
      photon_flux = SpectdNdE(E_test, flux_gamma);

      get_positron_flux_table(E_min, sigma_v, spectrum_eplus, flux_eplus);
      positron_flux = SpectdNdE(E_test, flux_eplus);

      get_antiproton_flux_table(E_min, sigma_v, spectrum_pbar, flux_pbar);
      antiproton_flux = SpectdNdE(E_test, flux_pbar);

      if (verbose_flag) {
         printf("Annihilation cross section = %.8E cm^3/s\n", sigma_v);
         printf("Photon flux = %.8E (cm^2 s GeV)^{-1} for E = %.4E GeV\n",
                photon_flux, E_test);
         printf("Positron flux = %.8E (cm^2 s GeV)^{-1} for E = %.4E GeV\n",
                positron_flux, E_test);
         printf("Antiproton flux = %.8E (cm^2 s GeV)^{-1} for E = %.4E GeV\n",
                antiproton_flux, E_test);
      }
   }

   if (direct_det_flag) {
      if (verbose_flag) {
         printf("###### Direct detection calculation ######\n");
      }

      double pA0[2];
      double pA5[2];
      double nA0[2];
      double nA5[2];

      double nucleon_mass = 0.939;

      error_code = calculate_nucleon_amplitudes(CDM1, FeScLoop, pA0, pA5,
                                                nA0, nA5);

      if (error_code) {
         fprintf(stderr, "Error: could not calculate scattering amplitudes! "
                 "(error code %d)\n", error_code);
         return error_code;
      }

      proton_SI_amp = pA0[0];
      proton_SI_anti_amp = pA0[1];
      proton_SD_amp = pA5[0];
      proton_SD_anti_amp = pA5[1];
      neutron_SI_amp = nA0[0];
      neutron_SI_anti_amp = nA0[1];
      neutron_SD_amp = nA5[0];
      neutron_SD_anti_amp = nA5[1];

      if (verbose_flag) {
         printf("CDM[antiCDM]-nucleon amplitudes for %s \n", CDM1);
         printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",
                proton_SI_amp, proton_SI_anti_amp,
                proton_SD_amp, proton_SD_anti_amp);
         printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",
                neutron_SI_amp, neutron_SI_anti_amp,
                neutron_SD_amp, neutron_SD_anti_amp);
      }

      double pb_to_cm2 = 1.0E-36;

      proton_SI_xsec = pb_to_cm2 *
         get_SI_cross_section(Mcdm, nucleon_mass, proton_SI_amp);
      proton_SI_anti_xsec = pb_to_cm2*
         get_SI_cross_section(Mcdm, nucleon_mass, proton_SI_anti_amp);
      proton_SD_xsec = pb_to_cm2*
         get_SD_cross_section(Mcdm, nucleon_mass, proton_SD_amp);
      proton_SD_anti_xsec = pb_to_cm2*
         get_SD_cross_section(Mcdm, nucleon_mass, proton_SD_anti_amp);
      neutron_SI_xsec = pb_to_cm2 *
         get_SI_cross_section(Mcdm, nucleon_mass, neutron_SI_amp);
      neutron_SI_anti_xsec = pb_to_cm2 *
         get_SI_cross_section(Mcdm, nucleon_mass, neutron_SI_anti_amp);
      neutron_SD_xsec = pb_to_cm2 *
         get_SD_cross_section(Mcdm, nucleon_mass, neutron_SD_amp);
      neutron_SD_anti_xsec = pb_to_cm2 *
         get_SD_cross_section(Mcdm, nucleon_mass, neutron_SD_anti_amp);

   }

   /* write to file */
   FILE *ofp;
   ofp = fopen(output_filename, "w");

   fprintf(ofp, "Block MICROMEGAS\n");

   if (relic_density_flag) {
      fprintf(ofp, "     1   % 16.8E   # relic density\n", relic_density);
   }

   if (indirect_det_flag) {
      fprintf(ofp, "     2   % 16.8E   # annihilation xsection (cm^3/s)\n", sigma_v);
      fprintf(ofp, "     3   % 16.8E   # test energy (GeV)\n", E_test);
      fprintf(ofp, "     4   % 16.8E   # photon flux ((cm^2 s GeV)^{-1})\n", photon_flux);
      fprintf(ofp, "     5   % 16.8E   # positron flux ((cm^2 s GeV)^{-1})\n", positron_flux);
      fprintf(ofp, "     6   % 16.8E   # antiproton flux ((cm^2 s GeV)^{-1})\n", antiproton_flux);
   }

   if (direct_det_flag) {
      fprintf(ofp, "     7   % 16.8E   # DM-proton SI scattering amplitude\n",
              proton_SI_amp);
      fprintf(ofp, "     8   % 16.8E   # antiDM-proton SI scattering amplitude\n",
              proton_SI_anti_amp);
      fprintf(ofp, "     9   % 16.8E   # DM-proton SD scattering amplitude\n",
              proton_SD_amp);
      fprintf(ofp, "    10   % 16.8E   # antiDM-proton SD scattering amplitude\n",
              proton_SD_anti_amp);
      fprintf(ofp, "    11   % 16.8E   # DM-neutron SI scattering amplitude\n",
              neutron_SI_amp);
      fprintf(ofp, "    12   % 16.8E   # antiDM-neutron SI scattering amplitude\n",
              neutron_SI_anti_amp);
      fprintf(ofp, "    13   % 16.8E   # DM-neutron SD scattering amplitude\n",
              neutron_SD_amp);
      fprintf(ofp, "    14   % 16.8E   # antiDM-neutron SD scattering amplitude\n",
              neutron_SD_anti_amp);
      fprintf(ofp, "    15   % 16.8E   # DM-proton SI xsection (cm^2)\n",
              proton_SI_xsec);
      fprintf(ofp, "    16   % 16.8E   # antiDM-proton SI xsection (cm^2)\n",
              proton_SI_anti_xsec);
      fprintf(ofp, "    17   % 16.8E   # DM-proton SD xsection (cm^2)\n",
              proton_SD_xsec);
      fprintf(ofp, "    18   % 16.8E   # antiDM-proton SD xsection (cm^2)\n",
              proton_SD_anti_xsec);
      fprintf(ofp, "    19   % 16.8E   # DM-neutron SI xsection (cm^2)\n",
              neutron_SI_xsec);
      fprintf(ofp, "    20   % 16.8E   # antiDM-neutron SI xsection (cm^2)\n",
              neutron_SI_anti_xsec);
      fprintf(ofp, "    21   % 16.8E   # DM-neutron SD xsection (cm^2)\n",
              neutron_SD_xsec);
      fprintf(ofp, "    22   % 16.8E   # antiDM-neutron SD xsection (cm^2)\n",
              neutron_SD_anti_xsec);
   }

   fclose(ofp);

   return EXIT_SUCCESS;
}

void parse_options(int* argcp, char*** argvp)
{
   static char* short_opts = "hf:";
   static struct option long_opts[] =
      {
         {"help", no_argument, 0, 'h'},
         {"verbose", no_argument, &verbose_flag, 1},
         {"output-file", required_argument, 0, 'f'},
         {"disable-relic-density", no_argument, &relic_density_flag, 0},
         {"enable-relic-density", no_argument, &relic_density_flag, 1},
         {"disable-indirect-detection", no_argument, &indirect_det_flag, 0},
         {"enable-indirect-detection", no_argument, &indirect_det_flag, 1},
         {"disable-direct-detection", no_argument, &direct_det_flag, 0},
         {"enable-direct-detection", no_argument, &direct_det_flag, 1},
         {0, 0, 0, 0}
      };

   int done;
   while ((done = getopt_long(*argcp, *argvp, short_opts, long_opts, 0)) != -1)
   {
      switch (done) {
      case 0:
         break;
      case 'h':
         print_usage();
         exit(EXIT_SUCCESS);
      case 'f':
         strcpy(output_filename, optarg);
         break;
      case '?':
         break;
      default:
         abort();
      }
   }
}

void print_usage(void)
{
   fprintf(stdout, "Usage: ./run_slha_file.x [OPTION]... FILE\n"
           "Example: ./run_slha_file.x suspect2_lha.out\n\n"
           "Read the given SLHA input file and calculate the\n"
           "resulting relic density in the model.  By default,\n"
           "only the relic density is calculated.\n\n"
           "Options:\n"
           "--enable-relic-density        calculate relic density\n"
           "--disable-relic-density       do not calculate relic density\n"
           "--enable-indirect-detection   "
           "calculate indirect detection observables\n"
           "--disable-indirect-detection  "
           "do not calculate indirect detection observables\n"
           "--enable-direct-detection     "
           "calculate direct detection observables\n"
           "--disable-direct-detection    "
           "do not calculate direct detection observables\n"
           "--output-file                 name of the file to output to\n"
           "--verbose                     print additional information\n"
           "-h, --help                    display this help and exit\n");
}

void set_standard_model_inputs(void)
{
  assignValW("alfSMZ", 0.1184);
  assignValW("aEWinv", 128.8445033);
  assignValW("Gf", 1.16637e-5);
}
