/**
 * @file CSE6SSM_utilities.h
 * @brief contains prototypes for model-specific functions
 */

#ifndef CSE6SSM_H
#define CSE6SSM_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief prints the mass and quantum numbers of the given particle
 *
 * @param candidate the name of the particle to print information for
 */
void print_DM_candidate_info(char *candidate);

/**
 * @brief returns 1 if two DM candidates are present in the model
 *
 * @return 1 if two DM candidates have been found, 0 otherwise
 */
int has_two_DM_candidates(void);

/**
 * @brief calculates the relic density for one or two DM candidates
 *
 * @param fast use a fast calculation (fast = 1) or high precision (fast = 0)
 * @param thresh minimum value of Boltzmann suppression factor to allow
 * @param freeze_out use freeze-out approximation
 * @return the total relic density
 */
double calculate_relic_density(int fast, double thresh, int freeze_out);

/**
 * @brief calculates the relic density for a single DM candidate
 *
 * @param fast use a fast calculation (fast = 1) or high precision (fast = 0)
 * @param thresh minimum value of Boltzmann suppression factor to allow
 * @param xf the ratio of the DM mass to the freeze-out temperature
 * @param freeze_out use freeze-out approximation
 * @return the total relic density
 */
double calculate_single_DM_relic_density(int fast, double thresh, double *xf,
                                         int freeze_out);

/**
 * @brief calculates spectra for DM annihilation at rest
 *
 * @param key integer value indicating which effects to include
 * @param gamma annihilation spectrum for photons
 * @param eplus annihilation spectrum for positrons
 * @param pbar annihilation spectrum for antiprotons
 * @param nu_e annihilation spectrum for electron neutrinos
 * @param nu_mu annihilation spectrum for muon neutrinos
 * @param nu_tau annihilation spectrum for tau neutrinos
 * @param error error code, non-zero if a problem is encountered
 * @return cross section times velocity in cm^3/s
 */
double calculate_annihilation_spectra(int key, double *gamma, double *eplus,
                                      double *pbar, double *nu_e,
                                      double *nu_mu, double *nu_tau,
                                      int *error);

/**
 * @brief calculate photon flux from photon spectrum
 *
 * @param fi angle of sight in radians
 * @param dfi half of cone angle in radians
 * @param sigma_v annihilation cross section in cm^3/s
 * @param spectrum_gamma calculated photon spectrum
 * @param flux_gamma table containing the resulting photon flux
 */
void get_photon_flux_table(double fi, double dfi, double sigma_v,
                           double *spectrum_gamma, double *flux_gamma);

/**
 * @brief calculate positron flux from positron spectrum
 *
 * @param E_min cut-off energy in GeV
 * @param sigma_v annihilation cross section in cm^3/s
 * @param spectrum_eplus calculated positron spectrum
 * @param flux_eplus table containing the resulting positron flux
 */
void get_positron_flux_table(double E_min, double sigma_v,
                             double *spectrum_eplus, double *flux_eplus);

/**
 * @brief calculate antiproton flux from antiproton spectrum
 *
 * @param E_min cut-off energy in GeV
 * @param sigma_v annihilation cross section in cm^3/s
 * @param spectrum_pbar calculated antiproton spectrum
 * @param flux_pbar table containing the resulting antiproton flux
 */
void get_antiproton_flux_table(double E_min, double sigma_v,
                               double *spectrum_pbar, double *flux_pbar);

/**
 * @brief calculate the amplitudes for DM-nucleon elastic scattering
 *
 * @param candidate DM candidate to calculate for
 * @param loop_order function for setting accuracy of amplitude calculation
 * @param pA0 spin-independent amplitude for scattering on protons
 * @param pA5 spin-dependent amplitude for scattering on protons
 * @param nA0 spin-independent amplitude for scattering on neutrons
 * @param nA5 spin-dependent amplitude for scattering on neutrons
 */
int calculate_nucleon_amplitudes(char *candidate,
                                 double (*loop_order)(double,double,
                                                      double,double),
                                 double *pA0, double *pA5, double *nA0,
                                 double *nA5);

/**
 * @brief calculate the DM-proton spin-independent cross-section
 *
 * @param mdm mass of the DM candidate
 * @param mp mass of the target
 * @param amp spin-independent amplitude for scattering on target
 * @return cross section in pb
 */
double get_SI_cross_section(double mdm, double mp, double amp);

/**
 * @brief calculate the spin-dependent cross-section
 *
 * @param mdm mass of the DM candidate
 * @param mp mass of the target
 * @param amp spin-dependent amplitude for scattering on target
 * @return cross section in pb
 */
double get_SD_cross_section(double mdm, double mp, double amp);

#ifdef __cplusplus
}
#endif

#endif
