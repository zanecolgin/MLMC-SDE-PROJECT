/*/////////////////////////////////////////////////////////////

  2-Assest exchange option
    (Single Threaded CPU Version)
  using antithetic multilevel Monte Carlo


  Zane Colgin
  Middle Tennessee State University
  March 2016

  DATE        AUTHOR  COMMENTS
  ----------  ------  ---------------
  ????-??-??  MG      Original version
  2016-01-30  ZC      Modified version
  2016-03-20  ZC      Modified version


  NOTE:
    "!!!" in comments indicates programming note 
      i.e. look here for debugging notes, optimizations, etc.
    "???" in comments indicates missing information or 
      unfinished

  !!! IMPORTANT:
    see Simulation_Parameters.h to adjust critical definitions
    including file system specifics for output file location


  REFERENCES:
    [1] Analytic Approximations for Multi-Asset Option Pricing
      Carol Alexander, ICMA Centre, University of Reading
      Aanand Venkatramanan, ICMA Centre, University of Reading
      Version: December 2009
      ICMA Centre Discussion Papers in Finance DP2009-05
      Copyright 2009 Alexander and Venkatramanan.
    [2] Antithetic Multilevel Monte Carlo Estimation for Multi-
      dimensional SDEs Without L´evy Area Simulation
      Michael B. Giles, University of Oxford
      Lukasz Szpruch, University of Oxford
      The Annals of Applied Probability
      2014, Vol. 24, No. 4, 1585–1620
      DOI: 10.1214/13-AAP957
      Copyright 2014 Institute of Mathematical Statistics

  THANKS:
    Mike Giles -
      This code can be found in its original form at
      http://people.maths.ox.ac.uk/gilesm/mlmc
      Slight modifications have been made throughout the
      associated files (mlmc.cpp, mlmc_test.cpp, mrg32k3a.h)
      and major modifications have beem made to this file
      especially in the function mlmc_l(...)
      (originally named mcqmc06_l)


  ORIGINAL COMMENT (Mike Giles):
    % These are similar to the MLMC tests for the MCQMC06 paper
    % using a Milstein discretisation with 2^l timesteps on 
    % level l
    %
    % The figures are slightly different due to
    % -- change in MSE split
    % -- change in cost calculation
    % -- different random number generation
    % -- switch to S_0=100

*//////////////////////////////////////////////////////////////

#pragma once
#include "../mlmc_test.h"                                       // MLMC test file
#include "../mrg32k3a.h"                                        // RNG header
#include "Simulation_Parameters.h"                              // !!! IMPORTANT: file system specifics! see output file names


///////////////////////////////////////////////////////////////
// Prototypes
///////////////////////////////////////////////////////////////

void mlmc_l(int, int, double *);


///////////////////////////////////////////////////////////////
// Global Variables (model parameters)
// !!! Bad practice.  Consider creating a parameters sturcture
// !!! and modifiy the appropriate functions to pass it along
///////////////////////////////////////////////////////////////

double S1, S2, r, T, sigma1, sigma2, rho;


///////////////////////////////////////////////////////////////
// MAIN
///////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  const int M = __M__;                                          // refinement cost factor !!! must be 2
  const int initN = __INITIAL_N__;                              // initial samples on each level
  const int Lmin = __LEVEL_MIN__;                               // minimum number of refinement levels
  const int Lmax = __LEVEL_MAX__;                               // maximum number of refinement levels
  const int convergeN = __CONVERGENCE_N__;                      // samples for convergence tests
  const int convergeL = __CONVERGENCE_L__;                      // levels for convergence tests 
  double Eps[] = __CONVERGENCE_EPS__;
  const int offset = 1234;                                      // for RNG initialization

  FILE *fp;                                                     // Write our results to a text file
  std::ofstream outFile;                                        // Write our results to a csv file

  uint  V1[] = { 1, 2, 3 };                                     // seeding
  uint  V2[] = { 1, 2, 3 };                                     // seeding
  CPU_mrg32k3a_init(V1, V2, offset);                            // initialize RNG

  fopen_s(&fp, __OUT_FILE_NAME1__, "w");
  outFile.open(__OUT_FILE_NAME2__, std::ios::app);
  printf("\n ---- option: Exchange ----\n");


  /////////////////////////////////////////////////////////////
  // TABLE 2 page 28 [1]
  /////////////////////////////////////////////////////////////
  S1 = __S1__;

  const int nt = 39;                                            // number of tests;
  const int np = 9;                                             // number of parameters per test;

  // TABLE 2: Exchange option prices ([1] page 28)
  const double Parameters[nt][np] = {
    // S2     r  T sig1 sig2  rho      GBM  StdDev   Approx
    {  80, 0.00, 1, 0.1, 0.1, 0.5, 20.0000, 0.0001, 19.9999 },
    {  80, 0.00, 1, 0.2, 0.2, 0.5, 20.0000, 0.0003, 19.9999 },
    {  80, 0.00, 1, 0.1, 0.1, 0.8, 20.0000, 0.0001, 20.0009 },
    {  80, 0.00, 1, 0.2, 0.1, 0.8, 20.0003, 0.0003, 20.0023 },
    {  80, 0.00, 6, 0.1, 0.1, 0.5, 20.0012, 0.0004, 20.0013 },
    {  80, 0.00, 6, 0.2, 0.2, 0.5, 20.3073, 0.0020, 20.3123 },
    {  80, 0.00, 6, 0.1, 0.1, 0.8, 20.0003, 0.0003, 20.0052 },
    {  80, 0.04, 1, 0.1, 0.1, 0.5, 20.0667, 0.0001, 20.0667 },
    {  80, 0.04, 1, 0.2, 0.1, 0.5, 20.0667, 0.0003, 20.0664 },
    {  80, 0.04, 6, 0.1, 0.1, 0.8, 20.4043, 0.0003, 20.4093 },
    {  80, 0.04, 6, 0.2, 0.1, 0.8, 20.4327, 0.0015, 20.4444 },
    {  90, 0.00, 1, 0.1, 0.1, 0.5, 10.0001, 0.0001, 10.0000 },
    {  90, 0.00, 1, 0.2, 0.2, 0.5, 10.0727, 0.0005, 10.0744 },
    {  90, 0.00, 1, 0.2, 0.1, 0.8, 10.0039, 0.0003, 10.0084 },
    {  90, 0.00, 6, 0.2, 0.1, 0.8, 10.6073, 0.0024, 10.6123 },
    {  90, 0.00, 6, 0.1, 0.1, 0.8, 10.0133, 0.0004, 10.0225 },
    {  90, 0.04, 1, 0.1, 0.1, 0.5, 10.0335, 0.0001, 10.0334 },
    {  90, 0.04, 1, 0.2, 0.1, 0.5, 10.0633, 0.0004, 10.0635 },
    {  90, 0.04, 6, 0.2, 0.1, 0.5, 11.4786, 0.0034, 11.4833 },
    {  90, 0.04, 6, 0.2, 0.2, 0.5, 12.0095, 0.0033, 12.0155 },
    {  90, 0.04, 6, 0.2, 0.1, 0.8, 10.8216, 0.0025, 10.8267 },
    { 100, 0.00, 1, 0.2, 0.1, 0.5,  1.9934, 0.0016,  1.9924 },
    { 100, 0.00, 1, 0.1, 0.1, 0.5,  1.1516, 0.0009,  1.1487 },
    { 100, 0.00, 6, 0.2, 0.1, 0.5,  4.8802, 0.0043,  4.8767 },
    { 100, 0.00, 6, 0.1, 0.1, 0.5,  2.8203, 0.0022,  2.8131 },
    { 100, 0.04, 1, 0.2, 0.2, 0.5,  2.3106, 0.0018,  2.3048 },
    { 100, 0.04, 6, 0.2, 0.1, 0.5,  4.9788, 0.0044,  4.9753 },
    { 100, 0.04, 6, 0.1, 0.1, 0.5,  2.8773, 0.0022,  2.8699 },
    { 110, 0.00, 1, 0.1, 0.1, 0.5,  0.0004, 0.0001,  0.0000 },
    { 110, 0.00, 1, 0.2, 0.1, 0.5,  0.0569, 0.0004,  0.0413 },
    { 110, 0.00, 1, 0.2, 0.2, 0.5,  0.1242, 0.0006,  0.0903 },
    { 110, 0.00, 1, 0.2, 0.2, 0.8,  0.0053, 0.0001,  0.0000 },
    { 110, 0.00, 1, 0.2, 0.1, 0.8,  0.0091, 0.0002,  0.0003 },
    { 110, 0.04, 1, 0.1, 0.1, 0.5,  0.0004, 0.0001,  0.0000 },
    { 110, 0.04, 1, 0.2, 0.1, 0.5,  0.0571, 0.0004,  0.0414 },
    { 110, 0.04, 1, 0.2, 0.2, 0.5,  0.1246, 0.0006,  0.0906 },
    { 110, 0.04, 1, 0.2, 0.2, 0.8,  0.0054, 0.0001,  0.0000 },
    { 110, 0.04, 6, 0.2, 0.1, 0.5,  1.6339, 0.0031,  1.5934 },
    { 110, 0.04, 6, 0.1, 0.1, 0.5,  0.3103, 0.0010,  0.2612 }
  };

  outFile                                                       // Output the header for a CSV
    << "S2,"
    << "r,"
    << "T,"
    << "sigma1,"
    << "sigma2,"
    << "rho,"
    << "GBM" << std::endl;


  /////////////////////////////////////////////////////////////
  // Main loop
  /////////////////////////////////////////////////////////////

  for (int loop = 0; loop < nt; loop++) {                       // Loop over each test

    S2 = Parameters[loop][0];
    r = Parameters[loop][1];
    T = Parameters[loop][2] / 12.0;                             // convert months to years
    sigma1 = Parameters[loop][3];
    sigma2 = Parameters[loop][4];
    rho = Parameters[loop][5];
    double GBM = Parameters[loop][6];

    outFile                                                     // Output the results to a CSV for analysis
      << S2 << ","
      << r << ","
      << T << ","
      << sigma1 << ","
      << sigma2 << ","
      << rho << ","
      << GBM << std::endl;

    mlmc_test(mlmc_l, M, convergeN, convergeL, initN, Eps, Lmin,//////// PRIMARY FUNCTION CALL
      Lmax, fp, &outFile, __FIRST_LEVEL__);

    if(loop<(nt-1))                                             // Output the header for a CSV
      outFile
        << std::endl
        << std::endl
        << "S2,"
        << "r,"
        << "T,"
        << "sigma1,"
        << "sigma2,"
        << "rho,"
        << "GBM" << std::endl;
  }

  fclose(fp);
  outFile.close();
  return __OKAY__;
}




void mlmc_l(int l, int N, double *sums) {

  const double crho = sqrt(1.0 - rho*rho);

  int   nf, nc;
  double 
    hf, hc, 
    Xa[2], Xf[2], Xc[2],
    dWc[2],
    Pf, Pc, dP,
    dWf[2*2];
  
  ull   v1[3], v2[3];                                           // needed for RNG
  double x1, x2 = nan("");                                      // needed for Normal RNG
  for (int m = 0; m<3; m++) {                                   // initialize seeds
    v1[m] = CPU_mrg32k3a_v1[m];
    v2[m] = CPU_mrg32k3a_v2[m];
  }
  
  nf = 1 << l;                                                  // number of steps for fine discretization
  nc = nf / 2;                                                  // number of steps for coarse discretization
  hf = T / ((double)nf);                                        // time step for fine discretization
  hc = T / ((double)nc);                                        // time step for coarse discretization

  for (int k = 0; k < 6; k++) sums[k] = 0.0;

  for (int np = 0; np<N; np++) {

    Xc[0] = __S1__;                                             // initial values, Eq1, coarse discretization
    Xc[1] = S2;                                                 // initial values, Eq2, coarse discretization
    Xf[0] = __S1__;                                             // initial values, Eq1, fine discretization
    Xf[1] = S2;                                                 // initial values, Eq2, fine discretization


    if (l == 0) {                                               // For initial level
      for (int n = 0; n < nf; n++) {

        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[0] = sqrt(hf)*x1;
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[1] = sqrt(hf)*x1;

        Xf[0] *= 1 + r*hf +     sigma1*dWf[0]                   // Euler and Maruyama terms
                  + 0.5*sigma1*sigma1*(dWf[0] * dWf[0] - hf);   // Milstein term
        Xf[1] *= 1 + r*hf                                       // Euler terms
                  + rho*sigma2*dWf[0] + crho*sigma2*dWf[1]      // Maruyama term
                  + 0.5*sigma2*sigma2*(                         // Milstein term
                    rho*rho*(dWf[0] * dWf[0] - hf) +
                    2 * rho*crho*dWf[0] * dWf[1] +
                    crho*crho*(dWf[1] * dWf[1] - hf)  
                    );

      }
    }
    else {                                                      // For correction levels
      for (int n = 0; n<nc; n++) { 

        dWc[0] = 0.0;
        dWc[1] = 0.0;
        
        Xa[0] = Xf[0];
        Xa[1] = Xf[1];

        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[0] = sqrt(hf)*x1;                                   // Weiner increment, Dim 1, fine discretization refinement step 1/2
        dWc[0] += dWf[0];
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[1] = sqrt(hf)*x1;                                   // Weiner increment, Dim 2, fine discretization refinement step 1/2
        dWc[1] += dWf[1];
        
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[2] = sqrt(hf)*x1;                                   // Weiner increment, Dim 1, fine discretization refinement step 2/2
        dWc[0] += dWf[2];                                       // Weiner increment, Dim 1, coarse discretization
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[3] = sqrt(hf)*x1;                                   // Weiner increment, Dim 2, fine discretization refinement step 2/2
        dWc[1] += dWf[3];                                       // Weiner increment, Dim 2, coarse discretization


        Xf[0] *= 1 + r*hf +     sigma1*dWf[0]                     // Euler and Maruyama terms, Eq 1, fine discretization, refinement step 1/2
                    + 0.5*sigma1*sigma1*(dWf[0] * dWf[0] - hf);   // Milstein term,            Eq 1, fine discretization, refinement step 1/2
        Xf[1] *= 1 + r*hf + rho*sigma2*dWf[0] + crho*sigma2*dWf[1]// Euler and Maruyama terms, Eq 2, fine discretization, refinement step 1/2
                    + 0.5*sigma2*sigma2*(                         // Milstein term,            Eq 2, fine discretization, refinement step 1/2
                      rho*rho*(dWf[0] * dWf[0] - hf) +
                      2 * rho*crho*dWf[0] * dWf[1] +
                      crho*crho*(dWf[1] * dWf[1] - hf)  
                      );
        Xf[0] *= 1 + r*hf +     sigma1*dWf[2]                     // Euler and Maruyama terms, Eq 1, fine discretization, refinement step 2/2
                    + 0.5*sigma1*sigma1*(dWf[2] * dWf[2] - hf);   // Milstein ter,             Eq 1, fine discretization, refinement step 2/2
        Xf[1] *= 1 + r*hf + rho*sigma2*dWf[2] + crho*sigma2*dWf[3]// Euler and Maruyama terms, Eq 2, fine discretization, refinement step 2/2
                    + 0.5*sigma2*sigma2*(                         // Milstein term,            Eq 2, fine discretization, refinement step 2/2
                      rho*rho*(dWf[2] * dWf[2] - hf) +
                      2 * rho*crho*dWf[2] * dWf[3] +
                      crho*crho*(dWf[3] * dWf[3] - hf)  
                      );

        
        Xa[0] *= 1 + r*hf +     sigma1*dWf[2]                     // Euler and Maruyama terms
                    + 0.5*sigma1*sigma1*(dWf[2] * dWf[2] - hf);   // Milstein term
        Xa[1] *= 1 + r*hf + rho*sigma2*dWf[2] + crho*sigma2*dWf[3]// Euler and Maruyama terms
                    + 0.5*sigma2*sigma2*(                         // Milstein term
                      rho*rho*(dWf[2] * dWf[2] - hf) +
                      2 * rho*crho*dWf[2] * dWf[3] +
                      crho*crho*(dWf[3] * dWf[3] - hf)  
                      );
        Xa[0] *= 1 + r*hf +     sigma1*dWf[0]                     // Euler and Maruyama terms
                    + 0.5*sigma1*sigma1*(dWf[0] * dWf[0] - hf);   // Milstein term
        Xa[1] *= 1 + r*hf + rho*sigma2*dWf[0] + crho*sigma2*dWf[1]// Euler and Maruyama terms
                    + 0.5*sigma2*sigma2*(                         // Milstein term
                      rho*rho*(dWf[0] * dWf[0] - hf) +
                      2 * rho*crho*dWf[0] * dWf[1] +
                      crho*crho*(dWf[1] * dWf[1] - hf)  
                      );


        Xf[0] = (Xf[0] + Xa[0])/2.0;
        Xf[1] = (Xf[1] + Xa[1])/2.0;


        Xc[0] *= 1 + r*hc +     sigma1*dWc[0]                      // Euler and Maruyama terms, Eq 1, coarse discretization
                    + 0.5*sigma1*sigma1*(dWc[0] * dWc[0] - hc);    // Milstein term,            Eq 1, coarse discretization
        Xc[1] *= 1 + r*hc + rho*sigma2*dWc[0] + crho*sigma2*dWc[1] // Euler and Maruyama terms, Eq 2, coarse discretization
                    + 0.5*sigma2*sigma2*(                          // Milstein term,            Eq 2, coarse discretization
                      rho*rho*(dWc[0] * dWc[0] - hc) +
                      2 * rho*crho*dWc[0] * dWc[1] +
                      crho*crho*(dWc[1] * dWc[1] - hc)  
                      );

      }
    }

    Pf = fmax(0.0, Xf[0] - Xf[1]);                              // Payoff function for exchange option, fine discretization
    Pc = fmax(0.0, Xc[0] - Xc[1]);                              // Payoff function for exchange option, coarse discretization

    if (l == 0) dP = Pf;                                        // if initial level, return the approximation as the QoI
    else dP = Pf - Pc;                                          // if correction level, correction is QoI

    sums[0] += dP;                                              // QoI
    sums[1] += dP*dP;                                           // for second moment/variance 
    sums[2] += dP*dP*dP;
    sums[3] += dP*dP*dP*dP;
    sums[4] += Pf;                                              // fine approximation ONLY at this level
    sums[5] += Pf*Pf;
  }
}