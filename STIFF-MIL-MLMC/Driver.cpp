/*/////////////////////////////////////////////////////////////

  Stiff Problem
    (Single Threaded CPU Version)
  using Milstein method and multilevel Monte Carlo

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

    [1] Split-step Adams–Moulton Milstein methods for systems of
      stiff stochastic differential equations.
    David A. Voss, Department of Mathematics, Western Illinois
      University, Macomb, IL 61455, USA.
    Abdul Q. M. Khaliq, Department of Mathematical Sciences
      and Center for Computational Science, Middle Tennessee
      State University, Murfreesboro, TN 37132, USA.
    International Journal of Computer Mathematics, 
      Vol. 92, No. 5, pp. 995–1011, (2015).
      http://dx.doi.org/10.1080/00207160.2014.915963

    [2] Split-step backward Milstein methods for stiff
      stochastic systems.
    Peng Wang, Institute of Mathematics, Jilin University,
      Changchun 130012, China.
    Yue-cai Han, College of Mathematics, Jilin University,
      Changchun 130012, China.
    Journal of Jilin University (Science Edition)
      Vol. 47, No. 6, pp. 52-56, (2009).

    [3] Stabilized Milstein Type Methods for Stiff
      Stochastic Systems.
    Peng Wang, Institute of Mathematics, Jilin University,
      Changchun 130012, China.
    Zhenxin Liu, College of Mathematics, Jilin University,
      Changchun 130012, China.
    Journal of Numerical Mathematics and Stochastics
      Vol. 1, No. 1, pp. 33-44, (2009).


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
#include "../mlmc_test.h"                                          // master MLMC file
#include "../mrg32k3a.h"                                           // RNG header
#include "Simulation_Parameters.h"                              //

#define WINDOWS
#ifdef WINDOWS
#include <Windows.h>
#endif


void mlmc_l(int, int, double*);

double c1, c2, T;

int main(int argc, char **argv) {

  #ifdef WINDOWS
  HWND consoleWindow = GetConsoleWindow();
  SetWindowPos( consoleWindow, 0, -10, 0, 0, 0, SWP_NOSIZE | SWP_NOZORDER );
  #endif

  int M = __M__;                                                // refinement cost factor
  int initN = __INITIAL_N__;                                    // initial samples on each level
  int Lmin = __LEVEL_MIN__;                                     // minimum refinement level
  int Lmax = __LEVEL_MAX__;                                     // maximum refinement level
  int convergeN = __CONVERGENCE_N__;                            // samples for convergence tests
  int convergeL = __CONVERGENCE_L__;                            // levels for convergence tests 
  double Eps[] = __CONVERGENCE_EPS__;

  FILE *fp;
  std::ofstream outFile;                                        // Write our results to a text file
  outFile.open(__OUT_FILE_NAME2__, std::ios::app);

  int offset = 1234;                                            // RNG initialization
  uint  V1[] = { 1, 2, 3 };                                     // seeds
  uint  V2[] = { 1, 2, 3 };                                     // seeds
  CPU_mrg32k3a_init(V1, V2, offset);

  fopen_s(&fp, __OUT_FILE_NAME1__, "w");
  printf("\n ---- Stiff System ----\n");


  /////////////////////////////////////////////////////////////
  // Page 1008-1009 (15-16/18) [1]
  /////////////////////////////////////////////////////////////
  T = __T__;

  const int nt = __N_TESTS_;                                    // number of tests;
  const int np = 2;                                             // number of parameters per test;

  const double Parameters[nt][np] = {
    __TESTS__
  };

  outFile                                                       // Output the header for a CSV
    << "c1,"
    << "c2,"
    << "T" << std::endl;

  for (int loop = 0; loop < nt; loop++) {                       // Loop over each test
    
    c1 = Parameters[loop][0];
    c2 = Parameters[loop][1];

    outFile                                                     // Output the results to a CSV for analysis
      << c1 << ","
      << c2 << ","
      << T << std::endl;

    mlmc_test(mlmc_l, M, convergeN, convergeL, initN, Eps, Lmin,
      Lmax, fp, &outFile, __FIRST_LEVEL__);

    if(loop<(nt-1))
      outFile
        << std::endl
        << std::endl                                            // Output the header for a CSV
        << "c1,"
        << "c2," 
        << "T" << std::endl;
  }

  fclose(fp);  
  outFile.close();
  return __OKAY__;
}




void mlmc_l(int l, int N, double *sums) {
  
  int   nf, nc;
  double 
    hf, hc, 
    Xf[2], Xc[2],
    dWc[2],
    Pf, Pc, dP,
    dWf[2];
  
  ull   v1[3], v2[3];                                           // needed for RNG
  double x1, x2 = nan("");                                      // needed for Normal RNG
  for (int m = 0; m<3; m++) {                                   // initialize seeds
    v1[m] = CPU_mrg32k3a_v1[m];
    v2[m] = CPU_mrg32k3a_v2[m];
  }
  
  nf = 1 << l + __FIRST_LEVEL__;
  nc = nf / 2;
  hf = T / ((double)nf);
  hc = T / ((double)nc);

  for (int k = 0; k < 6; k++) sums[k] = 0.0;

  for (int np = 0; np<N; np++) {

    Xc[0] = __Y10__;                                            //
    Xc[1] = __Y20__;                                            //
    Xf[0] = __Y10__;                                            //
    Xf[1] = __Y20__;                                            //

    
    
    if (l == 0) {
    //if (l == __FIRST_LEVEL__) {
      for (int n = 0; n < nf; n++) {

        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[0] = sqrt(hf)*x1;
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[1] = sqrt(hf)*x1;

        Xf[0] *= 1 + c1*hf;                                        // Euler and Maruyama terms
        Xf[1] += c2*Xf[1]*hf + 4.0*Xf[0]*dWf[0] + 0.5*dWf[1];      // Euler and Maruyama terms

      }

    }
    else {


      for (int n = 0; n<nc; n++) { 

        dWc[0] = 0.0;
        dWc[1] = 0.0;
          
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[0] = sqrt(hf)*x1;
        dWc[0] += dWf[0];
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[1] = sqrt(hf)*x1;
        dWc[1] += dWf[1];
        
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[2] = sqrt(hf)*x1;
        dWc[0] += dWf[2];
        CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
        dWf[3] = sqrt(hf)*x1;
        dWc[1] += dWf[3];

        Xf[0] *= 1 + c1*hf;                                        // Euler and Maruyama terms
        Xf[1] += c2*Xf[1]*hf + 4.0*Xf[0]*dWf[0] + 0.5*dWf[1];      // Euler and Maruyama terms

        Xf[0] *= 1 + c1*hf;                                        // Euler and Maruyama terms
        Xf[1] += c2*Xf[1]*hf + 4.0*Xf[0]*dWf[2] + 0.5*dWf[3];      // Euler and Maruyama terms

        Xc[0] *= 1 + c1*hc;                                        // Euler and Maruyama terms
        Xc[1] += c2*Xc[1]*hc + 4.0*Xc[0]*dWc[0] + 0.5*dWc[1];      // Euler and Maruyama terms

      }


    }
    
    
    Pf = Xf[1];
    Pc = Xc[1];

    if (l == 0) dP = Pf;
    else dP = Pf - Pc;

    sums[0] += dP;
    sums[1] += dP*dP;
    sums[2] += dP*dP*dP;
    sums[3] += dP*dP*dP*dP;
    sums[4] += Pf;
    sums[5] += Pf*Pf;
  }
}
