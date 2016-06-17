/*/////////////////////////////////////////////////////////////

  Gene Transcription Problem
    (Single Threaded CPU Version)
  using Milstein method and antithetic multilevel Monte Carlo


  Zane Colgin
  Middle Tennessee State University
  March 2016

  DATE        AUTHOR  COMMENTS
  ----------  ------  ---------------
  ????-??-??  MG      Original version
  2016-01-30  ZC      Modified version exchange option problem
  2016-03-07  ZC      Modified for gene transcription problem


  NOTE:
    "!!!" in comments indicates programming note 
      i.e. look here for debugging notes, optimizations, etc.
    "???" in comments indicates missing information or 
      unfinished

  !!! IMPORTANT:
    see Simulation_Parameters.h to adjust critical definitions
    including file system specifics for output file location


  REFERENCES:
    [1] D. F. Anderson and D. J. Higham. Multilevel Monte Carlo 
        for continuous time Markov chains, with applications in 
        biochemical kinetics, Multiscale Model. Simul., 10
        (2012), pp. 146-179.


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
#include "../mlmc_test.h"                                          // MLMC test file
#include "../mrg32k3a.h"                                           // RNG header
#include "Simulation_Parameters.h"                              // !!! IMPORTANT: file system specifics! see output file names
#include "../geneModel.h"

#define WINDOWS
#ifdef WINDOWS
#include <Windows.h>
#endif

///////////////////////////////////////////////////////////////
// Prototypes
///////////////////////////////////////////////////////////////

void mlmc_l(int, int, double *);
void solve_timestep(double* Y1p, double* Y2p, double* Y3p, double h,
  double dW1, double dW2, double dW3, double dW4, double dW5);

///////////////////////////////////////////////////////////////
// Global Variables (model parameters)
// !!! Bad practice.  Consider creating a parameters sturcture
// !!! and modifiy the appropriate functions to pass it along
///////////////////////////////////////////////////////////////

static const double T = __T__;
static const double S0[geneModel::kNumSpecies] = __S0__;
static const double G0 = __G0__;


///////////////////////////////////////////////////////////////
// MAIN
///////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  #ifdef WINDOWS
  HWND consoleWindow = GetConsoleWindow();
  SetWindowPos( consoleWindow, 0, -10, 0, 0, 0, SWP_NOSIZE | SWP_NOZORDER );
  #endif

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
  printf("\n ---- Gene Transcription Model ----\n");  
  
  mlmc_test(mlmc_l, M, convergeN, convergeL, initN, Eps, Lmin,//////// PRIMARY FUNCTION CALL
    Lmax, fp, &outFile, __FIRST_LEVEL__);
  
  fclose(fp);
  outFile.close();
  return __OKAY__;
}




void mlmc_l(int l, int N, double *sums) {

  int   nf, nc;
  double 
    hf, hc, 
    Xf[geneModel::kNumSpecies], 
    Xa[geneModel::kNumSpecies],
    Xc[geneModel::kNumSpecies],
    dWc[geneModel::kNumReactions],
    Pf, Pc, dP,
    dWf1[geneModel::kNumReactions],                             // 1st Weiner increment for fine level (all dim)
    dWf2[geneModel::kNumReactions];                             // 2st Weiner increment for fine level (all dim)
  
  ull   v1[3], v2[3];                                           // needed for RNG
  double x1, x2 = nan("");                                      // needed for Normal RNG
  for (int m = 0; m<3; m++) {                                   // initialize seeds
    v1[m] = CPU_mrg32k3a_v1[m];
    v2[m] = CPU_mrg32k3a_v2[m];
  }
  
  nf = 1 << l + __FIRST_LEVEL__;                                // number of steps for fine discretization
  nc = nf / 2;                                                  // number of steps for coarse discretization
  hf = T / ((double)nf);                                        // time step for fine discretization
  hc = T / ((double)nc);                                        // time step for coarse discretization

  for (int k = 0; k < 6; k++) sums[k] = 0.0;

  for (int np = 0; np<N; np++) {

    Xc[0] = S0[0];                                              // initial values, Eq1, coarse discretization
    Xc[1] = S0[1];                                              // initial values, Eq2, coarse discretization
    Xc[2] = S0[2];                                              // initial values, Eq3, coarse discretization
    Xf[0] = S0[0];                                              // initial values, Eq1, fine discretization
    Xf[1] = S0[1];                                              // initial values, Eq2, fine discretization
    Xf[2] = S0[2];                                              // initial values, Eq3, fine discretization


    if (l == 0) {                                               // For initial level
      for (int n = 0; n < nf; n++) {

        for (int r = 0; r < geneModel::kNumReactions; r++) {
          CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
          dWf1[r] = sqrt(hf)*x1;
        }

        solve_timestep( &(Xf[0]), &(Xf[1]), &(Xf[2]), hf, 
          dWf1[0], dWf1[1], dWf1[2], dWf1[3], dWf1[4]);

      }
    }
    else {                                                      // For correction levels
      for (int n = 0; n<nc; n++) { 

        for (int r = 0; r < geneModel::kNumReactions; r++) {    // Weiner increment, fine discretization refinement step 1/2
          dWc[r] = 0.0;
          CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
          dWf1[r] = sqrt(hf)*x1;
          dWc[r] += dWf1[r];
        }
        for (int r = 0; r < geneModel::kNumReactions; r++) {    // Weiner increment, fine discretization refinement step 2/2
          CPU_mrg32k3a_next_normal(v1, v2, x1, x2);
          dWf2[r] = sqrt(hf)*x1;
          dWc[r] += dWf2[r];
        }

        Xa[0] = Xf[0];
        Xa[1] = Xf[1];
        Xa[2] = Xf[2];
        
        solve_timestep( &(Xf[0]), &(Xf[1]), &(Xf[2]), hf,       // fine increment step 1/2
          dWf1[0], dWf1[1], dWf1[2], dWf1[3], dWf1[4]);

        solve_timestep( &(Xf[0]), &(Xf[1]), &(Xf[2]), hf,       // fine increment step 2/2
          dWf2[0], dWf2[1], dWf2[2], dWf2[3], dWf2[4]);

        solve_timestep( &(Xa[0]), &(Xa[1]), &(Xa[2]), hf,       // antithetic increment step 1/2
          dWf2[0], dWf2[1], dWf2[2], dWf2[3], dWf2[4]);

        solve_timestep( &(Xa[0]), &(Xa[1]), &(Xa[2]), hf,       // antithetic increment step 2/2
          dWf1[0], dWf1[1], dWf1[2], dWf1[3], dWf1[4]);

        Xf[0] = (Xf[0] + Xa[0])/2.0;
        Xf[1] = (Xf[1] + Xa[1])/2.0;
        Xf[2] = (Xf[2] + Xa[2])/2.0;

        solve_timestep( &(Xc[0]), &(Xc[1]), &(Xc[2]), hc,       // coarse increment
          dWc[0], dWc[1], dWc[2], dWc[3], dWc[4]);

      }
    }

    Pf = Xf[2];
    Pc = Xc[2];


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








void solve_timestep(double* Y1p, double* Y2p, double* Y3p, double h,
  double dW1, double dW2, double dW3, double dW4, double dW5){

  static const double c1 = geneModel::kC[0];                           // c1
  static const double c2 = geneModel::kC[1];                           // c2
  static const double c3 = geneModel::kC[2];                           // c3
  static const double c4 = geneModel::kC[3];                           // c4
  static const double c5 = geneModel::kC[4];                           // c5
  
  static const double sqc1 = sqrt(c1);                                 // sqrt(c1)
  static const double sqc2 = sqrt(c2);                                 // sqrt(c2)
  static const double sqc3 = sqrt(c3);                                 // sqrt(c3)
  static const double sqc4 = sqrt(c4);                                 // sqrt(c4)
  static const double sqc5 = sqrt(c5);                                 // sqrt(c5)

  static const double c1G0 = c1*G0;                                    // c1*G0
  static const double sqc1G0 =  sqrt(c1*G0);                           // c1*G0

  static const double sq2c3 = sqrt(2*c3);                              // sqrt(2*c3)
  static const double c3_2 = c3/2.0;                                   // c3/2
  static const double sqc3_2 = sqrt(c3_2);                             // sqrt(c3/2)
  static const double sq2 = sqrt(2.0);

  double Y1, Y2, Y3;
  double sqY1, sqY2, sqY3;
  double Y2Y21, sqY2Y21, sqY21, Y205;
  double dW11t,  dW12, dW13, dW14, dW15,
                dW22t, dW23, dW24, dW25,
                      dW33t, dW34, dW35,
                            dW44t, dW45,
                                  dW55t;

  double dMil1, dMil2, dMil3;

  Y1 = *Y1p;
  Y2 = *Y2p;
  Y3 = *Y3p;
  sqY1 = sqrt(Y1);
  sqY2 = sqrt(Y2);
  sqY3 = sqrt(Y3);

  Y2Y21 = Y2*(Y2 - 1.0);
  Y205 = Y2 - 0.5;
  //if ((Y2-1.0)>=0)
  //  sqY21 = sqrt(Y2 - 1.0);
  //else
  //  sqY21 = 1.0;
  //if (Y2Y21>=0)
  //  sqY2Y21 = sqrt(Y2Y21);
  //else
  //  sqY2Y21 = 1.0;


  sqY21 = sqrt(Y2 - 1.0);
  sqY2Y21 = sqrt(Y2Y21);


  //if (sqY21!=sqY21)
  //  printf("sqY21 is NaN\n"); 
  //if (sqY2Y21!=sqY2Y21)
  //  printf("sqY2Y21 is NaN\n");

  dW11t = dW1*dW1-h;
  dW12  = dW1*dW2;
  dW13  = dW1*dW3;
  dW14  = dW1*dW4;
  dW15  = dW1*dW5;
  dW22t = dW2*dW2-h;
  dW23  = dW2*dW3;
  dW24  = dW2*dW4;
  dW25  = dW2*dW5;
  dW33t = dW3*dW3-h;
  dW34  = dW3*dW4;
  dW35  = dW3*dW5;
  dW44t = dW4*dW4-h;
  dW45  = dW4*dW5;
  dW55t = dW5*dW5-h;

        
  Y1 += (c1G0 - c4*Y1)*h                                  // Euler term
    + sqc1G0*dW1 - sqc4*sqY1*dW4;                         // Maruyama term
  dMil1 = 0.5*(                                           // Milstein term
       c4 / 2.0 * dW44t            
      - sqc4*sqc1G0 / 2.0 / sqY1*dW14
   );
  Y2 += (c2*Y1 - c3*Y2Y21 - c5*Y2)*h                      // Euler term
    + sqc2*sqY1*dW2 - sq2c3*sqY2Y21*dW3 - sqc5*sqY2*dW5;  // Maruyama term
  dMil2 = 0.5*(                                           // Milstein term
        sqc1G0*sqc2 / 2.0 / sqY1*dW12                     //    1
      - sqc2*sqc4 / 2.0 * dW24                            //    2
      - sq2c3*sqc2*sqY1 / sqY2Y21*Y205*dW23               //    3
      + 2.0 * c3*Y205*dW33t                               //    4
      + sq2c3*sqc5 / sqY21*Y205*dW35                      //    5
      - sqc2*sqc5*sqY1 / 2.0 / sqY2*dW25                  //    6
      + sqc3*sqc5*sqY21 / sq2*dW35                        //    7
      + c5 / 2.0*dW55t                                    //    8
  );
  Y3 += c3_2*Y2Y21*h                                      // Euler term
    + sqc3_2*sqY2Y21*dW3;                                 // Maruyama term
  dMil3 = 0.5*(                                           // Milstein term
        sqc2*sqc3*sqY1 / sq2 / sqY2Y21*Y205*dW23
      - c3*Y205*dW33t
      - sqc3*sqc5 / sq2 / sqY21*Y205*dW35
    );

  //if (~std::isfinite(dMil1))
  //  dMil1 = 0;
  //if (~std::isfinite(dMil2))
  //  dMil2 = 0;
  //if (~std::isfinite(dMil3))
  //  dMil3 = 0;
        
  *Y1p = fmax(0.0, Y1+dMil1);
  *Y2p = fmax(0.0, Y2+dMil2);
  *Y3p = fmax(0.0, Y3+dMil3);

  return;

}


