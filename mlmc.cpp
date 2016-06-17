/*/////////////////////////////////////////////////////////////

  Multilevel Monte Carlo
    (Single Threaded CPU Version)

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


  REFERENCES:
    ...see Mike Giles MCQMC06...


  THANKS:
    Mike Giles -
    This code can be found in its original form at
    http://people.maths.ox.ac.uk/gilesm/mlmc
    Slight modifications have been made throughout the
    associated files (mlmc.cpp, mlmc_test.cpp, mrg32k3a.h)
    and major modifications have beem made to the main
    file, Driver.cpp
    

  ADDITIONAL PARAMETERS:
    L0 = initial level discritization to use to begin 
      (i.e. for L0 = 3, the initial MLMC level will have
      a discritization of 2^3 segments)


  ORIGINAL COMMENT (Mike Giles):

    P = mlmc(Lmin,Lmax,N0,eps, mlmc_l, alpha,beta,gamma, Nl)

    multilevel Monte Carlo control routine

    Lmin  = minimum level of refinement       >= 2
    Lmax  = maximum level of refinement       >= Lmin
    N0    = initial number of samples         > 0
    eps   = desired accuracy (rms error)      > 0

    alpha -> weak error is  O(2^{-alpha*l})
    beta  -> variance is    O(2^{-beta*l})
    gamma -> sample cost is O(2^{gamma*l})    > 0

    if alpha, beta are not positive then they will be estimated

    mlmc_l(l,N,sums)   low-level function
    l       = level
    N       = number of paths
    sums[0] = sum(Y)
    sums[1] = sum(Y.^2)
    where Y are iid samples with expected value:
    E[P_0]           on level 0
    E[P_l - P_{l-1}] on level l>0

    P     = value
    Nl    = number of samples at each level

*//////////////////////////////////////////////////////////////
#pragma once
#include "mlmc.h"


double mlmc(int L0, int Lmin, int Lmax, int N0, double eps,
  std::function<void(int, int, double*)> mlmc_l,
  double alpha_0, double beta_0, double gamma, int *Nl){

  double sums[6], suml[3][21];
  double  ml[21], Vl[21], Cl[21], x[21], y[21], alpha, beta, sum, theta;
  int    dNl[21], L, converged;

  //
  // check input parameters
  //

  if (Lmin<2) {
    fprintf(stderr, "error: needs Lmin >= 2 \n");
    exit(1);
  }
  if (Lmax<Lmin) {
    fprintf(stderr, "error: needs Lmax >= Lmin \n");
    exit(1);
  }

  if (N0 <= 0 || eps <= 0.0 || gamma <= 0.0) {
    fprintf(stderr, "error: needs N>0, eps>0, gamma>0 \n");
    exit(1);
  }

  //
  // initialization
  //

  //alpha = fmax(0.0, alpha_0);
  //beta = fmax(0.0, beta_0);
  alpha = 0.0;
  beta = 0.0;
  theta = 0.25;             // MSE split between bias^2 and variance

  L = Lmin;
  converged = 0;

  for (int l = 0; l <= Lmax; l++) {
    Nl[l] = 0;
    //Cl[l] = pow(2.0, (double)l*gamma);
    Cl[l] = pow(2.0, (double)((l+L0)*gamma));

    for (int n = 0; n<3; n++) 
      suml[n][l] = 0.0;
  }

  for (int l = 0; l <= Lmin; l++) 
    dNl[l] = N0;

  //
  // main loop
  //

  while (!converged) {

    //
    // update sample sums
    //

    for (int l = 0; l <= L; l++) {
      if (dNl[l]>0) {
        mlmc_l(l, dNl[l], sums);
        suml[0][l] += (double)dNl[l];
        suml[1][l] += sums[0];
        suml[2][l] += sums[1];
      }
    }

    //
    // compute absolute average and variance,
    // correct for possible under-sampling,
    // and set optimal number of new samples
    //

    sum = 0.0;

    for (int l = 0; l <= L; l++) {
      ml[l] = fabs(suml[1][l] / suml[0][l]);
      Vl[l] = fmax(suml[2][l] / suml[0][l] - ml[l] * ml[l], 0.0);
      if (l>1) {
        ml[l] = fmax(ml[l], 0.5f*ml[l - 1] / pow(2.0, alpha));
        Vl[l] = fmax(Vl[l], 0.5f*Vl[l - 1] / pow(2.0, beta));
      }

      sum += sqrt(Vl[l] * Cl[l]);
    }

    for (int l = 0; l <= L; l++) {
      dNl[l] = ceil(fmax(0.0,
        sqrt(Vl[l] / Cl[l])*sum / ((1.0 - theta)*eps*eps)
        - suml[0][l]));
    }

    //
    // use linear regression to estimate alpha, beta if not given
    //

    if (alpha_0 <= 0.0) {
      for (int l = 1; l <= L; l++) {
        x[l - 1] = l;
        y[l - 1] = -log2(ml[l]);
      }
      regression(L, x, y, alpha, sum);
    }

    if (beta_0 <= 0.0) {
      for (int l = 1; l <= L; l++) {
        x[l - 1] = l;
        y[l - 1] = -log2(Vl[l]);
      }
      regression(L, x, y, beta, sum);
    }

    //
    // if (almost) converged, estimate remaining error and decide 
    // whether a new level is required
    //

    sum = 0.0;
    for (int l = 0; l <= L; l++)
      sum += fmax(0.0, (double)dNl[l] - 0.01*suml[0][l]);

    if (sum == 0) {
      converged = 1;
      double rem = ml[L] / (pow(2.0, gamma) - 1.0);

      if (rem > sqrt(theta)*eps) {
        if (L == Lmax)
          printf("*** failed to achieve weak convergence *** \n");
        else {
          converged = 0;
          L++;

          sum = 0.0;
          for (int l = 0; l <= L; l++) sum += sqrt(Vl[l] * Cl[l]);
          for (int l = 0; l <= L; l++)
            dNl[l] = ceil(fmax(0.0,
              sqrt(Vl[l] / Cl[l])*sum / ((1.0 - theta)*eps*eps)
              - suml[0][l]));
        }
      }
    }
  }

  //
  // finally, evaluate multilevel estimator
  //

  double P = 0.0;
  for (int l = 0; l <= L; l++) {
    P += suml[1][l] / suml[0][l];
    Nl[l] = suml[0][l];
  }

  return P;
}



//
// linear regression routine
//

void regression(int N, double *x, double *y, double &a, double &b) {

  double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sumy0 = 0.0, sumy1 = 0.0;

  for (int i = 1; i<N; i++) {
    sum0 += 1.0;
    sum1 += x[i];
    sum2 += x[i] * x[i];

    sumy0 += y[i];
    sumy1 += y[i] * x[i];
  }

  a = (sum0*sumy1 - sum1*sumy0) / (sum0*sum2 - sum1*sum1);
  b = (sum2*sumy0 - sum1*sumy1) / (sum0*sum2 - sum1*sum1);
}