/*/////////////////////////////////////////////////////////////

  Multilevel Monte Carlo convergence tests
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
    outFile = 2nd file pointer for separtated output
    L0 = initial level discritization to use to begin 
      (i.e. for L0 = 3, the initial MLMC level will have
      a discritization of 2^3 segments)


  ORIGINAL COMMENT (Mike Giles):

    mlmc_test(mlmc_l, M, N,L, N0,Eps,Lmin,Lmax, fp)

    multilevel Monte Carlo test routine

    mlmc_l(l,N,sums)     low-level routine

    inputs:  l = level
    N = number of paths

    output: sums[0] = sum(Pf-Pc)
    sums[1] = sum((Pf-Pc).^2)
    sums[2] = sum((Pf-Pc).^3)
    sums[3] = sum((Pf-Pc).^4)
    sums[4] = sum(Pf)
    sums[5] = sum(Pf.^2)

    M      = refinement cost factor
    N      = number of samples for convergence tests
    L      = number of levels for convergence tests
    N0     = initial number of samples
    Eps    = desired accuracy array (terminated by value 0)
    Lmin   = minimum level of refinement
    Lmax   = maximum level of refinement

*//////////////////////////////////////////////////////////////
#pragma once
#include "mlmc_test.h"


// https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
// variadic macro to print to both file and stdout
#define PRINTF2(fp, ...) {printf(__VA_ARGS__);fprintf(fp,__VA_ARGS__);}

void mlmc_test(function<void(int, int, double*)> mlmc_l, int M, int N, int L,
  int N0, double *Eps, int Lmin, int Lmax, FILE *fp, 
  std::ofstream* outFile, int L0) {

  //
  // first, convergence tests
  //

  PRINTF2(fp, "\n");
  PRINTF2(fp, "**********************************************************\n");
  PRINTF2(fp, "*** Convergence tests, kurtosis, telescoping sum check ***\n");
  PRINTF2(fp, "**********************************************************\n");
  PRINTF2(fp, "\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)");
  PRINTF2(fp, "    kurtosis     check \n-------------------------");
  PRINTF2(fp, "--------------------------------------------------\n");

  *outFile
    << "l,"
    << "avg(Pf-Pc),"
    << "avg(Pf),"
    << "var(Pf-Pc),"
    << "var(Pf)"
    << std::endl;

  double sums[6];
  double *del1 = (double *)malloc((L + 1)*sizeof(double));
  double *del2 = (double *)malloc((L + 1)*sizeof(double));
  double *var1 = (double *)malloc((L + 1)*sizeof(double));
  double *var2 = (double *)malloc((L + 1)*sizeof(double));
  double *chk1 = (double *)malloc((L + 1)*sizeof(double));
  double *kur1 = (double *)malloc((L + 1)*sizeof(double));

  for (int l = 0; l <= L; l++) {
    mlmc_l(l, N, sums);

    for (int m = 0; m<6; m++) sums[m] = sums[m] / N;

    del1[l] = sums[0];
    del2[l] = sums[4];    
    //var1[l] = fmax(sums[1] - sums[0] * sums[0], 1e-10);       // !!! changed to support double precision floats
    //var2[l] = fmax(sums[5] - sums[4] * sums[4], 1e-10);       // !!! changing to 1e-20 may need adjustment
    var1[l] = fmax(sums[1] - sums[0] * sums[0], 1e-20);         // !!!
    var2[l] = fmax(sums[5] - sums[4] * sums[4], 1e-20);         // !!!

    kur1[l] = (sums[3]
      - 4.0*sums[2] * sums[0]
      + 6.0*sums[1] * sums[0] * sums[0]
      - 3.0*sums[0] * sums[0] * sums[0] * sums[0])
      / (var1[l] * var1[l]);

    if (l == 0)
      chk1[l] = 0.0;
    else
      chk1[l] = sqrt((double)N) *
      fabs(del1[l] + del2[l - 1] - del2[l])
      / (3.0*(sqrt(var1[l]) + sqrt(var2[l - 1]) + sqrt(var2[l])));

    PRINTF2(fp, "%2d   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n",
      l + L0, del1[l], del2[l], var1[l], var2[l], kur1[l], chk1[l]);
    //PRINTF2(fp, "%2d   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n",
    //  l, del1[l], del2[l], var1[l], var2[l], kur1[l], chk1[l]);
  
    *outFile
      << l + L0 << ","
      //<< l << ","
      << del1[l] << ","
      << del2[l] << ","
      << var1[l] << ","
      << var2[l]
      << std::endl;

  }

  //
  // print out a warning if kurtosis or consistency check looks bad
  //

  if (kur1[L] > 100.0) {
    PRINTF2(fp, "\n WARNING: kurtosis on finest level = %f \n", kur1[L]);
    PRINTF2(fp, " indicates MLMC correction dominated by a few rare paths; \n");
    PRINTF2(fp, " for information on the connection to variance of sample variances,\n");
    PRINTF2(fp, " see http://mathworld.wolfram.com/SampleVarianceDistribution.html \n");
  }

  double max_chk = 0.0;
  for (int l = 0; l <= L; l++) max_chk = fmax(max_chk, chk1[l]);
  if (max_chk > 1.0) {
    PRINTF2(fp, "\n WARNING: maximum consistency error = %f \n", max_chk);
    PRINTF2(fp, " indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied \n");
  }

  //
  // use linear regression to estimate alpha, beta
  //

  double alpha, beta, gamma, foo;
  double *x = (double *)malloc(L*sizeof(double));
  double *y = (double *)malloc(L*sizeof(double));

  for (int l = 1; l <= L; l++) {
    x[l - 1] = l;
    y[l - 1] = -log2(fabs(del1[l]));
  }
  regression(L, x, y, alpha, foo);

  for (int l = 1; l <= L; l++) {
    x[l - 1] = l;
    y[l - 1] = -log2(var1[l]);
  }
  regression(L, x, y, beta, foo);

  gamma = log2((double)M);

  PRINTF2(fp, "\n******************************************************\n");
  PRINTF2(fp, "*** Linear regression estimates of MLMC parameters ***\n");
  PRINTF2(fp, "******************************************************\n");
  PRINTF2(fp, "\n alpha = %f  (exponent for MLMC weak convergence)\n", alpha);
  PRINTF2(fp, " beta  = %f  (exponent for MLMC variance) \n", beta);
  PRINTF2(fp, " gamma = %f  (exponent for MLMC cost) \n", gamma);

  *outFile
    << "alpha,"
    << "beta"
    << std::endl;
  *outFile
    << alpha << ","
    << beta
    << std::endl;

  //
  // second, mlmc complexity tests
  //

  PRINTF2(fp, "\n");
  PRINTF2(fp, "***************************** \n");
  PRINTF2(fp, "*** MLMC complexity tests *** \n");
  PRINTF2(fp, "***************************** \n\n");
  PRINTF2(fp, "  eps   mlmc_cost   std_cost  savings     N_l \n");
  PRINTF2(fp, "----------------------------------------------- \n");
  *outFile
    << "eps,"
    << "mlmc_cost,"
    << "std_cost,"
    << "savings,"
    << "N_l"
    << std::endl;

  int i = 0;
  int *Nl = (int *)malloc((Lmax + 1)*sizeof(int));

  while (Eps[i]>0) {
    double eps = Eps[i++];

    mlmc(L0, Lmin, Lmax, N0, eps, mlmc_l, alpha, beta, gamma, Nl);

    double std_cost = 0.0, mlmc_cost = 0.0, theta = 0.25f;

    for (int l = 0; l <= Lmax; l++) {
      if (Nl[l]>0) {
        mlmc_cost += (1.0 + 1.0 / M)*Nl[l] * pow((double)M, (double)(l+L0));
        //mlmc_cost += (1.0 + 1.0 / M)*Nl[l] * pow((double)M, (double)l);
        if (l <= L) {
          std_cost += var2[l] * pow((double)M, (double)(l+L0));
          //std_cost += var2[l] * pow((double)M, (double)l);
        }
        else {          
          std_cost += var2[L] * pow((double)M, (double)(l+L0));
          //std_cost += var2[L] * pow((double)M, (double)l);
        }
      }
    }
    std_cost = std_cost / ((1.0 - theta)*eps*eps);

    PRINTF2(fp, "%.4f  %.3e  %.3e  %7.2f ",
      eps, mlmc_cost, std_cost, std_cost / mlmc_cost);
    *outFile
      << eps << ","
      << mlmc_cost << ","
      << std_cost << ","
      << std_cost / mlmc_cost << ",";
    for (int l = 0; Nl[l]>0; l++) {
      PRINTF2(fp, "%9d", Nl[l]);
      *outFile  << Nl[l] << ",";
    }
    PRINTF2(fp, "\n");    
    *outFile << std::endl;
  }
  PRINTF2(fp, "\n");  

}