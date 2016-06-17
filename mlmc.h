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

*//////////////////////////////////////////////////////////////

#pragma once

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <functional>

void regression(int, double *, double *, double &a, double &b);

double mlmc(int L0, int Lmin, int Lmax, int N0, double eps,
  std::function<void(int, int, double*)> mlmc_l,
  double alpha_0, double beta_0, double gamma, int *Nl);