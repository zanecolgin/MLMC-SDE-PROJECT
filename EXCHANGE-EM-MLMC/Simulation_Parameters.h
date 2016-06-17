/*/////////////////////////////////////////////////////////////

  2-Assest exchange option
    (Single Threaded CPU Version)
  using Euler-Maruyama method and multilevel Monte Carlo


  Zane Colgin
  Middle Tennessee State University
  January 2016

  DATE        AUTHOR  COMMENTS
  ----------  ------  ---------------
  2016-01-30  ZC      Modified version


  NOTE:
    "!!!" in comments indicates programming note 
      i.e. look here for debugging notes, optimizations, etc.
    "???" in comments indicates missing information or 
      unfinished

*//////////////////////////////////////////////////////////////

#pragma once


///////////////////////////////////////////////////////////////
// Definitions
///////////////////////////////////////////////////////////////

#define __OKAY__  0                                             /* Return values                                            */
#define __OUT_FILE_NAME1__ "../RESULTS/EXCHANGE-EM-MLMC.txt"    /* !!! file system dependent                                */
#define __OUT_FILE_NAME2__ "../RESULTS/EXCHANGE-EM-MLMC.csv"    /* !!! file system dependent                                */
#define __M__ 2                                                 /* !!! must be 2                                            */
#define __INITIAL_N__ 100
#define __CONVERGENCE_N__ 10000                                 /* number of samples for convergence tests                  */
#define __CONVERGENCE_L__ 6                                     /* number of levels for convergence tests                   */
#define __LEVEL_MIN__ 2                                         /* number of levels for epsilon tests                       */
#define __LEVEL_MAX__ 14                                        /* number of levels for epsilon tests                       */
#define __CONVERGENCE_EPS__ { 5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 0.0 } /* accuracy (terminated by value 0)           */
#define __S1__ 100                                              /* Initial value                                            */

#define __FIRST_LEVEL__ 0