/*/////////////////////////////////////////////////////////////

  Stiff Problem
    (Single Threaded CPU Version)
  using drifting split-step backwards Milstein method and
    multilevel Monte Carlo


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
#define __OUT_FILE_NAME1__ "../RESULTS/STIFF-DSSBM-MLMC.txt"
#define __OUT_FILE_NAME2__ "../RESULTS/STIFF-DSSBM-MLMC.csv"
#define __M__ 2

#define __T__             1                                     /*                                                          */
#define __Y10__           0.5                                   /*                                                          */
#define __Y20__           0.5                                   /*                                                          */
#define __LEVEL_MIN__     2

#define __INITIAL_N__     1000
#define __CONVERGENCE_N__ 10000                                 /* number of samples for convergence tests                  */
#define __CONVERGENCE_L__ 7                                     /* number of levels for convergence tests                   */
#define __LEVEL_MAX__     16
#define __FIRST_LEVEL__   5
#define __CONVERGENCE_EPS__ { 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 0.0 } /* desired accuracy array (terminated by value 0)    */

                          // c1     c2
#define __TESTS__         { - 4.0, - 4.0},\
                          { - 8.0, - 8.0},\
                          { -12.0, -12.0},\
                          { -16.0, - 1.0 }
#define __N_TESTS_        4