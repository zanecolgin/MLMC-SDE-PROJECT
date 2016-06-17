/*/////////////////////////////////////////////////////////////

  Gene Transcription Problem
    (Single Threaded CPU Version)
  using antithetic multilevel Monte Carlo


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

*//////////////////////////////////////////////////////////////

#pragma once


///////////////////////////////////////////////////////////////
// Definitions
///////////////////////////////////////////////////////////////

#define __OKAY__  0                                             /* Return values                                            */
#define __OUT_FILE_NAME1__ "../RESULTS/GENE-MIL-AMLMC.txt"      /* !!! file system dependent                                */
#define __OUT_FILE_NAME2__ "../RESULTS/GENE-MIL-AMLMC.csv"      /* !!! file system dependent                                */
#define __M__ 2                                                 /* !!! must be 2                                            */
#define __INITIAL_N__ 1000
#define __CONVERGENCE_N__ 10000                                 /* number of samples for convergence tests                  */
#define __T__ 1                                                 /*                                                          */
#define __G0__ 1.0                                              /* initial number of Gene molecules                         */
#define __S0__ { 10.0, 20.0, 0.0 }                              /* initial conditions for SDE                               */
#define __LEVEL_MIN__ 2                                         /* number of levels for epsilon tests                       */



#define __CONVERGENCE_L__     4                                 /* number of levels for convergence tests                   */
#define __LEVEL_MAX__         16                                /* number of levels for epsilon tests                       */
#define __FIRST_LEVEL__       4
#define __CONVERGENCE_EPS__ { 50.0, 20.0, 10.0, 5.0, 1.0, 0.0 } /* desired accuracy array (terminated by value 0)           */
