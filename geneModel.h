/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ABOUT

  FILE:
    geneModel.h

  AUTHOR:
    Zane Colgin

  VERSION:
    March 2016

  NOTE:
    "!!!" in comments indicates programming note
      i.e. look here for debugging notes, optimizations, etc.
    "???" in comments indicates missing information or 
      unfinished

  DESCRIPTION:  
    ???
    
  REFERENCES:
    [1] D. F. Anderson and D. J. Higham. Multilevel Monte Carlo
        for continuous time Markov chains, with applications in
        biochemical kinetics, Multiscale Model. Simul., 10
        (2012), pp. 146-179.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once


namespace geneModel {

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~ Declare Namespace Consants
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  static const int kNumReactions = 5;
  static const int kNumSpecies = 3;
  static const double kC[kNumReactions] = {
  // !!! IMPORTANT:  values in comments below
  // !!!   indicate value referenced in [1] page 169
  // !!!   DO NOT CHANGE
      25    , // c1 = 25
    1000    , // c2 = 1000
       0.001, // c3 = 0.001
       0.1  , // c4 = 0.1
       1      // c5 = 1
  };

};