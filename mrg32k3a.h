/*/////////////////////////////////////////////////////////////

  Random Number Generator

  Zane Colgin
  Middle Tennessee State University
  January 2016

  DATE        AUTHOR  COMMENTS
  ----------  ------  ---------------
  ????-??-??  MG      Original version
  2016-01-30  ZC      Modified version


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
    Only slight modifications have been made.

*//////////////////////////////////////////////////////////////

#pragma once

#include <cstdio>
#include <cstdlib>
#include <math.h>
//#include <amp_math.h>

typedef unsigned int uint;
typedef unsigned long long ull;

ull CPU_mrg32k3a_v1[3], CPU_mrg32k3a_v2[3];

#define mrg32k3a_pi     3.14159265358979323846
#define mrg32k3a_norm   2.3283065498378288e-10

#define mrg32k3a_2pow32 4294967296

#define mrg32k3a_m1     4294967087
#define mrg32k3a_m1c           209
#define mrg32k3a_m2     4294944443
#define mrg32k3a_m2c         22853

#define mrg32k3a_a12                 1403580
#define mrg32k3a_a13   (4294967087 -  810728)
#define mrg32k3a_a13n                 810728


#define mrg32k3a_a21                  527612
#define mrg32k3a_a23   (4294944443 - 1370589)
#define mrg32k3a_a23n                1370589

#define NORMCDF(a) (0.5f*erfcf(-(a)/sqrtf(2.0f)))


//using namespace Concurrency::precise_math;

void CPU_mrg32k3a_init(uint *v1, uint *v2, int offset)
{
  ull A1[9], A2[9], A1b[9], A2b[9], a1, a2, v1b[3], v2b[3];
  int i, j, k, n;

  for (n = 0; n<3; n++) {
    CPU_mrg32k3a_v1[n] = v1[n];
    CPU_mrg32k3a_v2[n] = v2[n];
  }

  // initialise skip-ahead matrices

  A1[0] = 0;            A1[3] = 1;            A1[6] = 0;
  A1[1] = 0;            A1[4] = 0;            A1[7] = 1;
  A1[2] = mrg32k3a_a13; A1[5] = mrg32k3a_a12; A1[8] = 0;

  A2[0] = 0;            A2[3] = 1;            A2[6] = 0;
  A2[1] = 0;            A2[4] = 0;            A2[7] = 1;
  A2[2] = mrg32k3a_a23; A2[5] = 0;            A2[8] = mrg32k3a_a21;

  while (offset>0) {

    // apply offset to seed vectors

    if (offset % 2 == 1) {
      for (i = 0; i<3; i++) {
        v1b[i] = ((A1[i + 3 * 0] * CPU_mrg32k3a_v1[0]) % mrg32k3a_m1
          + (A1[i + 3 * 1] * CPU_mrg32k3a_v1[1]) % mrg32k3a_m1
          + (A1[i + 3 * 2] * CPU_mrg32k3a_v1[2]) % mrg32k3a_m1) % mrg32k3a_m1;
        v2b[i] = ((A2[i + 3 * 0] * CPU_mrg32k3a_v2[0]) % mrg32k3a_m2
          + (A2[i + 3 * 1] * CPU_mrg32k3a_v2[1]) % mrg32k3a_m2
          + (A2[i + 3 * 2] * CPU_mrg32k3a_v2[2]) % mrg32k3a_m2) % mrg32k3a_m2;
      }
      for (i = 0; i<3; i++) {
        CPU_mrg32k3a_v1[i] = v1b[i];
        CPU_mrg32k3a_v2[i] = v2b[i];
      }
    }
    offset = offset / 2;

    // square skip-ahead matrices

    for (i = 0; i<3; i++) {
      for (j = 0; j<3; j++) {
        a1 = 0;
        a2 = 0;
        for (k = 0; k<3; k++) {
          a1 += (A1[i + 3 * k] * A1[k + 3 * j]) % mrg32k3a_m1;
          a2 += (A2[i + 3 * k] * A2[k + 3 * j]) % mrg32k3a_m2;
        }
        A1b[i + 3 * j] = a1%mrg32k3a_m1;
        A2b[i + 3 * j] = a2%mrg32k3a_m2;
      }
    }
    for (i = 0; i<9; i++) {
      A1[i] = A1b[i];
      A2[i] = A2b[i];
    }
  }
}

//
// Note: CPU version faster using 64-bit % operations
//

inline void CPU_mrg32k3a_next_step(ull *v1, ull *v2, ull &p)
{
  p = mrg32k3a_a12*v1[1] + mrg32k3a_a13n*(mrg32k3a_m1 - v1[0]);

  // alternative 1: 64-bit modulus
  p = p % mrg32k3a_m1;

  // alternative 2: avoids 64-bit modulus
  // p = (p & (mrg32k3a_2pow32-1)) + (p>>32)*mrg32k3a_m1c;
  // if (p>=mrg32k3a_m1) p -= mrg32k3a_m1;

  v1[0] = v1[1]; v1[1] = v1[2]; v1[2] = p;

  p = mrg32k3a_a21*v2[2] + mrg32k3a_a23n*(mrg32k3a_m2 - v2[0]);

  // alternative 1: 64-bit modulus
  p = p % mrg32k3a_m2;  // replaced by next 3 lines

                        // alternative 2: avoids 64-bit modulus
                        // p = (p & (mrg32k3a_2pow32-1)) + (p>>32)*mrg32k3a_m2c;
                        // p = (p & (mrg32k3a_2pow32-1)) + (p>>32)*mrg32k3a_m2c;
                        // if (p>=mrg32k3a_m2) p -= mrg32k3a_m2;

  v2[0] = v2[1]; v2[1] = v2[2]; v2[2] = p;

  p = v1[2] - v2[2];
  if (v1[2]<v2[2]) p += mrg32k3a_m1;
}


template <class FP>
inline void CPU_mrg32k3a_next_uniform(ull *v1, ull *v2, FP &x)
{
  ull p;

  CPU_mrg32k3a_next_step(v1, v2, p);
  x = (FP)mrg32k3a_norm*(FP)p;
}


void CPU_mrg32k3a_uniform(int N, float *P)
{
  ull v1[3], v2[3], p;
  int n;

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n<N; n++) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    P[n] = (float)mrg32k3a_norm*(float)p;
  }
}


void CPU_mrg32k3a_uniformD(int N, double *P)
{
  ull v1[3], v2[3], p;
  int n;

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n<N; n++) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    P[n] = mrg32k3a_norm*p;
  }
}


void CPU_mrg32k3a_exponential(int N, float *P)
{
  ull v1[3], v2[3], p;
  int n;

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n<N; n++) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    P[n] = logf((float)mrg32k3a_norm*(float)p);
  }
}


void CPU_mrg32k3a_exponentialD(int N, double *P)
{
  ull v1[3], v2[3], p;
  int n;

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n<N; n++) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    P[n] = log(mrg32k3a_norm*p);
  }
}


//void sincos(float arg, float *s, float *c) {  sincosf(arg, s, c); }

template <class FP>
inline void CPU_mrg32k3a_next_normal(ull *v1, ull *v2, FP &x0, FP &x1)
{
  ull p;
  FP  x2;

  if (isnan(x1)) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    x2 = sqrt(-(FP)2.0*log((FP)mrg32k3a_norm*(FP)p));

    CPU_mrg32k3a_next_step(v1, v2, p);
    //sincos((FP)2.0*(FP)mrg32k3a_pi*(FP)mrg32k3a_norm*(FP)p, &x0, &x1);
    x0 = sin((FP)2.0*(FP)mrg32k3a_pi*(FP)mrg32k3a_norm*(FP)p);
    x1 = cos((FP)2.0*(FP)mrg32k3a_pi*(FP)mrg32k3a_norm*(FP)p);
    x0 = x0*x2;
    x1 = x1*x2;
  }
  else {
    x0 = x1;
    x1 = (FP)nanf("");
  }
}


void CPU_mrg32k3a_normal(int N, float *P)
{
  ull   v1[3], v2[3], p;
  int   n;
  float x2 = 0.0f;

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n <= N; n++) {
    CPU_mrg32k3a_next_step(v1, v2, p);

    if (n % 2 == 0) {
      x2 = sqrtf(-2.0f*logf((float)mrg32k3a_norm*p));
    }
    else {
      // sincosf(2.0f*mrg32k3a_pi*mrg32k3a_norm*p,&x0,&x1);
      P[n - 1] = x2*sinf(2.0f*(float)mrg32k3a_pi*(float)mrg32k3a_norm*p);
      if (n<N)
        P[n] = x2*cosf(2.0f*(float)mrg32k3a_pi*(float)mrg32k3a_norm*p);
    }
  }
}


void CPU_mrg32k3a_normalD(int N, double *P)
{
  ull   v1[3], v2[3], p;
  int   n;
  double x2 = 0.0f;

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n <= N; n++) {
    CPU_mrg32k3a_next_step(v1, v2, p);

    if (n % 2 == 0) {
      x2 = sqrt(-2.0*log(mrg32k3a_norm*p));
    }
    else {
      // sincosf(2.0f*mrg32k3a_pi*mrg32k3a_norm*p,&x0,&x1);
      P[n - 1] = x2*sin(2.0*mrg32k3a_pi*mrg32k3a_norm*p);
      if (n<N)
        P[n] = x2*cos(2.0*mrg32k3a_pi*mrg32k3a_norm*p);
    }
  }
}


inline float CPU_mrg32k3a__next_uniform(ull *v1, ull *v2)
{
  ull p;
  CPU_mrg32k3a_next_step(v1, v2, p);
  return (float)mrg32k3a_norm*(float)p;
}

inline float CPU_mrg32k3a__next_normal(ull *v1, ull *v2, float &x2)
{
  ull   p;
  float x3;

  if (isnan(x2)) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    x3 = sqrt(-2.0f*logf((float)mrg32k3a_norm*(float)p));
    CPU_mrg32k3a_next_step(v1, v2, p);
    x2 = x3*cosf(2.0f*(float)mrg32k3a_pi*(float)mrg32k3a_norm*(float)p);
    x3 = x3*sinf(2.0f*(float)mrg32k3a_pi*(float)mrg32k3a_norm*(float)p);
  }
  else {
    x3 = x2;
    x2 = nanf(" ");
  }
  return x3;
}


inline float CPU_mrg32k3a_next_gamma(ull *v1, ull *v2, float alpha, float &x2)
{
  float prefix, d, c, x, v, u;

  if (alpha <= 1.0f) {
    prefix = powf(CPU_mrg32k3a__next_uniform(v1, v2), 1.0f / alpha);
    alpha += 1.0f;
  }
  else {
    prefix = 1.0f;
  }

  d = alpha - 1.0f / 3.0f;
  c = 1.0f / sqrt(9.0f*d);

  do {
    x = CPU_mrg32k3a__next_normal(v1, v2, x2);
    v = 1.0f + c*x;
    v = v*v*v;
    u = CPU_mrg32k3a__next_uniform(v1, v2);
    // printf(" CPU: x, u = %f %f \n",x,u);
  } while ((v <= 0.0f) || log(u) >= 0.5f*x*x + d*(1.0f - v + log(v)));

  return prefix*d*v;
}


void CPU_mrg32k3a_gamma(int N, float alpha, float *P)
{
  ull   v1[3], v2[3];
  int   n;
  float x2 = nanf(" ");

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n<N; n++) {
    P[n] = CPU_mrg32k3a_next_gamma(v1, v2, alpha, x2);
  }
}


inline double CPU_mrg32k3a__next_uniformD(ull *v1, ull *v2)
{
  ull p;
  CPU_mrg32k3a_next_step(v1, v2, p);
  return (double)mrg32k3a_norm*(double)p;
}

inline double CPU_mrg32k3a__next_normalD(ull *v1, ull *v2, double &x2)
{
  ull    p;
  double x3;

  if (isnan(x2)) {
    CPU_mrg32k3a_next_step(v1, v2, p);
    x3 = sqrt(-2.0*log((double)mrg32k3a_norm*(double)p));
    CPU_mrg32k3a_next_step(v1, v2, p);
    x2 = x3*cos(2.0*(double)mrg32k3a_pi*mrg32k3a_norm*(double)p);
    x3 = x3*sin(2.0*(double)mrg32k3a_pi*mrg32k3a_norm*(double)p);
  }
  else {
    x3 = x2;
    x2 = nan(" ");
  }
  return x3;
}


inline double CPU_mrg32k3a_next_gammaD(ull *v1, ull *v2, double alpha, double &x2)
{
  double prefix, d, c, x, v, u;

  if (alpha <= 1.0) {
    prefix = pow(CPU_mrg32k3a__next_uniformD(v1, v2), 1.0 / alpha);
    alpha += 1.0;
  }
  else {
    prefix = 1.0;
  }

  d = alpha - 1.0 / 3.0;
  c = 1.0 / sqrt(9.0*d);

  do {
    x = CPU_mrg32k3a__next_normalD(v1, v2, x2);
    v = 1.0 + c*x;
    v = v*v*v;
    u = CPU_mrg32k3a__next_uniformD(v1, v2);
    // printf(" CPU: x, u = %f %f \n",x,u);
  } while ((v <= 0.0) || log(u) >= 0.5*x*x + d*(1.0 - v + log(v)));

  return prefix*d*v;
}


void CPU_mrg32k3a_gammaD(int N, double alpha, double *P)
{
  ull    v1[3], v2[3];
  int    n;
  double x2 = nanf(" ");

  for (n = 0; n<3; n++) {
    v1[n] = CPU_mrg32k3a_v1[n];
    v2[n] = CPU_mrg32k3a_v2[n];
  }

  for (n = 0; n<N; n++) {
    P[n] = CPU_mrg32k3a_next_gammaD(v1, v2, alpha, x2);
  }
}

