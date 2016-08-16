/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn                   */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef COLOR_LIB_INCLUDED
#define COLOR_LIB_INCLUDED

#include <math.h>
#include "funct_lib.c"

/*---------------------------------------------------------------------------*/

void vector_to_RGB 
(
             /****************************************************************/
  float x,   /* in  : x-component                                            */
  float y,   /* in  : y-component                                            */
  int   *R,  /* out : red component                                          */
  int   *G,  /* out : green component                                        */
  int   *B   /* out : blue component                                         */
             /****************************************************************/
 ) 
{
                     /********************************************************/
  float Pi;          /* pi                                                   */
  float amp;         /* amplitude (magnitude)                                */
  float phi;         /* phase (angle)                                        */
  float alpha, beta; /* weights for linear interpolation                     */
                     /********************************************************/

  /* set pi */
  Pi = 2.0 * acos(0.0);

  /* determine amplitude and phase (cut amp at 1) */
  amp = sqrt (x * x + y * y);
  if (amp > 1) amp = 1;
  if (x == 0.0)
    if (y >= 0.0) phi = 0.5 * Pi;
    else phi = 1.5 * Pi;
  else if (x > 0.0)
    if (y >= 0.0) phi = atan (y/x);
    else phi = 2.0 * Pi + atan (y/x);
  else phi = Pi + atan (y/x);
  
  phi = phi / 2.0;

  // interpolation between red (0) and blue (0.25 * Pi)
  if ((phi >= 0.0) && (phi < 0.125 * Pi)) {
    beta  = phi / (0.125 * Pi);
    alpha = 1.0 - beta;
    *R = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
    *G = (int)floor(amp * (alpha *   0.0 + beta *   0.0));
    *B = (int)floor(amp * (alpha *   0.0 + beta * 255.0));
  }
  if ((phi >= 0.125 * Pi) && (phi < 0.25 * Pi)) {
    beta  = (phi-0.125 * Pi) / (0.125 * Pi);
    alpha = 1.0 - beta;
    *R = (int)floor(amp * (alpha * 255.0 + beta *  64.0));
    *G = (int)floor(amp * (alpha *   0.0 + beta *  64.0));
    *B = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
  }
  // interpolation between blue (0.25 * Pi) and green (0.5 * Pi)
  if ((phi >= 0.25 * Pi) && (phi < 0.375 * Pi)) {
    beta  = (phi - 0.25 * Pi) / (0.125 * Pi);
    alpha = 1.0 - beta;
    *R = (int)floor(amp * (alpha *  64.0 + beta *   0.0));
    *G = (int)floor(amp * (alpha *  64.0 + beta * 255.0));
    *B = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
  }
  if ((phi >= 0.375 * Pi) && (phi < 0.5 * Pi)) {
    beta  = (phi - 0.375 * Pi) / (0.125 * Pi);
    alpha = 1.0 - beta;
    *R = (int)floor(amp * (alpha *   0.0 + beta *   0.0));
    *G = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
    *B = (int)floor(amp * (alpha * 255.0 + beta *   0.0));
  }
  // interpolation between green (0.5 * Pi) and yellow (0.75 * Pi)
  if ((phi >= 0.5 * Pi) && (phi < 0.75 * Pi)) {
    beta  = (phi - 0.5 * Pi) / (0.25 * Pi);
    alpha = 1.0 - beta;
    *R = (int)floor(amp * (alpha * 0.0   + beta * 255.0));
    *G = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
    *B = (int)floor(amp * (alpha * 0.0   + beta * 0.0));
  }
  // interpolation between yellow (0.75 * Pi) and red (Pi)
  if ((phi >= 0.75 * Pi) && (phi <= Pi)) {
    beta  = (phi - 0.75 * Pi) / (0.25 * Pi);
    alpha = 1.0 - beta;
    *R = (int)floor(amp * (alpha * 255.0 + beta * 255.0));
    *G = (int)floor(amp * (alpha * 255.0 + beta *   0.0));
    *B = (int)floor(amp * (alpha * 0.0   + beta *   0.0));
  }

  /* check RGB range */
  *R = byte_range(*R);
  *G = byte_range(*G);
  *B = byte_range(*B);
}

/*--------------------------------------------------------------------------*/

#endif
