

#ifndef __LIBPNG_H
#define __LIBPNG_H


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <png.h>

#define umaxof(t) (((0x1ULL << ((sizeof(t) * 8ULL)-1ULL)) - 1ULL) |\
                  (0xFULL <<((sizeof(t) * 8ULL) - 4ULL)))

// Creates a test image for saving. Creates a Mandelbrot Set fractal of size width x height
float *createMandelbrotImage(int width, int height, float xS, float yS, float rad, int maxIteration);

// This takes the float value 'val', converts it to red, green & blue values, then
// sets those values into the image memory buffer location pointed to by 'ptr'
inline void setRGB(png_byte *ptr, float val);

// This function actually writes out the PNG image file. The string 'title' is
// also written into the image file
int writeImageRGB(char* filename, int width, int height, float *buffer, char* title);
int writeImageGray(char* filename, int width, int height, float *buffer, char* title);

float sRGB_to_linear(float x);
float linear_to_sRGB(float x);
float * convertRGBtoGray(int width, int height,float* origin);



#endif
