#ifndef __HDF5_TOOLS_H
#define __HDF5_TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>
#include "hdf5.h"
#include <hdf5_hl.h>
//#include "debug.h"

#define MAX_NAME 1024
float *hdf5_read_lf( const char* file_name, const char* dset_name,
size_t &W, size_t &H, size_t &S, size_t &T , size_t &C);
int hdf5_write_data_simple(const char* file_name, const char* dset_name,int* dimarray, int dimsize, float* data);


void do_dtype(hid_t);
void do_dset(hid_t);
void do_link(hid_t, char *);
void scan_group(hid_t);
void do_attr(hid_t);
void scan_attrs(hid_t);
void do_plist(hid_t);


#endif
