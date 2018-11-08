/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file ksvd.h
 *
 *  KSVD is a high performance software framework for computing 
 *  a dense SVD on distributed-memory manycore systems provided by KAUST
 *
 * @version 2.0.0
 * @author Dalal Sukkari
 * @author Hatem Ltaief
 * @date 2018-11-08
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpi.h>
#include "myscalapack.h"
#include "flops.h"
#include <elpa/elpa.h>

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

enum POLAR_ALGORITHM { POLAR_ALGORITHM_START = 0, KSVD_QDWH = 1, KSVD_ZOLOPD = 2, POLAR_ALGORITHM_END };
int getPolarAlgorithm( );
int setPolarAlgorithm( int );

int pdgeqsvd( char *jobu, char *jobvt, char *eigtype, 
              int m, int n, 
              double *A, int iA, int jA, int *descA, 
              double *S, 
              double *U,     int iU,     int jU, int *descU,
              double *VT,    int iVT,    int jVT, int *descVT,
              double *Work,  int lWork,
              int    *iWork, int liWork, int *info);
