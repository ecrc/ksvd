/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file common.h
 *
 *  KSVD is a high performance software framework for computing 
 *  a dense SVD on distributed-memory manycore systems provided by KAUST
 *
 * @version 1.0.0
 * @author Dalal Sukkari
 * @author Hatem Ltaief
 * @date 2017-11-13
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
