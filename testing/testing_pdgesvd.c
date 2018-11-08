/**
 *
 * Copyright (c) 2017, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 **/

/**
 *
 * @file testing_pdgesvd.c
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
#include "ksvd.h"

/* Default values of parameters */
int nprow         = 1;
int npcol         = 1;
int lvec          = 1;
int rvec          = 1;
int n             = 5120;
int nb            = 128;
int mode          = 4;
double cond       = 9.0072e+15;
int optcond       = 0;
int start         = 5120;
int stop          = 5120;
int step          = 1;
int niter         = 1;
int polarqdwh     = 0;
int polarzolopd   = 0;
int polarsvd      = 0;
int slsvd         = 0;
int ksvdmr          = 0;
int ksvddc          = 0;
int ksvdel          = 0;
int check         = 0;
int verbose       = 0;


void print_usage(void)
{
    fprintf(stderr,
            "======= KSVD testing using ScaLAPACK\n"
            " -p      --nprow         : Number of MPI process rows\n"
            " -q      --npcol         : Number of MPI process cols\n"
            " -jl     --lvec          : Compute left singular vectors\n"
            " -jr     --rvec          : Compute right singular vectors\n"
            " -n      --N             : Dimension of the matrix\n"
            " -b      --nb            : Block size\n"
            " -m      --mode          : [1:6] Mode from pdlatms used to generate the matrix\n"
            " -k      --cond          : Condition number used to generate the matrix\n"
            " -o      --optcond       : Estimate Condition number using QR\n"
            " -i      --niter             : Number of iterations\n"
            " -r      --n_range           : Range for matrix sizes Start:Stop:Step\n"
            " -polarqdwh --polarqdwh      : Find polar decomposition using QDWH A=UH \n"
            " -polarzolopd --polarzolopd  : Find polar decomposition using ZOLO-PD A=UH \n"
            " -polarsvd  --polarsvd       : Find the polar decomposition using scalapack-svd \n"
            " -s      --slsvd             : Run reference ScaLAPACK SVD\n"
            " -w      --ksvdmr            : Run KSVD with ScaLAPACK MRRR EIG\n"
            " -e      --ksvddc            : Run KSVD with ScaLAPACK DC EIG\n"
            " -l      --ksvdel            : Run KSVD with ScaLAPACK DC EIG\n"
            " -c      --check             : Check the solution\n"
            " -fksvd --profksvd           : Enable profiling KSVD\n"
            " -v      --verbose           : Verbose\n"
            " -h      --help              : Print this help\n" );
}

#define GETOPT_STRING "p:q:x:y:n:b:m:i:o:r:Q,Z,S:s:w:e:c:f:t:v:h"

static struct option long_options[] =
    {
        /* PaRSEC specific options */
        {"nprow",      required_argument,  0, 'p'},
        {"npcol",      required_argument,  0, 'q'},
        {"jl",         no_argument,        0, 'x'},
        {"lvec",       no_argument,        0, 'x'},
        {"jr",         no_argument,        0, 'y'},
        {"rvec",       no_argument,        0, 'y'},
        {"N",          required_argument,  0, 'n'},
        {"n",          required_argument,  0, 'n'},
        {"nb",         required_argument,  0, 'b'},
        {"b",          required_argument,  0, 'b'},
        {"mode",       required_argument,  0, 'm'},
        {"m",          required_argument,  0, 'm'},
        {"cond",       required_argument,  0, 'k'},
        {"k",          required_argument,  0, 'k'},
        {"optcond",    required_argument,  0, 'o'},
        {"o",          required_argument,  0, 'o'},
        {"i",          required_argument,  0, 'i'},
        {"niter",      required_argument,  0, 'i'},
        {"r",          required_argument,  0, 'r'},
        {"n_range",    required_argument,  0, 'r'},
        //{"polar",      no_argument,        0, 'u'},
        {"polarqdwh",   no_argument,        0, 'Q'},
        {"polarzolopd", no_argument,        0, 'Z'},
        {"polarsvd",    no_argument,        0, 'S'},
        {"slsvd",       no_argument,        0, 's'},
        {"ksvdmr",        no_argument,        0, 'w'},
        {"ksvddc",        no_argument,        0, 'e'},
        {"ksvdel",        no_argument,        0, 'l'},
        {"e",           no_argument,        0, 'c'},
        {"check",       no_argument,        0, 'c'},
        {"profksvd",   no_argument,        0, 'f'},
        {"fksvd",      no_argument,        0, 'f'},
        {"verbose",     no_argument,        0, 'v'},
        {"help",        no_argument,        0, 'h'},
        {"h",           no_argument,        0, 'h'},
        {0, 0, 0, 0}
    };

static void parse_arguments(int argc, char** argv)
{
    int opt = 0;
    int c;
    int myrank_mpi;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);

    do {
#if defined(HAVE_GETOPT_LONG)
        c = getopt_long_only(argc, argv, "",
                        long_options, &opt);
#else
        c = getopt(argc, argv, GETOPT_STRING);
        (void) opt;
#endif  /* defined(HAVE_GETOPT_LONG) */

        switch(c) {
        case 'p': nprow     = atoi(optarg); break;
        case 'q': npcol     = atoi(optarg); break;
        case 'n': n         = atoi(optarg); start = n; stop = n; step = 1; break;
        case 'b': nb        = atoi(optarg); break;
        case 'm': mode      = atoi(optarg); break;
        case 'k': cond      = atof(optarg); break;
        case 'o': optcond   = atof(optarg); break;
        case 'Q': polarqdwh   = 1; break;
        case 'Z': polarzolopd = 1; break;
        case 'S': polarsvd    = 1; break;
        case 's': slsvd       = 1; break;
        case 'w': ksvdmr        = 1; break;
        case 'e': ksvddc        = 1; break;
        case 'l': ksvdel        = 1; break;
        case 'i': niter       = atoi(optarg); break;
        case 'r': get_range( optarg, &start, &stop, &step ); break;
        case 'c': check     = 1; break;
        case 'v': verbose   = 1; break;
        case 'h':
            if (myrank_mpi == 0) print_usage(); MPI_Finalize(); exit(0);
            break;
        default:
            break;
        }
    } while(-1 != c);
}

int main(int argc, char **argv) {

    int myrank_mpi, nprocs_mpi;
    int ictxt, myrow, mycol;
    int mloc, nloc, mlocW;
    int mpi_comm_rows, mpi_comm_cols;
    int useQr, THIS_REAL_ELPA_KERNEL_API;
    int i, j, k, iter, size, info_facto, info, iseed;
    int my_info_facto;
    int i0 = 0, i1 = 1;
    int lwork, liwork, *iWloc=NULL, ldw;
    long int LDW;
    int descA[9], descAcpy[9], descU[9], descVT[9], descH[9], descSigma[9];

    double *A=NULL, *Acpy=NULL, *U=NULL, *VT=NULL, *S=NULL, *Wloc=NULL, *D=NULL, *C=NULL;
    double *H=NULL;
    double *Sigma=NULL;


    double eps = LAPACKE_dlamch_work('e');
    int iprepad, ipostpad, sizemqrleft, sizemqrright, sizeqrf, sizeqtq,
        sizechk, sizesyevx, isizesyevx,
        sizesubtst, isizesubtst, sizetst,
        isizetst;

    double flops, GFLOPS;

    double berr = 0.0, my_berr = 0.0, norm_sv;
    double my_berr_ksvdmr = 0.0, my_berr_ksvddc = 0.0, my_berr_ksvdel = 0.0;
    double my_acc_ksvdmr = 0.0, my_acc_ksvddc = 0.0, my_acc_ksvdel = 0.0, my_acc_slsvd = 0.0;
    double my_orthR_ksvdmr = 0.0, my_orthR_ksvddc = 0.0, my_orthR_ksvdel = 0.0, my_orthR_slsvd = 0.0;
    double my_orthL_ksvdmr = 0.0, my_orthL_ksvddc = 0.0, my_orthL_ksvdel = 0.0, my_orthL_slsvd = 0.0;

    double orth_Usvd, berr_UHsvd, frobA;

    int success;
    double alpha, beta;
    char *jobu, *jobvt, *eigtype;
    int vl, vu, il, iu, nbeigvals, nbeigvecs;

    jobu  = lvec ? "V" : "N";
    jobvt = rvec ? "V" : "N";


/**/

    if (verbose & myrank_mpi == 0) fprintf(stderr, "Program starts... \n");

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

    if (verbose & myrank_mpi == 0) fprintf(stderr, "MPI Init done\n");
    parse_arguments(argc, argv);
    if (verbose & myrank_mpi == 0) fprintf(stderr, "Checking arguments done\n");

    Cblacs_get( -1, 0, &ictxt );
    Cblacs_gridinit( &ictxt, "R", nprow, npcol );
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
    if (myrank_mpi == 0) printf("\n =================== nprow %d npcol %d \n", nprow, npcol);
    if (verbose & myrank_mpi == 0) fprintf(stderr, "BLACS Init done\n");

    if (myrank_mpi == 0) {
        fprintf(stderr, "# \n");
        fprintf(stderr, "# NPROCS %d P %d Q %d\n", nprocs_mpi, nprow, npcol);
        fprintf(stderr, "# niter %d\n", niter);
        fprintf(stderr, "# n_range %d:%d:%d mode: %d cond: %2.4e \n", start, stop, step, mode, cond);
        fprintf(stderr, "# \n");
    }

    /* to run only the polar decompsition */
    if(polarsvd )  {slsvd = 0;} 

    if ( ksvdmr )
       eigtype = "r";
    else if (ksvddc)
       eigtype = "d";
    else if (ksvdel)
       eigtype = "e";

    if (verbose & myrank_mpi == 0) fprintf(stderr, "Range loop starts\n");

    if (polarzolopd){
        setPolarAlgorithm( KSVD_ZOLOPD );
    }

    // Begin loop over range
    for (size = start; size <= stop; size += step) {
        while ( (int)((double)size / (double)nb) < ( max(nprow , npcol) )){
            if (myrank_mpi == 0) fprintf(stderr, " Matrix size is small to be facrorized using this number of processors \n");
            size += step;
        }
        n = size; ldw = 2*n, LDW = ldw*n; long int matsize = n*n;

	mloc  = numroc_( &n, &nb, &myrow, &i0, &nprow );
	nloc  = numroc_( &n, &nb, &mycol, &i0, &npcol );
	mlocW = numroc_( &ldw, &nb, &myrow, &i0, &nprow );

        if (verbose & myrank_mpi == 0) fprintf(stderr, "Desc Init starts %d\n", mloc);
	descinit_( descA, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	descinit_( descAcpy, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	descinit_( descH, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	descinit_( descSigma, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );

	descinit_( descU, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
	descinit_( descVT, &n, &n, &nb, &nb, &i0, &i0, &ictxt, &mloc, &info );
        if (verbose & myrank_mpi == 0) fprintf(stderr, "Desc Init ends %d\n", mloc);

	A     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	Acpy  = (double *)malloc(mloc*nloc*sizeof(double)) ;
	H     = (double *)malloc(mloc*nloc*sizeof(double)) ;
	D  = (double *)malloc(n*sizeof(double)) ;
	Sigma = (double *)calloc(mloc*nloc,sizeof(double)) ;


	U      = (double *)malloc(mloc*nloc*sizeof(double)) ;
	VT      = (double *)malloc(mloc*nloc*sizeof(double)) ;
	S      = (double *)malloc(n*sizeof(double)) ;

        
        /* Generate matrix by pdlatms */
        {
           char   *dist = "N"; /* NORMAL( 0, 1 )  ( 'N' for normal ) */
           int    iseed[4] = {1, 0, 0, 1};
           char   *sym = "P"; /* The generated matrix is symmetric, with
                                eigenvalues (= singular values) specified by D, COND,
                                MODE, and DMAX; they will not be negative.
                                "N" not supported. */
           //int    mode = 4; /* sets D(i)=1 - (i-1)/(N-1)*(1 - 1/COND) */
           //double cond = 1.0/eps;
           double dmax = 1.0;
           int    kl   = n;
           int    ku   = n;
           char   *pack = "N"; /* no packing */
           int    order = n;
           int    info;
         
           pdlasizesep_( descA, 
                         &iprepad, &ipostpad, &sizemqrleft, &sizemqrright, &sizeqrf, 
                         &lwork, 
                         &sizeqtq, &sizechk, &sizesyevx, &isizesyevx, &sizesubtst, 
                         &isizesubtst, &sizetst, &isizetst );
           if (verbose & myrank_mpi == 0) fprintf(stderr, "Setting lwork done\n");
           Wloc = (double *)calloc(lwork,sizeof(double)) ;

           pdlatms_(&n, &n, dist,
                    iseed, sym, D, &mode, &cond, &dmax,
                    &kl, &ku, pack, 
                    A, &i1, &i1, descA, &order, 
                    Wloc, &lwork, &info);
           if (verbose & myrank_mpi == 0) fprintf(stderr, "MatGen done\n");
           if (info != 0) {
               fprintf(stderr, "An error occured during matrix generation: %d\n", info );
                   return EXIT_FAILURE;
           }
           pdlacpy_( "All", &n, &n, 
                     A, &i1, &i1, descA, 
                     Acpy, &i1, &i1, descAcpy ); 
           frobA  = pdlange_ ( "f", &n, &n, A, &i1, &i1, descA, Wloc);
           beta = 0.0;
           pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);

           for (i = 1; i <= n; i++) {
               int idum1, idum2, iloc, jloc;
               if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                     &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                           iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                           jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                           Sigma[ (jloc-1)*mloc + (iloc-1) ] = D[i-1];
               }
           } 

           norm_sv    = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc);
           if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy to Acpy done\n");

           free( Wloc );
        }

        if (myrank_mpi == 0) fprintf(stderr, "\n\n");
        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

        // QDWH + EIG
        if ( ksvdmr || ksvddc || ksvdel ) {
            /*
             * SVD decomposition is (A *U)* S * U' = C * S * U'.
             */
            // loop over iteration
            for (iter = 0; iter < niter; iter++) {
               flops = 0.0;

               if( (ksvdmr || ksvddc || ksvdel) ){
                   pdlacpy_( "A", &n, &n, A, &i1, &i1, descA, Acpy, &i1, &i1, descAcpy );
               }

               /*
                * Find the SVD using QDWH + EIG 
                */
               if (verbose & myrank_mpi == 0) fprintf(stderr, "EIG starts...\n");
              /*
               * Find Workspace 
               */
               lwork  = -1; liwork = -1;
               Wloc   = (double *)calloc(1,sizeof(double));
               iWloc  = (int *)calloc(1,sizeof(int));
               pdgeqsvd( jobu, jobvt, eigtype,  
                          n, n, 
                          A, i1, i1, descA, 
                          S, 
                          U,     i1,     i1,  descU,
                          VT,    i1,     i1,  descVT,
                          Wloc,  lwork,
                          iWloc, liwork, &my_info_facto);
               lwork  = (int)Wloc[0];
               liwork = (int)iWloc[0];
	       Wloc   = (double *)calloc(lwork,sizeof(double)) ;
	       iWloc  = (int *)calloc(liwork,sizeof(int)) ;
               /*
                * QDWH + EIG
                */

               pdgeqsvd( jobu, jobvt, eigtype,
                          n, n, 
                          A, i1, i1, descA, 
                          S, 
                          U,     i1,     i1,  descU,
                          VT,    i1,     i1,  descVT,
                          Wloc,  lwork,
                          iWloc, liwork, &my_info_facto);

               if (verbose & myrank_mpi == 0) fprintf(stderr, "Compute left singular vectors end...\n");
  
	       MPI_Allreduce( &my_info_facto, &info_facto, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

               if (verbose & myrank_mpi == 0) fprintf(stderr, "\nQDWH + ScaLAPACK EIG done\n");
         
              /*
               * Checking the SVD decomposition
               */
               if (check ) {
                    if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing starts...\n");
                    alpha = 1.0; beta = 0.0;
                   /*
                    * Set the singular values on the main diagonal 
                    * |A - U*Sigma*V'|
                    */
                    pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);

                    for (i = 1; i <= n; i++) {
                          int idum1, idum2, iloc, jloc;
                          if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                          &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                  iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                  jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                  Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                          }
                    } 


                    pdlacpy_( "All", &n, &n, 
                           Acpy, &i1, &i1, descAcpy, 
                           H, &i1, &i1, descH ); 
                    alpha = 1.0; beta = 0.0;
                    pdgemm_( "N", "N", &n, &n, &n, 
                           &alpha, 
                           U   , &i1, &i1, descU, 
                           Sigma, &i1, &i1, descSigma, 
                           &beta, 
                           A, &i1, &i1, descA);
                    beta = -1.0;
                    pdgemm_( "N", "T", &n, &n, &n, 
                           &alpha, 
                           A, &i1, &i1, descA, 
                           VT, &i1, &i1, descVT, 
                           &beta, 
                           H, &i1, &i1, descH);
                    my_berr_ksvdmr = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc) / (frobA * n);
                        
                   /* 
                    * Accuracy of singular values 
                    */
                    //for(i=0; i < n ; i++ )
                    //    D[i] = fabs(D[i]);
                    dlasrt_( "D", &n, S, &info );
                    dlasrt_( "D", &n, D, &info );
                    for(i=0; i < n ; i++ )
                          S[i] = S[i] - D[i];
                    for (i = 1; i <= n; i++) {
                          int idum1, idum2, iloc, jloc;
                          if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                          &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                               iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                               jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                               Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                          }
                    } 
                    my_acc_ksvdmr = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc) / norm_sv;

                   /* 
                    * Orthogonality of Left singular vectors 
                    */
                    alpha = 0.0; beta = 1.0;
                    pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                    alpha = 1.0; beta = -1.0;
                    pdgemm_( "T", "N", &n, &n, &n, 
                           &alpha, 
                           U, &i1, &i1, descU, 
                           U, &i1, &i1, descU, 
                           &beta, 
                           Sigma, &i1, &i1, descSigma);
                    my_orthL_ksvdmr = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                   /* 
                    * Orthogonality of Right singular vectors 
                    */
                    alpha = 0.0; beta = 1.0;
                    pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                    alpha = 1.0; beta = -1.0;
                    pdgemm_( "T", "N", &n, &n, &n, 
                           &alpha, 
                           VT, &i1, &i1, descVT, 
                           VT, &i1, &i1, descVT, 
                           &beta, 
                           Sigma, &i1, &i1, descSigma);
                    my_orthR_ksvdmr = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                    if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing ends...\n");
                    if (  myrank_mpi == 0) {
                          fprintf(stderr, "# QDWH + EIG\n"); 
                          fprintf(stderr, "#\n");
	                  fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ   \t\tinfo     \tAcc-sv    \tOrth-Rsv    \tOrth-Lsv     \tBerr  \n");
	                  fprintf(stderr, "  %6d \t%4d \t%3d \t%3d \t%3d", n, nb, nprocs_mpi, nprow, npcol);
	                  fprintf(stderr, "\t\t%d \t\t%2.4e \t%2.4e \t%2.4e \t%2.4e \n", 
                                       info_facto, my_acc_ksvdmr, my_orthR_ksvdmr, my_orthL_ksvdmr, my_berr_ksvdmr);
                    }
               }
               free(Wloc); free(iWloc);
               if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

               /*
                * Save copy of A in Acpy  
                */
               pdlacpy_( "All", &n, &n, 
                         Acpy, &i1, &i1, descAcpy, 
                         A, &i1, &i1, descA ); 
               if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy back to A done\n");
            }
            // loop over iteration
        }
        // QDWH + EIG

        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

       /*
        * ScaLAPACK SVD
        */
	if ( slsvd || polarsvd) {
            lwork = -1;
            Wloc  = (double *)calloc(1,sizeof(double));
            if( polarsvd){
               jobu  = "V";
               jobvt = "V";
            }
            pdgesvd_( jobu, jobvt, &n, &n, A, &i1, &i1, descA, 
                      S, 
                      U, &i1, &i1, descU, 
                      VT, &i1, &i1, descVT, 
                      Wloc, &lwork, &my_info_facto );
            lwork = (int)Wloc[0];
	    Wloc  = (double *)calloc(lwork,sizeof(double)) ;
 
            // loop over iteration
            for (iter = 0; iter < niter; iter++) {
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "\nScaLAPACK dgesvd starts...\n");
                  if( polarsvd){
                      jobu  = "V";
                      jobvt = "V";
                  }

                  pdlacpy_( "All", &n, &n, 
                            Acpy, &i1, &i1, descAcpy, 
                            A, &i1, &i1, descA ); 

                  pdgesvd_( jobu, jobvt, &n, &n, A, &i1, &i1, descA, 
                            S, 
                            U, &i1, &i1, descU, 
                            VT, &i1, &i1, descVT, 
                            Wloc, &lwork, &my_info_facto );
	          MPI_Allreduce( &my_info_facto, &info_facto, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "\nScaLAPACK dgesvd ends\n");

                  if (polarsvd && !slsvd){
                      /*
                       */ 
                       alpha = 1.0; beta = 0.0;
                       pdgemm_( "N", "N", &n, &n, &n, 
                                &alpha, 
                                U   , &i1, &i1, descU, 
                                VT,   &i1, &i1, descVT, 
                                &beta, 
                                A, &i1, &i1, descA);
                       alpha = 1.0; beta = 0.0;
                       pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);
                       for (i = 1; i <= n; i++) {
                           int idum1, idum2, iloc, jloc;
                           if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                           &&     ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                    iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                    jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                    Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                           }
                       } 
                       pdgemm_( "T", "N", &n, &n, &n, 
                                &alpha, 
                                VT   , &i1, &i1, descVT, 
                                Sigma,   &i1, &i1, descSigma, 
                                &beta, 
                                U, &i1, &i1, descU);
                       pdgemm_( "N", "N", &n, &n, &n, 
                                &alpha, 
                                U   , &i1, &i1, descU, 
                                VT,   &i1, &i1, descVT, 
                                &beta, 
                                Sigma, &i1, &i1, descSigma);
	               MPI_Allreduce( &my_info_facto, &info_facto, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                  }
                  /*
                   * Checking the polar factorization
                   */

                  if(polarsvd && check && !slsvd ){
                       /*
                        * checking orthogonality of Up
                        */ 
                        alpha = 0.0; beta = 1.0;
                        pdlaset_( "G", &n, &n, &alpha, &beta, H, &i1, &i1, descH);
                        alpha = 1.0; beta = -1.0;
                        pdgemm_( "T", "N", &n, &n, &n, 
                                 &alpha, 
                                 A, &i1, &i1, descA, 
                                 A, &i1, &i1, descA, 
                                 &beta, 
                                 H, &i1, &i1, descH);
                        orth_Usvd  = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc)/frobA;

                       /*
                        * checking the factorization |A-Up*H|
                        */ 
                        pdlacpy_( "A", &n, &n, Acpy, &i1, &i1, descAcpy, H, &i1, &i1, descH );
                        pdgemm_( "N", "N", &n, &n, &n, 
                                 &alpha, 
                                 A, &i1, &i1, descA, 
                                 Sigma, &i1, &i1, descSigma, 
                                 &beta, 
                                 H, &i1, &i1, descH);
                        berr_UHsvd  = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc)/frobA;
                        if ( polarsvd && myrank_mpi == 0) {
                            fprintf(stderr, "# Polar decomposition using ScaLAPACK DGESVD \n"); 
                            fprintf(stderr, "#\n");
	                    fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ   \tBerr_UpH  \tOrth_Up  \tinfo     \n");
	                    fprintf(stderr, "   %6d \t%4d \t%4d \t%3d \t%3d ", n, nb, nprocs_mpi, nprow, npcol);
	                    fprintf(stderr, "\t%2.4e \t%2.4e \t%d \n", berr_UHsvd, orth_Usvd,info_facto);
                            fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
                        }
                  }

                  if (check && slsvd && !polarsvd ) {
                     if ( verbose & myrank_mpi == 0) printf( "\n check Only on first iter %d \n", iter);
                     if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing starts...\n");
                     alpha = 1.0; beta = 0.0;
                     pdlaset_( "G", &n, &n, &beta, &beta, Sigma, &i1, &i1, descSigma);

                     for (i = 1; i <= n; i++) {
                             int idum1, idum2, iloc, jloc;
                             if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                             &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                     iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                     jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                     Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                             }

                     } 
                     pdgemm_( "N", "N", &n, &n, &n, 
                              &alpha, 
                              U   , &i1, &i1, descU, 
                              Sigma, &i1, &i1, descSigma, 
                              &beta, 
                              A, &i1, &i1, descA);
                     pdlacpy_( "All", &n, &n, 
                              Acpy, &i1, &i1, descAcpy, 
                              H, &i1, &i1, descH ); 
                     beta = -1.0;
                     pdgemm_( "N", "N", &n, &n, &n, 
                              &alpha, 
                              A, &i1, &i1, descA, 
                              VT, &i1, &i1, descVT, 
                              &beta, 
                              H, &i1, &i1, descH);
                     my_berr = pdlange_ ( "f", &n, &n, H, &i1, &i1, descH, Wloc) / (frobA * n);

                     /* Accuracy of singular values */
                     for(i=0; i < n ; i++ )
                          D[i] = fabs(D[i]);
                     dlasrt_( "D", &n, D, &info );
                     dlasrt_( "D", &n, S, &info );
                     for(i=0; i < n ; i++ )
                         S[i] = S[i] - D[i];
                     for (i = 1; i <= n; i++) {
                          int idum1, idum2, iloc, jloc;
                          if ( ( myrow == indxg2p_( &i, &nb, &idum1, &i0, &nprow ) )
                          &&   ( mycol == indxg2p_( &i, &nb, &idum1, &i0, &npcol ) ) ){
                                 iloc = indxg2l_( &i, &nb, &idum1, &idum2, &nprow );
                                 jloc = indxg2l_( &i, &nb, &idum1, &idum2, &npcol );
                                 Sigma[ (jloc-1)*mloc + (iloc-1) ] = S[i-1];
                          }
                     } 
                     my_acc_slsvd = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc) / norm_sv;

                     /* Orthogonality of Left singular vectors */
                     alpha = 0.0; beta = 1.0;
                     pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                     alpha = 1.0; beta = -1.0;
                     pdgemm_( "T", "N", &n, &n, &n, 
                             &alpha, 
                             U, &i1, &i1, descU, 
                             U, &i1, &i1, descU, 
                             &beta, 
                             Sigma, &i1, &i1, descSigma);
                     my_orthL_slsvd = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;

                     /* Orthogonality of Right singular vectors */
                     alpha = 0.0; beta = 1.0;
                     pdlaset_( "G", &n, &n, &alpha, &beta, Sigma, &i1, &i1, descSigma);
                     alpha = 1.0; beta = -1.0;
                     pdgemm_( "N", "T", &n, &n, &n, 
                             &alpha, 
                             VT, &i1, &i1, descVT, 
                             VT, &i1, &i1, descVT, 
                             &beta, 
                             Sigma, &i1, &i1, descSigma);
                     my_orthR_slsvd = pdlange_ ( "f", &n, &n, Sigma, &i1, &i1, descSigma, Wloc)/n;
	             if (slsvd && myrank_mpi == 0){
                         fprintf(stderr, "# ScaLAPACK DGESVD\n"); 
                         fprintf(stderr, "#\n");
	                 fprintf(stderr, "# \tN     \tNB   \tNP   \tP   \tQ    \t\tinfo     \tAcc-sv    \tOrth-Rsv    \tOrth-Lsv     \tBerr \n");
	                 fprintf(stderr, "   %6d \t%4d \t%3d \t%3d \t%3d", n, nb, nprocs_mpi, nprow, npcol);
                         fprintf(stderr, "\t\t%d \t\t%2.4e \t%2.4e \t%2.4e \t%2.4e\n", info_facto, my_acc_slsvd, my_orthR_slsvd, my_orthL_slsvd, my_berr/niter);
                     }

                     if (verbose & myrank_mpi == 0) fprintf(stderr, "Testing ends...\n");
                  } 
                  if (verbose & myrank_mpi == 0) fprintf(stderr, "Copy back to A done\n");
            }
            // loop over iteration
	    free( Wloc );
        }
        // ScaLAPACK SVD

        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");
        if (myrank_mpi == 0) fprintf(stderr, "/////////////////////////////////////////////////////////////////////////\n");

	free( Sigma );
	free( A );
	free( Acpy );
	free( H );
	free( D );

	free( U );
	free( VT );
	free( S );
        //if( ksvdmr || ksvddc || ksvdel )
        if (verbose & myrank_mpi == 0) fprintf(stderr, "Free matrices done\n");
    } // End loop over range


    if (verbose & myrank_mpi == 0) fprintf(stderr, "Range loop ends\n");

    blacs_gridexit_( &i0 );
    MPI_Finalize();
    if (verbose & myrank_mpi == 0) fprintf(stderr, "Program ends...\n");
    return 0;
}



