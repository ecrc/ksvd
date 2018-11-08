Installation
============

Installation requires at least **CMake** of version 3.2.3. To build KSVD,
follow these instructions:

1.  Get KSVD from git repository

        git clone git@github.com:ecrc/ksvd

2.  Go into KSVD folder

        cd ksvd

3.  Create build directory and go there

        mkdir build && cd build

4.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/

5.  To build the testing binaries (optional)

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/ -DKSVD_TESTING:BOOL=ON 

6.  Use CMake to build KSVD based on existing installations of the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/ -DKSVD_TESTING:BOOL=ON -DPOLAR_DIR=/path/to/polar/installation -DSCALAPACK_DIR=/path/to/scalapack/installation -DSLTMG_LIBRARIES=/path/to/scalapack/installation/lib/libsltmg.a

7.  Build KSVD

        make -j

8.  Install KSVD

        make install

9. Add line

        export PKG_CONFIG_PATH=/path/to/install/lib/pkgconfig:$PKG_CONFIG_PATH

    to your .bashrc file.

Now you can use pkg-config executable to collect compiler and linker flags for
KSVD.

Testing and Timing
==================

The directories testing and timing contain an example 
to test the accuracy and the performance of KSVD using
ill/well-conditioned random matrices.

   The complete list of options is available below with -h option:
  
  ```
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
```
     On Cray systems, the launching command typically looks like:
    
       srun --ntasks=nT --hint=nomultithread ./main --nprow p --npcol q --b 64 --cond 1e16 --niter 1 --n_range start:stop:step --check --qwmr --qwdc --qwel --slsvd

     1. The number of the nodes is N, the number of tasks (nT) = N * (number_of_cores per node ). The programming model is pure MPI (no OpenMP, i.e., sequential BLAS).
     2. PxQ is the process grid configuration, where (nT - PxQ = 0)
     3. To compute the SVD decomposition using KSVD, the polar decomposition is calculated first, then followed by MRRR (--qwmr) or 
	 DC (--qwdc) or ELPA-DC (--qwel), as various alternatives for the symmetric eigensolvers.
     4. To use the regular bidiagonal reduction SVD from ScaLAPACK PDGESVD: --slsvd


