pipeline {
/*
 * Defining where to run
 */
//// Any:
// agent any
//// By agent label:
//      agent { label 'sandybridge' }

    agent { label 'jenkinsfile' }
    triggers {
        pollSCM('H/10 * * * *')
    }

    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }

    stages {
        stage ('build') {
            steps {
                sh '''#!/bin/bash -le
                    # loads modules
                    module purge
                    module load old-modules
                    module load intel/15
                    module load mpi-impi/5.0.1-intel-15
                    module load spack/git-morse-modules
                    module load cmake-3.5.2-gcc-4.8.5-6kea3zp
                    module load elpa/2015.11.001-intel-15

                    set -x
                    module list

                    # QDWH
                    wget -q "http://ecrcwiki.kaust.edu.sa:8080/job/ecrcrepo/job/qdwh-dev/job/master/lastSuccessfulBuild/artifact/build/QDWH-2.0.0-Linux.tar.gz"
                    tar -zxf QDWH-2.0.0-Linux.tar.gz
                    QDWHDIR=$PWD/QDWH-2.0.0-Linux

                    # build
                    mkdir -p build
                    cd build && rm -rf ./*
                    export I_MPI_CC="icc"
                    export I_MPI_F90="ifort"
                    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DKSVD_TESTING:BOOL=ON -DEXTRA_LIBS="ifcore" -DQDWH_DIR=$QDWHDIR

                    # build
                    make

                    # install
                    make install

                '''
            }
        }
        stage ('test') {
            steps {
                sh '''#!/bin/bash -le
                    # loads modules
                    module purge
                    module load old-modules
                    module load intel/15
                    module load mpi-impi/5.0.1-intel-15
                    module load elpa/2015.11.001-intel-15

                    set -x

                    module list

                    # Delete previous CTest results and run tests
                    rm -rf $WORKSPACE/build/Testing
                    cd $WORKSPACE/build
                    # need to fix the tests
                    sed -i s/:512:/:2048:/ testing/CTestTestfile.cmake
                    sed -i s/:512:/:2048:/ timing/CTestTestfile.cmake
                    export PATH=$PATH:. 
                    ctest --no-compress-output -T Test
                '''
            }
        }
    }
    // Post build actions
    post {
        //always {
        //}
        //success {
        //}
        //unstable {
        //}
        //failure {
        //}
        unstable {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
        failure {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
    }
}
