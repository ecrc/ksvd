pipeline {
/*
 * Defining where to run
 */
//// Any:
// agent any
//// By agent label:
//      agent { label 'sandybridge' }
    // no agents, each stage must declare it
    agent none
    triggers {
        pollSCM('H/10 * * * *')
    }

    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }

    stages {
        stage ('polar-build') {
            agent { label "jenkinsfile" }
            steps {
                sh '''#!/bin/bash -le
                    # loads modules
                    module purge
                    module load cmake/3.9.6
                    module load intel/2017
                    module load intelmpi/2017/intel-2017
                    module load elpa/2015.11.001-intel-2017

                    set -x
                    module list

                    export CC=icc # just in case
                    export FC=ifort # just in case
                    export F90=ifort # just in case
                    export I_MPI_CC="$CC"
                    export I_MPI_FC="$FC"
                    export I_MPI_F90="$F90"

                    # POLAR
                    wget -q "http://ecrcwiki.kaust.edu.sa:8080/job/ecrcrepo/job/qdwh-dev/job/polar/lastSuccessfulBuild/artifact/build/POLAR-2.0.0-Linux.tar.gz" -O - | tar -zx
                    POLARDIR=$PWD/POLAR-2.0.0-Linux

                    # build
                    mkdir -p build
                    cd build && rm -rf ./*
                    export I_MPI_CC="icc"
                    export I_MPI_F90="ifort"
                    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DKSVD_TESTING:BOOL=ON -DEXTRA_LIBS="ifcore" -DPOLAR_DIR=$POLARDIR

                    # build
                    make

                    # install
                    make install

                '''
                stash name: "build-polar", includes: "build/**"
            }
        }
        stage ('polar-test') {
            agent { label "jenkinsfile" }
            steps {
                unstash 'build-polar'
                sh '''#!/bin/bash -le
                    # loads modules
                    module purge
                    module load cmake/3.9.6
                    module load intel/2017
                    module load intelmpi/2017/intel-2017
                    module load elpa/2015.11.001-intel-2017

                    set -x
                    module list

                    # Delete previous CTest results and run tests
                    rm -rf $WORKSPACE/build/Testing
                    cd $WORKSPACE/build
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
