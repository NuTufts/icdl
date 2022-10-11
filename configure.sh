#!/bin/bash

# setup the environment variables for all the components

# save the folder where script is called
__icdl_configure_workdir__=$PWD

# set the basedir
export ICDL_BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $ICDL_BASEDIR

# larlite
cd larlite
source config/setup.sh
cd $ICDL_BASEDIR

# Geo2D
#cd Geo2D
#source config/setup.sh
#cd $ICDL_BASEDIR

# LArOpenCV
cd laropencv
source setup_laropencv.sh
cd $ICDL_BASEDIR

# # LArCV
cd larcv
source configure.sh
cd $ICDL_BASEDIR

# # Cilantro (3rd party)
cd cilantro
source dllee_setup.sh
cd $ICDL_BASEDIR

# # UB LArCV app
cd ublarcvapp
source configure.sh
cd $ICDL_BASEDIR

# # LArFlow
#cd larflow
#source configure.sh
#cd larmatchnet && source set_pythonpath.sh # setup larmatch environment variables
#cd $ICDL_BASEDIR

# # LArdly viewing tools
cd lardly
source setenv.sh
cd $ICDL_BASEDIR

cd ${__icdl_configure_workdir__}
