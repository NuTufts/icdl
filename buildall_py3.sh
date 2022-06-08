#!/bin/bash

alias python=python3
alias python-config=python3-config

__icdl_buildall_py3_workdir__=$PWD
build_log=${__icdl_buildall_py3_workdir__}/build.log

echo "<<< BUILD LARLITE >>>"
cd larlite
mkdir build
cd build
cmake -DUSE_PYTHON3=ON ../
make install -j4
cd $__icdl_buildall_py3_workdir__

#echo "<<< BUILD GEO2D >>>"
#cd Geo2D
#source config/setup.sh
#make -j4 >> ${build_log} 2>&1
#make -j4
#cd $__icdl_buildall_py3_workdir__

#echo "<<< BUILD LAROPENCV >>>"
#cd LArOpenCV
#make -j4 >> ${build_log} 2>&1
#make -j4
#cd $__icdl_buildall_py3_workdir__

echo "<<< BUILD LARCV >>>"
cd larcv
mkdir -p build
cd build
cmake -DUSE_PYTHON3=ON -DUSE_OPENCV=OFF -DUSE_FNAL=ON -DUSE_TORCH=OFF ../
#make install -j4 >> ${build_log} 2>&1
make install -j4
cd $__icdl_buildall_py3_workdir__

echo "<<< BUILD CILANTRO >>>"
cd cilantro
mkdir -p build
cd build
cmake ../
#make >> $build_log 2>&1
make
cd $__icdl_buildall_py3_workdir__

echo "<<< BUILD UBLARCVAPP >>>"
mkdir -p ublarcvapp/build
cd ublarcvapp
source configure.sh
cd build
cmake -DUSE_OPENCV=OFF ../
#make install -j4 >> ${build_log} 2>&1
make install -j4 
cd $__icdl_buildall_py3_workdir__

echo "<<< BUILD LARFLOW >>>"
mkdir -p larflow/build
cd larflow
source configure.sh
cd build
cmake -DUSE_PYTHON3=ON ../
#make install -j4 >> ${build_log} 2>&1
make install -j4
cd $__icdl_buildall_py3_workdir__

echo "built icdl modules"
