#!/bin/bash

# WE NEED TO SETUP ENV VARIABLES FOR ROOT, CUDA, OPENCV
alias python=python3
MACHINE=`uname --nodename`

if [ $MACHINE == "trex" ]
then
    echo "SETUP TREX"
    source /usr/local/root/6.16.00_py3/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-10.0
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/opencv/opencv-3.4.6/include
    export OPENCV_LIBDIR=/usr/local/opencv/opencv-3.4.6/lib

elif [ $MACHINE == "meitner" ]
then
    echo "SETUP MEITNER"

    source /usr/local/root/v6.24.06/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-10.0
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/include/opencv4
    export OPENCV_LIBDIR=/usr/local/lib    

elif [ $MACHINE == "goeppert" ]
then
    echo "SETUP GOEPPERT"
    source /usr/local/root/6.16.00_py3/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/include
    export OPENCV_LIBDIR=/usr/local/lib    

elif [ $MACHINE == "mayer" ]
then
    echo "SETUP MAYER"
    source /usr/local/root/6.16.00_python3/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-10.0
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/local/include
    export OPENCV_LIBDIR=/usr/local/lib    

elif [ $MACHINE == "blade" ]
then
    echo "SETUP TARITREE's RAZER BLADE"
    source /home/twongjirad/software/root6/py3_build/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include
    export OPENCV_LIBDIR=/usr/local/lib    
    
elif [ $MACHINE == "pop-os" ]
then
    echo "SETUP TARITREE's RAZER BLADE PRO"
    source /home/twongjirad/software/root/build_gcc11/bin/thisroot.sh

    export CUDA_HOME=/usr/lib/cuda/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include/opencv4
    export OPENCV_LIBDIR=/usr/lib/x86_64-linux-gnu/
    
elif [ $MACHINE == "mmr-Alienware-x15-R1" ]
then
    echo "SETUP MATT's ALIENWARE x15 R1"
    source /home/matthew/software/root/install/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda-11.3/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include/opencv4
    export OPENCV_LIBDIR=/usr/lib/x86_64-linux-gnu/

    export LIBTORCH_DIR="/home/matthew/software/pythonVEnvs/myPy3.8.10/lib/python3.8/site-packages/torch"
    export LIBTORCH_LIBDIR=${LIBTORCH_DIR}/lib
    export LIBTORCH_INCDIR=${LIBTORCH_DIR}/include
    [[ ":$LD_LIBRARY_PATH:" != *":${LIBTORCH_LIBDIR}:"* ]] && \
        export LD_LIBRARY_PATH="${LIBTORCH_LIBDIR}:${LD_LIBRARY_PATH}"
   
else
    echo "DEFAULT SETUP (COMPAT WITH SINGULARITY CONTAINER)"
    source /usr/local/root/bin/thisroot.sh
    #source /usr/local/root/root-6.22.06/bin/thisroot.sh

    export CUDA_HOME=/usr/local/cuda/
    [[ ":$LD_LIBRARY_PATH:" != *":${CUDA_HOME}/lib64:"* ]] && export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

    export OPENCV_INCDIR=/usr/include
    export OPENCV_LIBDIR=/usr/local/lib
    
fi

# LIBTORCH
# location below is typically where running `pip install torch` will put pytorch
# export LIBTORCH_DIR="/usr/local/lib/python2.7/dist-packages/torch"
#export LIBTORCH_DIR="/home/jmills/.local/lib/python3.5/site-packages/torch"
#export LIBTORCH_LIBDIR=${LIBTORCH_DIR}/lib
#export LIBTORCH_INCDIR=${LIBTORCH_DIR}/lib/include
#export LIBTORCH_DIR="/usr/local/torchlib"
#export LIBTORCH_LIBDIR=${LIBTORCH_DIR}/lib
#export LIBTORCH_INCDIR=${LIBTORCH_DIR}/include
#[[ ":$LD_LIBRARY_PATH:" != *":${LIBTORCH_LIBDIR}:"* ]] && \
#    export LD_LIBRARY_PATH="${LIBTORCH_LIBDIR}:${LD_LIBRARY_PATH}"
