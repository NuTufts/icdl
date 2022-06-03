#!/bin/bash
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode v09_52_01 -q e20:prof
setup cmake v3_22_2
setup opencv v4_2_0b -q e20:p392

export OPENCV_INCDIR=${OPENCV_INC}/opencv4/
export OPENCV_LIBDIR=${OPENCV_LIB}

