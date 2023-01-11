#!/bin/bash
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode v09_64_01 -q e20:prof
unsetup larcv2
setup cetmodules v3_20_00
#setup cmake v3_22_2
