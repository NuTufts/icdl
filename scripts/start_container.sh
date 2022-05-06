#!/bin/bash
cvmfs_config probe
singularity shell -B /cvmfs:/cvmfs ~/working/larbys/icarus/icarus-containers/sl7_base.simg
