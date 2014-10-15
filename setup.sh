#!/bin/bash

ulimit -S -c 0 # no core-dumps
source ~/scripts/set_STools.sh 09-31-01
#source ~/scripts/set_STools.sh 09-28-00
export PYTHONPATH=/u/gl/zimmer/storage/lib/python:$PYTHONPATH
export PYTHONPATH=$PWD/fermi-virgo:$PYTHONPATH
export TEMPDIR=/afs/slac/g/glast/users/zimmer/sandbox/trash/
export LOGDIR=/afs/slac/g/glast/users/zimmer/logs/
export MYTMPDIR=$TEMPDIR

mkdir -p $TMPDIR $LOGDIR

export GAL=/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/ring_2year_P76_v0.fits
export EGAL=/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/isotrop_2year_P76_source_v0.txt
export PYTHONPATH=/nfs/farm/g/glast/u55/zimmer/FermiData/ClusterLines/dev/pointlike/python:$PYTHONPATH
