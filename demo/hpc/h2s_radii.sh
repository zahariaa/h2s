#!/bin/sh
# 2018-11-09 AZ Created
# first argument is distribution type
# second argument is experiment type
# e.g.,
# ./h2s_radii.sh Uniform dimensions
#
# To release held jobs:
#
#for (( v = 1448325; v<=1448366; v++)) do qrls $v; done

# parameters
matlabver="2018b"
walltime="02:59:00"
memory="31GB" 
numnodes=12
parfor=false #true #false

# set up directories
basedir="~/results/"
function="h2s_radii"
scriptdir=$basedir$function/

# if [ ! -d "$scriptdir" ]
# then
#          mkdir $scriptdir
#          mkdir $scriptdir/e
#          mkdir $scriptdir/m
#          mkdir $scriptdir/o
#          mkdir $scriptdir/r
# fi

### BUILD MATLAB JOB
jobname="fit$jextra"

##### write matlab runfile
mfile="addpath(genpath('~/msubs'));addpath(genpath('~/h2s'));"
if $parfor;       then
         # PARFOR VERSION
         mfile="$mfile\n\nparpool('local',$numnodes)"
fi
if [[ "${2,,}" == "samples" ]]; then
#if [ (string lower $2)="samples"]; then
   mfile="$mfile\n\n$function(log2space(1,10,10),200,'$1',true,true);\n"
#if [ "${2,,}"="dimensions" ]; then
else
   mfile="$mfile\n\n$function(200,log2space(1,12,12),'$1',true,true);\n"
fi
##### END write matlab code

echo $mfile

