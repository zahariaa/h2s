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
matlabver="2017b" #2018b doesn't have parpool for whatever reason
walltime="02:59:00"
memory="2GB" # PER CPU
numnodes=12
parfor=true #false

# set up directories
basedir=~/results/
function=h2s_radii
scriptdir=$basedir$function

if [ ! -d "$scriptdir" ]
then
   mkdir $scriptdir
   mkdir $scriptdir/m
   mkdir $scriptdir/r
   mkdir $scriptdir/e
   mkdir $scriptdir/o
fi

### BUILD MATLAB JOB
jobname="$1$2"

##### write matlab runfile
mcode="addpath(genpath('~/msubs'));addpath(genpath('~/h2s'));"
if $parfor;       then
         # PARFOR VERSION
         mcode="$mcode\n\nparpool('local',$numnodes)"
fi
if [[ "${2,,}" == "samples" ]]; then
#if [ (string lower $2)="samples"]; then
   mcode="$mcode\n\n$function(log2space(1,10,10),200,'$1',true,true);\n"
#if [ "${2,,}"="dimensions" ]; then
else
   mcode="$mcode\n\n$function(200,log2space(1,12,12),'$1',true,true);\n"
fi
mfile="${scriptdir}/m/$jobname.m"
echo -e "$mcode" > $mfile
##### END write matlab code

######### construct SLURM (job array) script
jfile="${scriptdir}/$jobname.s"

pbtxt1="#!/bin/bash\n#"
pbtxt2="#SBATCH --account=nklab\n#SBATCH --job-name=$jobname"
pbtxt3="#SBATCH --nodes=1\n#SBATCH --cpus-per-task=$numnodes"
pbtxt4="#SBATCH --time=$walltime"
pbtxt5="#SBATCH --mem-per-cpu=$memory"
pbtxt8="#SBATCH --output=${scriptdir}/o/${jobname}%a.txt"
pbtxt9="#SBATCH --error=${scriptdir}/e/${jobname}%a.txt"
pbtxt9a="\nmodule load matlab/$matlabver\n"

# MAIN MATLAB EXECUTION PBS line
#Command to execute Matlab code
pbtxt11="matlab -nodisplay < $mfile > ${scriptdir}/r/$jobname.txt"
pbtxt10=""
pbtxt12=""

# WRITE SLURM SCRIPT
echo -e "${pbtxt1}\n${pbtxt2}\n${pbtxt3}\n${pbtxt4}\n${pbtxt5}\n${pbtxt8}\n${pbtxt9}\n${pbtxt9a}\n\n${pbtxt10}\n${pbtxt11}\n${pbtxt12}\n" > ${jfile}

# submit array job
jid=`sbatch $jfile`
echo $jid

