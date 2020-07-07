#!/bin/bash
# 2020-07-07 AZ Created from allcells.sh
# first argument is name of function. it will be run on instances from 1 to second argument
# third argument (and beyond) are passed to function
# e.g.,
# ./hpcarray.sh simulate_stats 75 1
# ./hpcarray.sh simulate_stats 50-75 2
#
# To release held jobs:
#
#for (( v = 1448325; v<=1448366; v++)) do qrls $v; done

# parameters
matlabver="2017b" #2018b doesn't have parpool for whatever reason
walltime="00:19:00"
memory="1500MB" # PER CPU
numnodes=12
parfor=true #false

# set up directories
function=$1
firstcell=1
if [ "$#" -eq 1 ]; then
   ncells=1
fi
basedir=~/results/
scriptdir=$basedir$function/

if [ ! -d "$scriptdir" ]
then
   mkdir $scriptdir
   mkdir $scriptdir/m
   mkdir $scriptdir/r
   mkdir $scriptdir/e
   mkdir $scriptdir/o
fi

# distribute extra input arguments to matlab function
extraargs=""
for var in "$@"
do
   case $var in
      $1) # do nothing
         ;;
      $2)
         if [ $(echo $var | grep -e "-") ]; then
            firstcell=${var%-*}
            ncells=${var##*-}
         else
            ncells=$var
         fi
         ;;
      *)
         extraargs="$extraargs,$var"
         ;;
   esac
   #echo $extraargs
done
#echo $firstcell-$ncells

######### LOOP ALL CELLS
jobname=$function$3
mfilebase="${scriptdir}m/$jobname"
for ((  icell = $firstcell ;  icell <= $ncells;  icell++  ))
do
   mfile="$mfilebase$icell.m"
   
   ##### write matlab runfile code
   mcode="addpath(genpath('~/msubs'));addpath(genpath('~/h2s'));"
   if $parfor;       then
      # PARFOR VERSION
      mcode="$mcode\n\nparpool('local',$numnodes)"
   fi
   mcode="$mcode\n\n$function"
   if [ "$#" -ne 1 ]; then
      mcode="$mcode($icell$extraargs)"
   fi
   ##### END write matlab runfile code
   echo -e "$mcode" > $mfile
done

######### construct SLURM (job array) script
jfile="${scriptdir}/$jobname.s"

pbtxt1="#!/bin/bash\n#"
pbtxt2="#SBATCH --account=nklab\n#SBATCH --job-name=$jobname"
pbtxt3="#SBATCH --nodes=1\n#SBATCH --cpus-per-task=$numnodes"
pbtxt4="#SBATCH --time=$walltime"
pbtxt5="#SBATCH --mem-per-cpu=$memory"
pbtxt8="#SBATCH --output=${scriptdir}o/${jobname}%a.txt"
pbtxt9="#SBATCH --error=${scriptdir}e/${jobname}%a.txt"
pbtxt9a="\nmodule load matlab/$matlabver\n"

# MAIN MATLAB EXECUTION PBS line
#Command to execute Matlab code
pbtxt11="matlab -nodisplay < $mfilebase\$SLURM_ARRAY_TASK_ID.m > ${scriptdir}r/$jobname\$SLURM_ARRAY_TASK_ID.txt"
pbtxt10=""
pbtxt12=""

# WRITE SLURM SCRIPT
echo -e "${pbtxt1}\n${pbtxt2}\n${pbtxt3}\n${pbtxt4}\n${pbtxt5}\n${pbtxt8}\n${pbtxt9}\n${pbtxt9a}\n\n${pbtxt10}\n${pbtxt11}\n${pbtxt12}\n" > ${jfile}

# submit array job
jid=`sbatch --array=$firstcell-$ncells $jfile`
echo $jid

