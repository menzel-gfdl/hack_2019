#!/bin/bash


function usage() {
    printf "usage: $0 [-h] [-n <name>] [-u] [-w <directory>] [-d] executable walltime ranks gpus\n"
}


function help_mesg() {
    usage
    cat <<EOF

Positional arguments:
executable                         Program to run.
walltime                           Wall clock time for run in minutes.
ranks                              Number of MPI ranks to use.
gpus                               Number of GPUs to use.

Optional arguments:
-d,--debug                         Attempt to the run the program in DDT.
-h,--help                          Prints this help message.
-n,--name <name>                   Name of the job.  Defaults to "job".
-u,--unique                        Append timestamp on runscript.
-w,--work-directory <directory>    Directory to run the executable in (will be created
                                   if it doesn't exist).  Defaults to 
                                   /gpfs/wolf/gen127/scratch/$USER/princeton_gpu_hackathon.
EOF
}


#Parse command line args.
pos_args=("" "" "" "")
c=0
jobname="job"
workdir="/gpfs/wolf/gen127/scratch/$USER/princeton_gpu_hackathon"
while [[ $# -gt 0 ]]; do
    arg="$1"
    case $arg in
        -d|--debug)
            debug="ddt --connect"
            shift
            ;;
        -h|--help)
            help_mesg
            exit 0
            ;;
        -n|--name)
            jobname="$2"
            if [ -z "$jobname" ]; then
                usage
                exit 1
            fi
            shift
            shift
            ;;
        -u|--unique)
            unique="true"
            shift
            ;;
        -w|--work-directory)
            workdir="$2"
            if [ -z "$workdir" ]; then
                usage
                exit 1
            fi
            shift
            shift
            ;;
        *)
            pos_args[$c]=$1
            c=$((c+1))
            shift
            ;;
    esac
done
if [ $c -ne 4 ]; then
    usage
    exit 1
fi
executable=${pos_args[0]}
walltime=${pos_args[1]}
ranks=${pos_args[2]}
gpus=${pos_args[3]}

#Create submission script.
execname=`basename $executable`
script="run_${execname}"
if [ ! -z "$unique" ]; then
    now=`date +%s`
    script="${script}.${now}"
fi
if [ -z "$debug" ]; then
    cmd="jsrun --nrs 1 --tasks_per_rs $ranks --cpu_per_rs $ranks --gpu_per_rs $gpus ./$execname"
else
    cmd="source $MODULESHOME/init/bash && module load forge/19.0.2 && $debug jsrun --nrs 1 --tasks_per_rs $ranks --cpu_per_rs $ranks --gpu_per_rs $gpus ./$execname"
fi
cat > $script << EOF
#!/bin/bash -xe
#BSUB -P GEN127
#BSUB -W 00:${walltime}
#BSUB -nnodes 1
#BSUB -J ${jobname}
#BSUB -o ${jobname}.%J
#BSUB -e ${jobname}.%J

mkdir -p $workdir
cp -f $executable $workdir
cd $workdir
$cmd
EOF

#Submit the script
bsub $script
