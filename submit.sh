#!/bin/bash


function usage() {
    printf "usage: $0 [-h] [-n <name>] [-u] executable walltime ranks gpus\n"
}


function help_mesg() {
    usage
    cat <<EOF
Positional arguments:
executable          Program to run.
walltime            Wall clock time for run in minutes.
ranks               Number of MPI ranks to use.
gpus                Number of GPUs to use.
Optional arguments:
-h,--help           Prints this help message.
-n,--name <name>    Name of the job.
-u,--unique         Append timestamp on runscript
EOF
}


#Parse command line args.
pos_args=("" "" "" "")
c=0
jobname="job"
while [[ $# -gt 0 ]]; do
    arg="$1"
    case $arg in
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
script="run_$executable"
if [ ! -z "$unique" ]; then
    now=`date +%s`
    script="${script}.${now}"
fi
cat > $script << EOF
#BSUB -P GEN127
#BSUB -W 00:${walltime}
#BSUB -nnodes 1
#BSUB -J ${jobname}
#BSUB -o ${jobname}.%J
#BSUB -e ${jobname}.%J

cd /gpfs/wolf/gen127/scratch/$USER
jsrun --nrs 1 --tasks_per_rs $ranks --cpu_per_rs $ranks --gpu_per_rs $gpus $executable
EOF

#Submit the script
bsub $script
