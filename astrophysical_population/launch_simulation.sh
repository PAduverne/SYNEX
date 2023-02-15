#!/bin/bash

# q=(1.1 3 10)
q=(1.1)
n_q=${#q[*]}
echo "il y a $n_q mass ratios"
echo "le premier mass ratio est ${q[0]}"

#m_tot=(4 5 6 7)
m_tot=(4)
n_m_tot=${#m_tot[*]}
echo "il y a $n_m_tot masses totales"
echo "la premiere masse totale est ${m_tot[0]}"

n_sim_start=$1
n_sim_end=$2
n_sim=$n_sim_start-$n_sim_end
cut=259200

path="/Users/duverne/Documents/LISA-Athena/SYNEX/astrophysical_population/"


#SBATCH --job-name=NZ_b20_tanh_res #job name

# 80 MPI tasks in total
# Stallo has 16 or 20 cores/node and therefore we take
# a number that is divisible by both
#SBATCH --ntasks=80
##SBATCH --ntasks-per-node=40
#SBATCH --partition quiet #queue names: quiet, bigmem, debug

# run for five minutes format: d-hh:mm:ss #day, hour, minutes, second
#SBATCH --time=30-00:00:00

# 500MB memory per core : this is a hard limit
#SBATCH --mem-per-cpu=500MB

# turn on all mail notification
#SBATCH --mail-user=duverne@apc.in2p3.fr
#SBATCH --mail-type=ALL
jobname=0
for ((i=0; i<$n_q; i++)); do
    for ((j=0; j<$n_m_tot; j++)); do
        for ((k=$n_sim_start; k<$n_sim_end; k++)); do
            rep="q_${q[$i]}_Mtot_${m_tot[$j]}"
            echo $rep
            full_path=$path$rep
            echo $full_path
            # define and create a unique scratch directory
            # SCRATCH_DIRECTORY=/work/rmignonr/runs/${SLURM_JOBID}
            SCRATCH_DIRECTORY=full_path/$k/logfiles #default path to run the simu
            mkdir -p ${SCRATCH_DIRECTORY}
            #echo ${SCRATCH_DIRECTORY}
            #cd ${SCRATCH_DIRECTORY}
            ((jobname = $jobname + 1 ))
            echo $jobname
            sbatch -J $jobname\
            --export jobname=$jobname \
            python SMBBH_inference.py \
            --path $full_path/$k/inference_param_files \
            --skip-inference \
            --cut $cut
     #sleep 15
        done
    done
done

