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

path="/Users/duverne/Documents/LISA-Athena/SYNEX/astrophysical_population"


jobname=0
for ((i=0; i<$n_q; i++)); do
    for ((j=0; j<$n_m_tot; j++)); do
        for ((kk=$n_sim_start; kk<$n_sim_end; kk++)); do
            rep="q_${q[$i]}_Mtot_${m_tot[$j]}"
            echo $rep
            fullpath=${path}/${rep}/${kk}/inference_param_files
            echo "Path to the parameter files : $fullpath"
            echo "jobname = $jobname"
            ((jobname = $jobname + 1 ))
            # sbatch -J $jobname \
            # --export jobname=$jobname,fullpath=${fullpath},cut=${cut} \
            # ./run_slurm.sh
            sleep 10
        done
    done
done

