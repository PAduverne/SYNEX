# Makefile for properly lauching the SMBBH simulation on the DANTE cluster

make_population: population_generation.py
	python3 population_generation.py \
    	--N_BBH 500 \
    	--path ./test

setup_test: setup_simulations.py
	python3 setup_simulations.py \
    	--n-sim 5 \
    	--injection-file ./parameters_cbc.dat \
    	--test
run_test:
	./launch_simulation.sh 0 2
