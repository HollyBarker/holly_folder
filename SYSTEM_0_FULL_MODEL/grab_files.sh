#!/bin/bash
# First line expands aliases (e.g. PATH, EDITOR, etc.). Second line loads in
# aliases/functions from the .bashrc file
shopt -s expand_aliases 
source ~/.bashrc

diffusivity=(2 4 6 8 10);
steps=(0);
for stepno in ${steps[@]}
do
    cp RESLT_DIFFUSION_0._05/step$stepno.dat STRESS_TIMESTEPS/diff_0._05"step"$stepno.dat
    for diff in ${diffusivity[@]}
    do
	cp RESLT_DIFFUSION_0._$diff/step$stepno.dat STRESS_TIMESTEPS/diff_0._$diff"step"$stepno.dat
    done
done
