#!/bin/bash

file=temp_file.txt;
rm -f $file
for i in `seq 1 100`
do
    ./one_d_act_inhibit > $file
    grep "Newton Step" $file | tail -1 | cut -d' ' -f3
done
