#!/bin/sh

for seed in $(seq 1 1 1000)
do
sed -i '1d' Tau_Sm1_Sm_$seed.dat
done

