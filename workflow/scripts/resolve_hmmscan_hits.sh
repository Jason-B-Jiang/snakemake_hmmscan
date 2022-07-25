#!/bin/bash

################################################################################

# Goal: Use cath-resolve-hits to resolve overlapping hmmscan hits

# Input: $1 to n - 1: All orthogroup architectures
#        $n: cath-resolve-hits executable

################################################################################

args=( "$@" )
crh=${args[-1]}
results_dir=$(dirname $(dirname ${args[0]}))

# make directory for cath-resolve-hits outputs
mkdir ${results_dir}/resolved_domain_architectures

# iterate up until last argument, which is the crh executable
length=$(expr ${#args[@]} - 1)
echo $length

for (( i=0; i<${length}; i++ ));
do
    orthogroup=$(basename ${args[$i]})

    ${crh} --input-format hmmscan_out --worst-permissible-bitscore 0.1 \
    ${args[$i]} > ${results_dir}/resolved_domain_architectures/${orthogroup}
done