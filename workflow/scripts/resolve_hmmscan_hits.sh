#!/bin/bash

################################################################################

# Goal: Use cath-resolve-hits to resolve overlapping hmmscan hits

# Input: $1 to n - 1: All orthogroup architectures
#        $n: cath-resolve-hits executable

################################################################################

args=( "$@" )
hmmscan_output_dir=${args[0]}
crh=${args[1]}
results_dir=$(dirname $args)

# # make directory for cath-resolve-hits outputs
mkdir ${results_dir}/resolved_domain_architectures

for hmmscan_file in $hmmscan_output_dir/*;
do
    orthogroup=$(basename ${hmmscan_file})
    ${crh} --input-format hmmscan_out --worst-permissible-bitscore 0.1 \
    ${hmmscan_file} > ${results_dir}/resolved_domain_architectures/${orthogroup}
done