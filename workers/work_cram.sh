#!/usr/bin/env bash
eval $(python2 - <<EOF
import sys, json, itertools

ARGS = json.load(open("$1"))
PARAMS = ARGS.pop("params")
IN = ARGS.pop("in")
OUT = ARGS.pop("out")

T_PARAMS = "declare -A PARAMS=( " + " ".join(['[%s]="%s"'] * len(PARAMS)) + " )"
F_PARAMS =  T_PARAMS % tuple(itertools.chain(*zip(PARAMS.keys(), PARAMS.values())))
T_IN = "declare -A IN=( " + " ".join(['[%s]="%s"'] * len(IN)) + " )"
F_IN = T_IN % tuple(itertools.chain(*zip(IN.keys(), IN.values())))
T_OUT = "declare -A OUT=( " + " ".join(['[%s]="%s"'] * len(OUT)) + " )"
F_OUT = T_OUT % tuple(itertools.chain(*zip(OUT.keys(), OUT.values())))

print "; ".join((F_PARAMS, F_IN, F_OUT))
EOF
)
set -euo pipefail
IFS=$'\n\t'

date > ${OUT["main"]}

bam squeeze --binMid --binQualS ${PARAMS["qbin"]} --keepDups --in ${IN["alig_csort"]} --out -.ubam | \
    samtools view -@4 -T ${PARAMS["genome"]} -C -o ${OUT["alig_csort"]} -
samtools index ${OUT["alig_csort"]}

bam squeeze --binMid --binQualS ${PARAMS["qbin"]} --keepDups --in ${IN["chim_se_csort"]} --out -.ubam | \
    samtools view -@4 -T ${PARAMS["genome"]} -C -o ${OUT["chim_se_csort"]} -
samtools index ${OUT["chim_se_csort"]}

bam squeeze --binMid --binQualS ${PARAMS["qbin"]} --keepDups --in ${IN["chim_pe_csort"]} --out -.ubam | \
    samtools view -@4 -T ${PARAMS["genome"]} -C -o ${OUT["chim_pe_csort"]} -
samtools index ${OUT["chim_pe_csort"]}

date >> ${OUT["main"]}
