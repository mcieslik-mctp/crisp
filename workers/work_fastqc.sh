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

mkdir -p ${OUT["report"]}
cd ${OUT["report"]}
if [[ -n ${IN["fastq1"]} ]]
then
    ln -sfn ${IN["fastq1"]} .
fi
if [[ -n ${IN["fastq2"]} ]]
then
    ln -sfn ${IN["fastq2"]} .
fi
fastqc *.fq
rm *.fq

date >> ${OUT["main"]}
