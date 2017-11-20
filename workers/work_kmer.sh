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

jellyfish count --disk --text -o ${OUT['kmer1']}.tmp -s 10M -m ${PARAMS['k']} ${IN['fastq1']}
sed 's/}/\n/g' ${OUT['kmer1']}.tmp | tail -n +1 > ${OUT['kmer1']}
rm ${OUT['kmer1']}.tmp

jellyfish count --disk --text -o ${OUT['kmer2']}.tmp -s 10M -m ${PARAMS['k']} ${IN['fastq2']}
sed 's/}/\n/g' ${OUT['kmer2']}.tmp | tail -n +1 > ${OUT['kmer2']}
rm ${OUT['kmer2']}.tmp

date >> ${OUT["main"]}
