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

## cut
hpcut -l ${PARAMS["min_len"]} -3 ${PARAMS["seq1"]} ${IN["fastq1"]} > ${OUT["cutfq1"]} 2>> ${OUT["cut"]}
if [ ${PARAMS["prune"]} == "True" ]
then
    rm -f $(readlink -e -f ${IN["fastq1"]})
fi
hpcut -l ${PARAMS["min_len"]} -3 ${PARAMS["seq2"]} ${IN["fastq2"]} > ${OUT["cutfq2"]}  2>> ${OUT["cut"]}
if [ ${PARAMS["prune"]} == "True" ]
then
    rm -f $(readlink -e -f ${IN["fastq2"]})
fi

## merge
STUMP=${IN["fastq1"]}
hpmerge -o $STUMP ${OUT["cutfq2"]} ${OUT["cutfq1"]} > ${OUT["mrg"]}
mv ${STUMP}_2.fq ${OUT["mrgfq1"]}
mv ${STUMP}_1.fq ${OUT["mrgfq2"]}
mv ${STUMP}_0.fq ${OUT["mrgfq3"]}

## clean fq.gz links
rm ${IN["fastq2"]} ${IN["fastq1"]}
echo "# not supported" > ${OUT["isize"]}

date >> ${OUT["main"]}
