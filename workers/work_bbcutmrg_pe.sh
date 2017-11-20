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

mkdir -p ${OUT["dir"]}
cd ${OUT["dir"]}

bbduk2.sh in=${IN["fastq1"]} in2=${IN["fastq2"]} fref=${PARAMS["rrna"]} \
          stats=${OUT["stats"]} out=${OUT["cutfq1"]} out2=${OUT["cutfq2"]} \
          minlength=${PARAMS["min_len"]} ${PARAMS["cutxargs"]} \
          threads=${PARAMS["p"]} -Xmx${PARAMS["xmx"]}g &> ${OUT["cut"]}

if [ ${PARAMS["prune"]} == "True" ]
then
    rm -f $(readlink -e -f ${IN["fastq1"]})
    rm -f $(readlink -e -f ${IN["fastq2"]})
fi

bbmerge.sh in1=${OUT["cutfq2"]} in2=${OUT["cutfq1"]} out=${OUT["mrgfq3"]} \
                outu1=${OUT["mrgfq2"]} outu2=${OUT["mrgfq1"]} ihist=${OUT["isize"]} \
                ${PARAMS["mrgxargs"]} \
                threads=${PARAMS["p"]} usejni=t -Xmx${PARAMS["xmx"]}g &> ${OUT["mrg"]}

## clean fq.gz links
rm ${IN["fastq2"]} ${IN["fastq1"]}

date >> ${OUT["main"]}
