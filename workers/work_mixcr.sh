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

mkdir -p ${OUT["dir"]}
cd ${OUT["dir"]}

TRA=$(tempfile -d .)
TRB=$(tempfile -d .)
TRG=$(tempfile -d .)
TRX=$(tempfile -d .)
TRXN=$(tempfile -d .)
TRX_1=$(tempfile -d .)
TRX_2=$(tempfile -d .)

IGH=$(tempfile -d .)
IGL=$(tempfile -d .)
IGK=$(tempfile -d .)
IGX=$(tempfile -d .)
IGXN=$(tempfile -d .)
IGX_1=$(tempfile -d .)
IGX_2=$(tempfile -d .)
SEQ_1=$(tempfile -d . -s fq.gz)
SEQ_2=$(tempfile -d . -s fq.gz)

## extract mapped reads
sambamba slice ${IN["alig_csort"]} ${PARAMS["TRA"]} > $TRA
sambamba slice ${IN["alig_csort"]} ${PARAMS["TRB"]} > $TRB
sambamba slice ${IN["alig_csort"]} ${PARAMS["TRG"]} > $TRG
sambamba merge -t ${PARAMS["p"]} $TRX $TRA $TRB $TRG
samtools sort -n -@ ${PARAMS["p"]} $TRX > $TRXN
bedtools bamtofastq -i $TRXN -fq $TRX_1 -fq2 $TRX_2 &> /dev/null 

sambamba slice ${IN["alig_csort"]} ${PARAMS["IGH"]} > $IGH
sambamba slice ${IN["alig_csort"]} ${PARAMS["IGL"]} > $IGL
sambamba slice ${IN["alig_csort"]} ${PARAMS["IGK"]} > $IGK
sambamba merge -t ${PARAMS["p"]} $IGX $IGH $IGL $IGK
samtools sort -n -@ ${PARAMS["p"]} $IGX > $IGXN
bedtools bamtofastq -i $IGXN -fq $IGX_1 -fq2 $IGX_2 &> /dev/null 

## compress
pigz -p ${PARAMS["p"]} $TRX_1
pigz -p ${PARAMS["p"]} $TRX_2
pigz -p ${PARAMS["p"]} $IGX_1
pigz -p ${PARAMS["p"]} $IGX_2

cat $TRX_1.gz $IGX_1.gz ${IN["unmap_1"]} > $SEQ_1
cat $TRX_2.gz $IGX_2.gz ${IN["unmap_2"]} > $SEQ_2

## cleanup
rm $TRA $TRB $TRG $TRX $TRXN $IGH $IGL $IGK $IGX $IGXN $TRX_1.gz $TRX_2.gz $IGX_1.gz $IGX_2.gz

## MiXCR
TEMP_3=$(tempfile -d .)
TEMP_3B=$(tempfile -d .)
TEMP_3T=$(tempfile -d .)
TEMP_3F=$(tempfile -d .)
TEMP_4=$(tempfile -d .)
TEMP_5=$(tempfile -d .)

mixcr align -f -t ${PARAMS["p"]} -g -s hsa -p rna-seq -OallowPartialAlignments=true -r ${OUT["aln.rep"]} $SEQ_1 $SEQ_2 $TEMP_3
mixcr filterAlignments -f -n 250000 -c IG  $TEMP_3 $TEMP_3B
mixcr filterAlignments -f -n 250000 -c TCR $TEMP_3 $TEMP_3T
rm $TEMP_3F
mixcr mergeAlignments $TEMP_3B $TEMP_3T $TEMP_3F
mixcr assemblePartial -f -r ${OUT["fix1.rep"]} $TEMP_3F $TEMP_4
mixcr assemblePartial -f -r ${OUT["fix2.rep"]} $TEMP_4  $TEMP_5
mixcr extendAlignments -f -r ${OUT["ext.rep"]} $TEMP_5 ${OUT["alig"]}
mixcr assemble -f -t ${PARAMS["p"]} -i ${OUT["index"]} -r ${OUT["asm.rep"]} ${OUT["alig"]} ${OUT["clone"]}

rm $SEQ_1 $SEQ_2 $TEMP_3 $TEMP_3B $TEMP_3T $TEMP_3F $TEMP_4 $TEMP_5

date >> ${OUT["main"]}
