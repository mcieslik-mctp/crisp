#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys, json, logging

from subprocess import Popen, PIPE
ARGS = json.load(open(sys.argv[1]))
PARAMS = ARGS.pop("params")
IN = ARGS.pop("in")
OUT = ARGS.pop("out")

log = OUT.get("main", None)
if log:
    out = open(log, "wb")
else:
    out = sys.stdout

def run_cmd(cmd):
    app = Popen(cmd, stdin=None, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = app.communicate(None)
    code = app.wait()
    return (code, stdout, stderr, cmd)

def log_run(code_stdout_stderr_cmd):
    code, stdout, stderr, cmd = code_stdout_stderr_cmd
    if not code:
        out.write("Success: %s\n" % cmd)
    else:
        out.write("Failure (code: %s): %s\n" % (code, cmd))
        out.write("stdout: \n" + stdout)
        out.write("stderr: \n" + stderr)
    return code

###

import numpy as np


for ii in ("alig_csort", "chim_se_csort", "chim_pe_csort"):
    CMD1 = "samtools idxstats %s" % IN[ii]
    RUN1 = run_cmd(CMD1)
    log_run(RUN1)

    cargo = RUN1[1].split("\n")
    counts = [map(int, l.split("\t")[2:4]) for l in cargo if l.strip()]
    counts = np.array(counts)
    MAPPED, _ = counts.sum(axis=0)
    
    TMP = IN[ii] + ".cov"
    AWK = '{print $1 "\t" $2 "\t" $3 "\t" ($4 * 1.0e6 / %s)}' % MAPPED
    PAR = (IN[ii], PARAMS["chrom_length"], AWK, TMP)

    BAM2COV = "bedtools genomecov -bg -split -ibam %s -g %s | awk '%s' - > %s" % PAR
    COV2BGR = "bedGraphToBigWig %s %s %s" % (TMP, PARAMS["chrom_length"], OUT[ii])
    RM = "rm %s" % TMP

    CMD = " && ".join([BAM2COV, COV2BGR, RM])
    log_run(run_cmd(CMD))
###

sys.exit()
