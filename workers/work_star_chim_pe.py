#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os, sys, json, logging, hashlib
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

try:
    os.mkdir(OUT["dir"])
except OSError, e:
    if e.errno != 17:
        raise
os.chdir(OUT["dir"])

KEEP_FASTQ = PARAMS.pop("keep_fastq", False)

## reverse fastq files for proper order of junction in Chimeric output
files = " ".join([
    IN["mrgfq2"],
    IN["mrgfq1"],
])

# required
PARAMS["genomeDir"] = PARAMS["genomeDir"]

# optional
PARAMS["genomeLoad"] = PARAMS.pop("genomeLoad", "NoSharedMemory")
PARAMS["runThreadN"] = PARAMS.pop("runThreadN", 8)

# forbidden
fields = os.path.basename(OUT["dir"]).split('-')
fields[-1] = fields[-1].split("_")[0]
PARAMS["outSAMattrRGline"] = "ID:%s SM:%s.%s LB:%s PU:%s.%s PL:ILLUMINA CN:MCTP" % \
                             (hashlib.md5("-".join(fields)).hexdigest()[:6], # 'unique' hash-id
                              fields[0], fields[1], fields[2], fields[3], fields[4])
PARAMS["outStd"] = "SAM"
PARAMS["alignTranscriptsPerReadNmax"] = "100000" # fix segfault on some samples
PARAMS["outSAMheaderHD"] = "@HD\tVN:1.4\tSO:queryname"

## STAR PAIRED
PARAMS["readFilesIn"] = files
CMD = "STAR " + \
      " ".join(["--%s %s" % (param, value) for param, value in PARAMS.iteritems() if value is not None]) + \
      " > /dev/null"
log_run(run_cmd(CMD))
if not KEEP_FASTQ:
    CMD = "rm %s %s" % (IN["mrgfq2"], IN["mrgfq1"])
    log_run(run_cmd(CMD))
CMD = "samtools view -@ 8 -1 -o %s %s" % ("Paired.out.bam", "Chimeric.out.sam")
log_run(run_cmd(CMD))
CMD = "mv Log.final.out Log.paired.out"
log_run(run_cmd(CMD))
CMD = "mv Chimeric.out.junction Chimeric.paired.junction"
log_run(run_cmd(CMD))


## STAR SINGLE
PARAMS["readFilesIn"] = IN["mrgfq3"]
CMD = "STAR " + \
      " ".join(["--%s %s" % (param, value) for param, value in PARAMS.iteritems() if value is not None]) + \
      " > /dev/null"
log_run(run_cmd(CMD))


if not KEEP_FASTQ:
    CMD = "rm %s" % (IN["mrgfq3"],)
    log_run(run_cmd(CMD))
CMD = "samtools view -@ 8 -1 -o %s %s" % ("Merged.out.bam", "Chimeric.out.sam")
log_run(run_cmd(CMD))
CMD = "mv Log.final.out Log.merged.out"
log_run(run_cmd(CMD))
CMD = "mv Chimeric.out.junction Chimeric.merged.junction"
log_run(run_cmd(CMD))
CMD = "rm Chimeric.out.sam"
log_run(run_cmd(CMD))

sys.exit(0)
