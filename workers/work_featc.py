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

CMD = "samtools view %s | head -n 2 | cut -f 1 | uniq | wc -l" % IN["alig_nsort"]
run = run_cmd(CMD)
log_run(run)
PAIRED_END = (run[1].strip() == "1")

OPTS = ["-T", str(PARAMS["n"]), PARAMS["xargs"]]
if PAIRED_END:
    OPTS.insert(2, "-p")
if PARAMS["stranded"]:
    OPTS.insert(2, "-s %s" % PARAMS["stranded"])
if PARAMS["duplicates"]:
    OPTS.insert(2, PARAMS["duplicates"])
opts = " ".join(OPTS)

TEMPLATE = "featureCounts -a %s -o %s " + opts + " %s"

CMD = TEMPLATE % (PARAMS["gtf"], OUT["cts"], IN["alig_nsort"],)
run = run_cmd(CMD)
log_run(run)

log_run(run_cmd("pigz -p 4 %s" % (OUT["cts"],)))

with open("%s.report" % (OUT["cts"],), "wb") as fh:
    fh.write(run[2]) # stderr

###
sys.exit()
