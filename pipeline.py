#!/usr/bin/env python2
import sys, os, glob
import numpy as np
import pandas as pd
import logging
import multiprocessing
from collections import defaultdict
from math import ceil, floor
from moke import *
from numap import NuMap
from papy.core import Plumber, Piper, Worker
from papy.util.script import script
from papy.util.func import ipasser, npasser

VERSION = "2.5.1"
CODE_DIR = path(os.getenv("CODE") or "/code")
REFS_DIR = path(os.getenv("REFS") or "/refs")
WORK_DIR = CODE_DIR / "workers"

EVALUATOR = "bash"
PREAMBLE = """
#!/usr/bin/env bash
""" 

MERGES = {
     "read":np.array([False, False, False, False, False, True, True]),
     "lane":np.array([False, False, False, False,  True, True, True]),
     "flow":np.array([False, False, False,  True,  True, True, True]),
     "lib": np.array([False, False,  True,  True,  True, True, True]),
     "rep": np.array([False,  True,  True,  True,  True, True, True]),
    }

RESOURCES = {
    # RESOURCES[number_of_cores][GB_per_core]
    # job1 - number of single-core jobs
    # jobn - number of multi-core jobs
    # n - numbers of cores per-jobn
    # memn - number of GBs of memory per jobn
    8: {
       12:{"job1": 8,"jobn":1,"n": 8,"memn":50},
       16:{"job1": 8,"jobn":1,"n": 8,"memn":75}},
    16: {
        6:{"job1":14,"jobn":1,"n":16,"memn":50},
        8:{"job1":14,"jobn":1,"n":16,"memn":75}},
    32: {
        4:{"job1":24,"jobn":1,"n":32,"memn":75},
        8:{"job1":24,"jobn":2,"n":16,"memn":90}},
    64: {
        4:{"job1":40,"jobn":2,"n":32,"memn":90},
        8:{"job1":40,"jobn":4,"n":16,"memn":100}},
}

def set_log(name):
    global LOG
    LOG = logging.getLogger(name)

DEBUG = False
def set_debug(proc_dir):
    global DEBUG
    DEBUG = True
    global DEBUG_DIR
    DEBUG_DIR = proc_dir / "debug"
    mkdir(DEBUG_DIR)

def plog(msg):
    log(msg, logger=LOG)

def dlog(msg):
    if DEBUG:
        with open(DEBUG_DIR / ("debug-%s.log" % os.getpid()), "a") as fh:
            fh.write(msg + "\n")

def dlog_run(code_stdout_stderr_cmd):
    code, stdout, stderr, cmd = code_stdout_stderr_cmd
    dlog("[%s] %s" % (code, cmd))
    dlog(stdout)
    dlog(stderr)

def merge_sheet(lib_sheet, lib_merge):
    dlog("parsing lib sheet: %s merging: %s" % (lib_sheet, lib_merge))
    tbl = np.array(pd.read_csv(lib_sheet, sep="\t", header=-1), dtype="object")
    tbl[tbl == "nan"] = "NA"
    tbl = np.array(sorted(tbl.tolist()))
    ## single-end
    tblse = np.delete(tbl[tbl[:,1] == "se"], 1,1)
    merges = MERGES[lib_merge]
    ids = tblse[:,np.where(~merges)[0]].tolist()
    se_fns = tblse[:,6].tolist()
    se_ids = ["-".join(id_) for id_ in ids]
    ## paired-end
    tblpe = np.delete(tbl[tbl[:,1] == "pe"], 1,1)
    merges = MERGES[lib_merge]
    ids = tblpe[:,np.where(~merges)[0]].tolist()
    pe_fns = tblpe[:,6].tolist()
    pe_ids = ["-".join(id_) for id_ in ids]
    return ((se_ids, se_fns), (pe_ids, pe_fns))

def input_fq(lib_sheet, orig_dir, proc_dir):
    dlog("preparing fastq input from: %s" % lib_sheet)
    (se_ids, se_fns), (pe_ids, pe_fns) = merge_sheet(lib_sheet, "read")
    if se_ids or se_fns:
        chk_exit(1, "single-end libraries are not supported.")
    pe_fqs = defaultdict(lambda: {"fastq1":None, "fastq2":None})
    for i,(id_,fn) in enumerate(zip(pe_ids, pe_fns)):
        read = "fastq%s" % ((i % 2) + 1)
        pe_fqs[id_][read] = orig_dir / fn
    pe_fqs_boxs = [{"rg":k,"fastq1":v["fastq1"],"fastq2":v["fastq2"]} for k,v in pe_fqs.iteritems()]
    return (pe_fqs_boxs, dict(pe_fqs))

def link_fq(inbox, proc_dir):
    fqs = {"1":inbox[0]["fastq1"], "2":inbox[0]["fastq2"]}
    read_group = inbox[0]["rg"]
    dlog("linking fastqs rg:%s 1:%s 2:%s" % (read_group, fqs["1"], fqs["2"]))
    box = {"fastq1":"", "fastq2":""}
    for read, fq in fqs.iteritems():
        if fq:
            out_fq = "%s_%s.fq.gz" % (proc_dir / read_group, read)
            cmd = "ln -sfn %s %s"  % (fq, out_fq)
            dlog_run(run_app(cmd))
            box["fastq%s" % read] = out_fq
    return box

def sample_fq(inbox, proc_dir, qc_frags):
    fqs = {"1":inbox[0]["fastq1"], "2":inbox[0]["fastq2"]}
    dlog("sampling from input files %s %s" % (fqs["1"], fqs["2"]))
    box = {"fastq1":"", "fastq2":""}
    for read, fq in fqs.iteritems():
        if fq:
            out_fq = "%s/%s_head_%s.fq" % (proc_dir, inbox[0]['rg'], read)
            CMD = "seqtk sample %s %s > %s" % (fq, qc_frags, out_fq)
            dlog_run(run_app(CMD))
            box["fastq%s" % read] = out_fq
    return box

def input_bam(lib_sheet, lib_merge, orig_dir, proc_dir):
    dlog("preparing bam by merging: %s" % lib_merge)
    (se_ids_rg, _), (pe_ids_rg, _) = merge_sheet(lib_sheet, "read")
    (se_ids_sm, _), (pe_ids_sm, _) = merge_sheet(lib_sheet, lib_merge)
    if se_ids_sm or se_ids_rg:
        chk_exit(1, "single-end libraries are not supported.")
    pe_bams = defaultdict(lambda:{"alig_nsort":[], "chim_pe_nsort":[], "chim_se_nsort":[],
                                  "alig_csort":[], "chim_pe_csort":[], "chim_se_csort":[],
                                  "unmap_1":[], "unmap_2":[], "sj":[], "junc_pe":[], "junc_se":[]})
    for idr, ids in zip(pe_ids_rg[::2], pe_ids_sm[::2]):
        pe_bams[ids]["alig_nsort"].append(proc_dir / (idr + "_alig_nsort.bam"))
        pe_bams[ids]["chim_pe_nsort"].append(proc_dir / (idr + "_chim_pe_nsort.bam"))
        pe_bams[ids]["chim_se_nsort"].append(proc_dir / (idr + "_chim_se_nsort.bam"))
        pe_bams[ids]["alig_csort"].append(proc_dir / (idr + "_alig_csort.bam"))
        pe_bams[ids]["chim_pe_csort"].append(proc_dir / (idr + "_chim_pe_csort.bam"))
        pe_bams[ids]["chim_se_csort"].append(proc_dir / (idr + "_chim_se_csort.bam"))
        pe_bams[ids]["unmap_1"].append(proc_dir / (idr + "_unmapped_1.fq.gz"))
        pe_bams[ids]["unmap_2"].append(proc_dir / (idr + "_unmapped_2.fq.gz"))
        pe_bams[ids]["sj"].append(proc_dir / (idr + "_alig.sj"))
        pe_bams[ids]["junc_pe"].append(proc_dir / (idr + "_chim_pe.jnc"))
        pe_bams[ids]["junc_se"].append(proc_dir / (idr + "_chim_se.jnc"))
    pe_bams_boxs = [{"alig_nsort":v["alig_nsort"],"chim_pe_nsort":v["chim_pe_nsort"],"chim_se_nsort":v["chim_se_nsort"],
                     "alig_csort":v["alig_csort"],"chim_pe_csort":v["chim_pe_csort"],"chim_se_csort":v["chim_se_csort"],
                     "unmap_1":v["unmap_1"],"unmap_2":v["unmap_2"],"junc_pe":v["junc_pe"],"junc_se":v["junc_se"],
                     "sj":v["sj"],"sample":k} for \
                    k,v in pe_bams.iteritems()]
    return (pe_bams_boxs, dict(pe_bams))

def move_alig_star(inbox, proc_dir):
    dlog("moving STAR files to output directory")
    inp_dir = path(inbox[0]["dir"])
    for fn, sfx in (("Aligned.out.bam", "_alig_nsort.bam"),
                    ("SJ.out.tab", "_alig.sj"),
                    ("Log.final.out", "_alig.log"),
                    ("Unmapped.out.mate1.gz", "_unmapped_1.fq.gz"),
                    ("Unmapped.out.mate2.gz", "_unmapped_2.fq.gz")):
        inp = inp_dir / fn
        out = proc_dir / (inp_dir.splitpath()[-1].split("_1_cutmrg-cutfq1_alig-dir")[0] + sfx)
        cmd = "mv %s %s"  % (inp, out)
        dlog_run(run_app(cmd))

def move_chim_star(inbox, proc_dir):
    dlog("moving STAR files to output directory")
    inp_dir = path(inbox[0]["dir"])
    for fn, sfx in (("Merged.out.bam", "_chim_se_nsort.bam"),
                    ("Paired.out.bam", "_chim_pe_nsort.bam"),
                    ("Chimeric.paired.junction", "_chim_pe.jnc"),
                    ("Chimeric.merged.junction", "_chim_se.jnc"),
                    ("Log.paired.out", "_chim_pe.log"),
                    ("Log.merged.out", "_chim_se.log")):
        inp = inp_dir / fn
        out = proc_dir / (inp_dir.splitpath()[-1].split("_1_cutmrg-mrgfq1_chim-dir")[0] + sfx)
        cmd = "mv %s %s"  % (inp, out)
        dlog_run(run_app(cmd))

def sambamba_csort(inbox, proc_dir, n, csort_mem):
    for k,v in inbox[0].items():
        if "nsort" in k:
            kk = k.replace("nsort", "csort")
            for p in v:
                pp = p.replace("nsort", "csort")
                td = pp.replace(".bam","-tmp")
                dlog_run(run_app("mkdir -p %s" % td))
                CMD = "sambamba sort --tmpdir %s -m %s -t %s -o %s %s" % (
                    td, csort_mem, n, pp, p)
                dlog_run(run_app(CMD))
                dlog_run(run_app("rm -rf %s" % td))
                
def merge_star(inbox, proc_dir, n):
    dlog("merging STAR files for sample: %s" % (inbox[0]["sample"],))
    out = {"sj":           proc_dir / (inbox[0]["sample"] + "_alig.sj"),
           "junc_pe":      proc_dir / (inbox[0]["sample"] + "_chim_pe.jnc"),
           "junc_se":      proc_dir / (inbox[0]["sample"] + "_chim_se.jnc"),
           "alig_nsort":   proc_dir / (inbox[0]["sample"] + "_alig_nsort.bam"),
           "alig_csort":   proc_dir / (inbox[0]["sample"] + "_alig_csort.bam"),
           "chim_pe_nsort":proc_dir / (inbox[0]["sample"] + "_chim_pe_nsort.bam"),
           "chim_se_nsort":proc_dir / (inbox[0]["sample"] + "_chim_se_nsort.bam"),
           "chim_pe_csort":proc_dir / (inbox[0]["sample"] + "_chim_pe_csort.bam"),
           "chim_se_csort":proc_dir / (inbox[0]["sample"] + "_chim_se_csort.bam"),
           "unmap_1":      proc_dir / (inbox[0]["sample"] + "_unmapped_1.fq.gz"),
           "unmap_2":      proc_dir / (inbox[0]["sample"] + "_unmapped_2.fq.gz"),
    }
    for k in ("alig_nsort", "alig_csort",
              "chim_pe_nsort", "chim_se_nsort", "chim_pe_csort", "chim_se_csort"):
        if (len(inbox[0]["alig_nsort"]) > 1):
            CMD = "sambamba merge -t %s %s %s" % (n, out[k], " ".join(inbox[0][k]))
        else:
            CMD = "ln %s %s" % (inbox[0][k][0], out[k])
        dlog_run(run_app(CMD))
        if "csort" in k:
            dlog_run(run_app("sambamba index %s" % (out[k],)))
    for k in ("sj", "junc_pe", "junc_se", "unmap_1", "unmap_2"):
        if (len(inbox[0][k]) > 1):
            CMD = "cat %s > %s" % (" ".join(inbox[0][k]), out[k])
        else:
            CMD = "ln %s %s" % (inbox[0][k][0], out[k])
        dlog_run(run_app(CMD))
    return out

## PREQC
def preqc(proc_dir, job1, jobn, n, stride):


    kmer_cfg = {
        "id":"kmer",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_kmer.sh"),
        "in":(
            "fastq1",
            "fastq2",
        ),
        "out":(
            ("main", "log"),
            ("kmer1", "txt"),
            ("kmer2", "txt"),
        ),
        "params":{
            "k":6,
        }
    }
    
    fastqc_cfg = {
        "id":"fastqc",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_fastqc.sh"),
        "in":(
            "fastq1",
            "fastq2",
        ),
        "out":(
            ("main", "log"),
            ("report", None),
        ),
        "params":{
        }
    }

    map_job1 = NuMap(worker_num=job1, ordered=False, stride=stride, buffer=10000)
    map_jobn = NuMap(worker_num=jobn, ordered=False, stride=stride, buffer=10000)
    p1 = Piper(Worker(sample_fq, (proc_dir, 1000000,)), parallel=map_job1)
    p2 = Piper(Worker(script, (kmer_cfg,)), parallel=map_job1)
    p3 = Piper(Worker(script, (fastqc_cfg,)), parallel=map_job1)
    p4 = Piper(Worker(npasser), parallel=map_job1)
    # topology
    pipeline = Plumber()
    pipeline.add_pipe((p1, p2, p4))
    pipeline.add_pipe((p1, p3, p4))
    return pipeline

## ALIGN
def align(proc_dir, job1, jobn, memn, stride, n, unstranded, genome_idx, full_model,
          genome_seq, rrna_seq, merge_mode, star_mem, alig_star_params, chim_star_params, 
          prune, keep_fastq):


    cutmrg_script = "work_bbcutmrg_pe.sh" if merge_mode == "bb" else "work_cutmrg_pe.sh"
    cutmrg_cfg = {
        "id":"cutmrg",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, cutmrg_script),
        "in":(
            "fastq1",
            "fastq2",
        ),
        "out":(
            ("main", "log"),
            ("cut", "log"),
            ("mrg", "log"),
            ("cutfq1", "fq"),
            ("cutfq2", "fq"),
            ("mrgfq1", "fq"),
            ("mrgfq2", "fq"),
            ("mrgfq3", "fq"),
            ("isize", "txt"),
            ("stats", "txt"),
            ("dir", None),
        ),
        "params":{
            "min_len":25,
            "cutxargs":"k=31 qskip=3 rieb=t tbo=t tpe=t",
            "mrgxargs":"k=31 prefilter=2 minoverlap=10 extend2=20 iterations=5",
            "seq1":"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "seq2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
            "rrna":rrna_seq,
            "prune":prune,
            "xmx":memn,
            "p":n,
        }        
    }
    
    star_alig_params = {
        "keep_fastq":keep_fastq,
        "genomeDir":genome_idx,
        "runThreadN":n,
        "genomeLoad":star_mem,
        ## spliced alignment
        "outFilterType":"BySJout",
        "outSAMstrandField":"intronMotif" if unstranded else "None",
        "alignSJoverhangMin":8,
        "alignSJDBoverhangMin":3,
        "scoreGenomicLengthLog2scale":0,
        "alignIntronMin":20,
        "alignIntronMax":1000000,
        "alignMatesGapMax":1000000,
    }
    star_alig_params.update(alig_star_params)

    star_chim_params = {
        "keep_fastq":keep_fastq,
        "genomeDir":genome_idx,
        "runThreadN":n,
        "genomeLoad":star_mem,
        "outFilterType":"Normal",
        ## chimeric alignment
        "alignIntronMax":150000,
        "alignMatesGapMax":150000,
        "chimSegmentMin":10,
        "chimJunctionOverhangMin":1,
        "chimScoreSeparation":0,
        "chimScoreJunctionNonGTAG":0,
        "chimScoreDropMax":1000,
        "chimScoreMin":1,
    }
    star_chim_params.update(chim_star_params)
    
    star_alig_cfg = {
        "id":"alig",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"python2",
        "script":"%s/%s" % (WORK_DIR, "work_star_alig_pe.py"),
        "in":(
            "cutfq1",
            "cutfq2",
    ),
        "out":(
            ("main", "log"),
            ("dir", None),
    ),
        "params":star_alig_params
    }
    
    star_chim_cfg = {
        "id":"chim",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"python2",
        "script":"%s/%s" % (WORK_DIR, "work_star_chim_pe.py"),
        "in":(
            "mrgfq1",
            "mrgfq2",
            "mrgfq3",
    ),
        "out":(
            ("main", "log"),
            ("dir", None),
    ),
        "params":star_chim_params
    }
    
    pipeline = Plumber()
    map_job1 = NuMap(worker_num=job1, ordered=False, stride=stride, buffer=10000)
    map_jobn = NuMap(worker_num=jobn, ordered=False, stride=stride, buffer=10000)
    p1 = Piper(Worker(link_fq, (proc_dir,)), parallel=map_job1)
    p2 = Piper(Worker(script, (cutmrg_cfg,)), parallel=(map_jobn if merge_mode=="bb" else map_job1))
    p3a = Piper(Worker(script, (star_alig_cfg,)), parallel=map_jobn)
    p3b = Piper(Worker(script, (star_chim_cfg,)), parallel=map_jobn)
    p4a = Piper(Worker(move_alig_star, (proc_dir,)), parallel=map_job1)
    p4b = Piper(Worker(move_chim_star, (proc_dir,)), parallel=map_job1)
    p5 = Piper(Worker(npasser), parallel=map_jobn)
    pipeline.add_pipe((p1, p2, p3a, p4a, p5))
    pipeline.add_pipe((p1, p2, p3b, p4b, p5))
    return pipeline


## CSORT
def csort(proc_dir, job1, jobn, memn, n, stride):
    csort_mem="%sG" % memn
    pipeline = Plumber()
    map_jobn = NuMap(worker_num=jobn, ordered=False, stride=stride, buffer=10000)
    p1 = Piper(Worker(sambamba_csort, (proc_dir, n, csort_mem)), parallel=map_jobn)
    p2 = Piper(Worker(npasser), parallel=map_jobn)
    pipeline.add_pipe((p1, p2))
    return pipeline

## MERGE
def merge(proc_dir, job1, jobn, n, stride):
    pipeline = Plumber()
    map_jobn = NuMap(worker_num=jobn, ordered=False, stride=stride, buffer=10000)
    p1 = Piper(Worker(merge_star, (proc_dir, n)), parallel=map_jobn)
    p2 = Piper(Worker(npasser), parallel=map_jobn)
    pipeline.add_pipe((p1, p2))
    return pipeline


## QUANT
def quant(proc_dir, job1, jobn, n, stride, unstranded, prot_model, full_model, linc_model):

    ## gene counting
    prot_cfg = {
        "id":"prot",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"python2",
        "script":"%s/work_featc.py" % WORK_DIR,
        "in":(
            "alig_nsort",
        ),
        "out":(
            ("main", "log"),
            ("cts", "cts"),
            ("tmp", None)
        ),
        "params":{
            "paired_end":True,
            "stranded":"0", # always unstranded
            "duplicates":"", # count duplicates
            "gtf":prot_model,
            "n":n,
            "xargs":""
        }
    }
    
    full_cfg = {
        "id":"full",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"python2",
        "script":"%s/work_featc.py" % WORK_DIR,
        "in":(
            "alig_nsort",
        ),
        "out":(
            ("main", "log"),
            ("cts", "cts"),
            ("tmp", None)
        ),
        "params":{
            "paired_end":True,
            "stranded":"0" if unstranded else "1", # first read is on the transcript strand
            "duplicates":"", # count duplicates
            "gtf":full_model,
            "n":n,
            "xargs":"",
        }
    }
    
    linc_cfg = {
        "id":"linc",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"python2",
        "script":"%s/work_featc.py" % WORK_DIR,
        "in":(
            "alig_nsort",
        ),
        "out":(
            ("main", "log"),
            ("cts", "cts"),
            ("tmp", None)
        ),
        "params":{
            "paired_end":True,
            "stranded":"0" if unstranded else "1", # first read is on the transcript strand
            "duplicates":"", # count duplicates
            "gtf":linc_model,
            "n":n,
            "xargs":""
        }
    }

    
    cfgs = [prot_cfg, full_cfg, linc_cfg]
    ## topology
    pipeline = Plumber()
    map_job1 = NuMap(worker_num=job1,     ordered=False, stride=stride, buffer=10000)
    map_job4 = NuMap(worker_num=job1 / 4, ordered=False, stride=stride, buffer=10000)

    p1 = Piper(Worker(ipasser), parallel=map_job1)
    p2 = Piper(Worker(npasser), parallel=map_job1)
    for cfg in cfgs:
        p = Piper(Worker(script, (cfg,)), parallel=map_job4)
        pipeline.add_pipe((p1, p, p2))
    return pipeline


## FINAL
def final(proc_dir, job1, jobn, n, stride, full_model, genome_idx, genome_seq, skip_mixcr, skip_cover):

    ## bamstat
    bamstat_cfg = {
        "id":"bamstat",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_bamstat.sh"),
        "in":(
            "alig_csort",
        ),
        "out":(
            ("main", "log"),
            ("idxstat", "txt"),
            ("flgstat", "txt"),
        ),
        "params":{
        }
    }

    virus_cfg = {
        "id":"virus",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_virus.sh"),
        "in":(
            "alig_csort",
        ),
        "out":(
            ("main", "log"),
            ("call", "txt"),
        ),
        "params":{
        }
    }
    
    gzip_cfg = {
        "id":"gzip",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_gzip.sh"),
        "in":(
            "sj",
            "junc_se",
            "junc_pe",
        ),
        "out":(
            ("main", "log"),
        ),
        "params":{
        }
    }
    
    cram_cfg = {
        "id":"pack",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_cram.sh"),
        "in":(
            "alig_csort",
            "chim_pe_csort",
            "chim_se_csort",
        ),
        "out":(
            ("main", "log"),
            ("alig_csort", "cram"),
            ("chim_pe_csort", "cram"),
            ("chim_se_csort", "cram"),
        ),
        "params":{
            "genome":genome_seq,
            "qbin":"2,10,20,25,30,35,40,42"
        }
    }

    cover_cfg = {
        "id":"cover",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"python2",
        "script":"%s/%s" % (WORK_DIR, "work_coverage.py"),
        "in":(
            "alig_csort",
            "chim_pe_csort",
            "chim_se_csort",
        ),
        "out":(
            ("main", "log"),
            ("alig_csort", "bw"),
            ("chim_pe_csort", "bw"),
            ("chim_se_csort", "bw"),
        ),
        "params":{
            "chrom_length":genome_idx / "chrNameLength.txt"
        }
    }

    mixcr_cfg = {
        "id":"mixcr",
        "evaluator":EVALUATOR,
        "preamble":PREAMBLE,
        "dir":proc_dir,
        "executable":"bash",
        "script":"%s/%s" % (WORK_DIR, "work_mixcr.sh"),
        "in":(
            "alig_csort",
            "unmap_1",
            "unmap_2",
        ),
        "out":(
            ("main", "log"),
            ("aln.rep", "txt"),
            ("fix1.rep", "txt"),
            ("fix2.rep", "txt"),
            ("ext.rep", "txt"),
            ("asm.rep", "txt"),
            ("index", "bin"),
            ("alig", "vdjca"),
            ("clone", "clns"),
            ("dir", None),
        ),
        "params":{
            "p":n,
            "TRA":"chr14:21543538-22556037",
            "TRB":"chr7:142290000-142820000",
            "TRG":"chr7:38237000-38382000",
            "IGK":"chr2:88789991-90313268",
            "IGH":"chr14:105580000-106880000",
            "IGL":"chr22:22020000-22927000",
        }
    }
    
    ## topology
    pipeline = Plumber()
    map_job1 = NuMap(worker_num=job1, ordered=False, stride=stride, buffer=10000)
    map_jobn = NuMap(worker_num=jobn, ordered=False, stride=stride, buffer=10000)
    p1 = Piper(Worker(ipasser), parallel=map_job1)
    p2b = Piper(Worker(script, (cram_cfg,)), parallel=map_job1)
    p2c = Piper(Worker(script, (gzip_cfg,)), parallel=map_job1)
    p2d = Piper(Worker(script, (bamstat_cfg,)), parallel=map_job1)
    p2e = Piper(Worker(script, (virus_cfg,)), parallel=map_job1)
    p3 = Piper(Worker(npasser), parallel=map_job1)
    pipeline.add_pipe((p1, p2b, p3))
    pipeline.add_pipe((p1, p2c, p3))
    pipeline.add_pipe((p1, p2d, p3))
    pipeline.add_pipe((p1, p2e, p3))
    if not skip_mixcr:
        p2f = Piper(Worker(script, (mixcr_cfg,)), parallel=map_jobn)
        pipeline.add_pipe((p1, p2f, p3))
    if not skip_cover:
        p2g = Piper(Worker(script, (cover_cfg,)), parallel=map_job1)
        pipeline.add_pipe((p1, p2g, p3))
    return pipeline

@task
def main(lib_sheet, inp=path("."), out=path("."), prune=False,
         lib_merge="lane", unstranded=False, merge_mode="bb",
         alig_star_params="", chim_star_params="",         
         star_mem="LoadAndKeep", ncores=8, memory=8,
         rrna_seq=REFS_DIR/"genomes/Hsapiens_rRNA.fa",
         genome_idx=REFS_DIR/"indices/star/hg38.rna-motr.v2",
         genome_seq=REFS_DIR/"genomes/hg38.rna.fa",
         prot_model=REFS_DIR/"gtf/motr.v2/motr.v2-prot.gtf",
         full_model=REFS_DIR/"gtf/motr.v2/motr.v2-full.gtf",
         linc_model=REFS_DIR/"gtf/fantom5/fantom5-robust-hg38.gtf",
         keep_tmp=False, keep_fastq=False, keep_align=False, keep_merge=False,
         skip_preqc=False, skip_align=False, skip_csort=False, skip_merge=False,
         skip_quant=False, skip_clean=False, skip_final=False, skip_split=False,
         skip_mixcr=False, skip_cover=False, debug=False
     ):
    """CRISP

    input and output:
      - lib_sheet(``path``) library sheet
      - inp(``path``) input path
      - out(``path``) output path
      - prune(``bool``) remove input FASTQ files (careful!)
    input settings:
      - lib_merge(``str``) library merge rule: none, read, lane, flow, lib, rep
      - unstranded(``bool``) set if the library is unstranded
    additional settings:
      - merge_mode(``str``) hp (fast) or bb (memory)
      - alig_star_params(``str``) custom star alignment parameters
      - chim_star_params(``str``) custom star chimeric parameters
    resources:
      - star_mem(``str``) STAR genomeLoad: LoadForever, LoadAndKeep, NoSharedMemory
      - ncores(``int``) number of available cores
      - memory(``int``) total available memory per core (in GB)
    reference files:
      - rrna_seq(``path``) Sequence file for ribosomal (rRNA) sequences
      - genome_idx(``path``) STAR genome index
      - genome_seq(``path``) FASTA genome sequence
      - prot_model(``path``) GTF file w/ protein gene model (unstranded)
      - full_model(``path``) GTF file w/ full gene model (stranded)
      - linc_model(``path``) GTF file w/ lincRNA model (stranded)
    keep settings:
      - keep_tmp(``bool``) keep temporary files and directories
      - keep_fastq(``bool``) keep processed FASTQ files from alignment
      - keep_align(``bool``) keep unmerged BAM files after alignment
      - keep_merge(``bool``) keep merged BAM files
    workflow settings:
      - skip_preqc(``bool``) skip preqc pipeline
      - skip_align(``bool``) skip align pipeline
      - skip_csort(``bool``) skip csort pipeline
      - skip_merge(``bool``) skip merge pipeline
      - skip_quant(``bool``) skip quant pipeline
      - skip_final(``bool``) skip final pipeline
      - skip_mixcr(``bool``) skip mixcr step
      - skip_cover(``bool``) skip cover step
      - skip_clean(``bool``) skip clean step
      - skip_split(``bool``) skip split step
      - debug(``bool``) debugging mode

    """
    set_log("pipeline")
    plog("@crisp pipeline (begin)")

    ORIG_DIR = inp.abspath()
    PROC_DIR = out.abspath()
    chk_exit(*inp_file(lib_sheet))
    chk_exit(*inp_file(genome_idx / "Genome"))
    chk_exit(*inp_file(rrna_seq))
    chk_exit(*inp_file(genome_seq))
    chk_exit(*inp_file(prot_model))
    chk_exit(*inp_file(full_model))
    chk_exit(*inp_file(linc_model))
    chk_exit(*inp_dir(ORIG_DIR))

    alig_star_params = dict([kv.split(":") for kv in alig_star_params.split(",")]) if \
                       alig_star_params else {}
    chim_star_params = dict([kv.split(":") for kv in chim_star_params.split(",")]) if \
                       chim_star_params else {}
    
    if debug:
        set_debug(PROC_DIR)
    TMP_DIR = PROC_DIR / "tmp"
    mkdir(TMP_DIR)
    run_app("cp %s %s/lib_sheet.tsv" % (lib_sheet, PROC_DIR))
    run_app(("echo CRISP %s > " % VERSION) + ("%s/VERSION" % (PROC_DIR,)))

    ## PREPARE
    fq_iboxs, readgs = input_fq(lib_sheet, orig_dir=ORIG_DIR, proc_dir=PROC_DIR)
    rg_iboxs, samples = input_bam(lib_sheet, lib_merge, orig_dir=ORIG_DIR, proc_dir=PROC_DIR)
    qf_iboxs = [{"alig_nsort":PROC_DIR / rgi["sample"] + "_alig_nsort.bam",
                  "alig_csort":PROC_DIR / rgi["sample"] + "_alig_csort.bam",
                  "chim_se_nsort":PROC_DIR / rgi["sample"] + "_chim_se_nsort.bam",
                  "chim_se_csort":PROC_DIR / rgi["sample"] + "_chim_se_csort.bam",
                  "chim_pe_nsort":PROC_DIR / rgi["sample"] + "_chim_pe_nsort.bam",
                  "chim_pe_csort":PROC_DIR / rgi["sample"] + "_chim_pe_csort.bam",
                  "sj":PROC_DIR / rgi["sample"] + "_alig.sj",
                  "junc_pe":PROC_DIR / rgi["sample"] + "_chim_pe.jnc",
                  "junc_se":PROC_DIR / rgi["sample"] + "_chim_se.jnc",
                  "unmap_1":PROC_DIR / rgi["sample"] + "_unmapped_1.fq.gz",
                  "unmap_2":PROC_DIR / rgi["sample"] + "_unmapped_2.fq.gz",
    } for rgi in rg_iboxs]
    
    with open(PROC_DIR / "lib_merge.tsv", "wb") as fh:
        for sample in samples:
            for bam in samples[sample]["alig_nsort"]:
                rg = bam.basename().split("_alig_nsort.bam")[0]
                fq1 = readgs[rg]["fastq1"]
                fq2 = readgs[rg]["fastq2"]
                fh.write("%s\t%s\t%s\t%s\n" % (sample, rg, fq1, fq2))


    ## resources
    res = RESOURCES[ncores][memory]
    ## PREQC
    if not skip_preqc:
        plog("@preqc (begin) job1=%s,jobn=%s,n=%s" % (res["job1"], res["jobn"], res["n"]))
        pipe_preqc = preqc(proc_dir=PROC_DIR,
                           job1=res["job1"], jobn=res["jobn"], stride=len(fq_iboxs), n=res["n"])
        pipe_preqc.start([fq_iboxs])
        pipe_preqc.run()
        pipe_preqc.wait()
        CMD = "mv %s/*/*_fastqc.html %s" % (PROC_DIR, PROC_DIR)
        dlog_run(run_app(CMD))
        CMD = "mv %s/*/*_fastqc.zip %s" % (PROC_DIR, PROC_DIR)
        dlog_run(run_app(CMD))
        plog("@preqc (end)")
            
    ## ALIGN
    if not skip_align:
        if star_mem == "LoadAndKeep":
            run_app("cd /tmp; STAR --genomeLoad Remove --genomeDir %s" % (genome_idx,))
            dlog_run(run_app("cd /tmp; STAR --genomeLoad LoadAndExit --genomeDir %s" % (genome_idx,)))
        plog("@alig (begin) job1=%s,jobn=%s,n=%s" % (res["job1"], res["jobn"], res["n"]))
        pipe_align = align(proc_dir=PROC_DIR, job1=res["job1"], jobn=res["jobn"],
                           memn=res["memn"], stride=res["job1"], n=res["n"],
                           unstranded=unstranded, genome_idx=genome_idx, 
                           full_model=full_model, genome_seq=genome_seq,
                           rrna_seq=rrna_seq, merge_mode=merge_mode, star_mem=star_mem,
                           alig_star_params=alig_star_params,
                           chim_star_params=chim_star_params,
                           prune=prune, keep_fastq=keep_fastq)
        pipe_align.start([fq_iboxs])
        pipe_align.run()
        pipe_align.wait()
        if star_mem == "LoadAndKeep":
            dlog_run(run_app("cd /tmp; STAR --genomeLoad Remove --genomeDir %s" % (genome_idx,)))
        plog("@alig (end)")
        
    ## CSORT
    if not skip_csort:
        plog("@csort (begin) job1=%s,jobn=%s,n=%s" % (res["job1"], res["jobn"], res["n"]))
        pipe_csort = csort(proc_dir=PROC_DIR,
                           job1=res["job1"], jobn=res["jobn"], stride=len(rg_iboxs), n=res["n"],
                           memn=res["memn"])
        pipe_csort.start([rg_iboxs])
        pipe_csort.run()
        pipe_csort.wait()
        plog("@csort (end)")

    ## MERGE
    if not skip_merge:
        plog("@merge (begin) job1=%s,jobn=%s,n=%s" % (res["job1"], res["jobn"], res["n"]))
        pipe_merge = merge(proc_dir=PROC_DIR,
                           job1=res["job1"], jobn=res["jobn"], stride=len(rg_iboxs), n=res["n"])
        pipe_merge.start([rg_iboxs])
        pipe_merge.run()
        pipe_merge.wait()
        plog("@merge (end)")
    
    ## QUANT    
    if not skip_quant:
        plog("@quant (begin) job1=%s,jobn=%s,n=%s" % (res["job1"], res["jobn"], res["n"]))
        pipe_quant = quant(proc_dir=PROC_DIR,
                           job1=res["job1"], jobn=res["jobn"], stride=len(qf_iboxs), n=res["n"],
                           unstranded=unstranded, prot_model=prot_model, full_model=full_model,
                           linc_model=linc_model)
        pipe_quant.start([qf_iboxs])
        pipe_quant.run()
        pipe_quant.wait()
        plog("@quant (end)")

    ## FINAL
    if not skip_final:
        plog("@final (begin) job1=%s,jobn=%s,n=%s" % (res["job1"], res["jobn"], res["n"]))
        pipe_final = final(proc_dir=PROC_DIR,
                           job1=res["job1"], jobn=res["jobn"], stride=len(qf_iboxs), n=res["n"],
                           full_model=full_model, genome_idx=genome_idx, genome_seq=genome_seq,
                           skip_mixcr=skip_mixcr, skip_cover=skip_cover)
        pipe_final.start([qf_iboxs])
        pipe_final.run()
        pipe_final.wait()
        plog("@final (end)")

    ## CLEAN
    if not skip_clean:
        plog("@clean (begin)")
        ## move directories
        CMD = "mv %s/*cut*dir %s" % (PROC_DIR, TMP_DIR)
        dlog_run(run_app(CMD))
        CMD = "mv %s/*-report %s" % (PROC_DIR, TMP_DIR)
        dlog_run(run_app(CMD))
        CMD = "mv %s/*mixcr-dir %s" % (PROC_DIR, TMP_DIR)
        dlog_run(run_app(CMD))
        if not keep_tmp:
            CMD = "rm %s/*head*.fq" % (PROC_DIR,)
            dlog_run(run_app(CMD))
        if not keep_align:
            for ibox in rg_iboxs:
                CMD = "rm %s" % (
                    " ".join(ibox["sj"] + ibox["junc_pe"] + ibox["junc_se"] + \
                             ibox["chim_pe_nsort"] + ibox["chim_se_nsort"] + ibox["alig_nsort"] + \
                             ibox["chim_pe_csort"] + ibox["chim_se_csort"] + ibox["alig_csort"] + \
                             [bam + ".bai" for bam in ibox["chim_pe_csort"]] + \
                             [bam + ".bai" for bam in ibox["chim_se_csort"]] + \
                             [bam + ".bai" for bam in ibox["alig_csort"]] + \
                             ibox["unmap_1"] + ibox["unmap_2"]),)
                dlog_run(run_app(CMD))
        if not keep_merge:
            for ibox in qf_iboxs:
                for key in ('alig_nsort', 'chim_pe_nsort', 'chim_se_nsort'):
                    CMD = "rm %s" % (ibox[key],)
                    dlog_run(run_app(CMD))
                for key in ('alig_csort', 'chim_pe_csort', 'chim_se_csort'):
                    CMD = "rm %s %s.bai" % (ibox[key], ibox[key])
                    dlog_run(run_app(CMD))
        plog("@clean (end)")
        
    ## SPLIT
    if not skip_split:
        plog("@split (begin)")
        for rgi in rg_iboxs:
            sample_path = PROC_DIR / rgi["sample"]
            sample_path.mkdir_p()
            dlog_run(run_app("mv %s/%s?* %s" % (PROC_DIR, rgi["sample"], sample_path)))
        plog("@split (end)")
        
    ## FINISH
    plog("@crisp pipeline (end)")
    sys.exit(0)


if __name__ == "__main__":
    task()
