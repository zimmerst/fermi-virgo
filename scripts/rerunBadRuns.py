'''
Created on Dec 6, 2013

@author: zimmer
'''
import lsf
import sys,glob,os
from time import sleep

if __name__ == '__main__':
    from optparse import OptionParser
    usage       = "Usage: %prog  [options] subdir logdir"
    description = "python script"
    parser      = OptionParser(usage=usage,description=description)
    parser.add_option("-f","--force",dest="force",action="store_true",default=False)
    parser.add_option("-d","--dry",dest="dry",action="store_true",default=False)
    parser.add_option("-q","--queue",dest="queue",default="bullet-medium")
    parser.add_option("-t","--tail",dest="tail",default=None,type=int,help="if set, look through n lines of file-tail")
    (opts, args)= parser.parse_args()
    # idea, try to figure out which jobs have failed.
    # first pass: just re-run job-by-job
    # input1: script dir
    # input2: log dir
    scriptDir = sys.argv[1]
    logDir = sys.argv[2]
    # each script should have a log file.
    matchedScripts = {}
    for sc in glob.glob(os.path.join(scriptDir,"*.*")):
        script = os.path.basename(sc) 
        logs = os.listdir(logDir)
        lf = os.path.join(os.path.abspath(logDir),script)
        rerun = False
        if script in logs:
            logf = logs[logs.index(script)]
            lf = os.path.join(os.path.abspath(logDir),logf)
        else:
            rerun = True # per def. if the script doesn't have an assoc. logfile, re-run
        matchedScripts[script]={"sub":sc,"log":lf,"rerun":rerun,"label":script.split(".")[0]}
    # done collecting, now check status
    for key in matchedScripts:
        proc = matchedScripts[key]
        if not proc["rerun"]:
            # check if log is okay
            #print proc["log"]
            logF = proc['log']
            if lsf.isComplete(logF,lines=opts.tail):
                #print '*found good run*',proc
                proc["rerun"]=False
            elif lsf.isExited(logF,lines=opts.tail):
                #print '*found bad run*',proc
                proc["rerun"]=True
            # check a few more instances of failed logs
            elif lsf.check_log(logF, "cannot stat", exists=True,lines=opts.tail):
                proc["rerun"]=True
            elif lsf.check_log(logF, "cannot access", exists=True,lines=opts.tail):
                proc["rerun"]=True
            else: proc["rerun"]=True
    # sort those that are bad
    badRuns = [matchedScripts[key] for key in matchedScripts if matchedScripts[key]["rerun"]==True]
    # done collecting all runs
    print "Ready to re-run, found %i bad runs, proceed?"%len(badRuns)
    for run in badRuns:
        lsf.bsub(run["label"], "bash %s"%run["sub"], run["log"], submit=not(opts.dry), sleep='30s',q=opts.queue)
        sleep(3)
