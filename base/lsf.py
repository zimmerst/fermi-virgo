#!/usr/bin/env python
import subprocess as sub
import os, shlex, time
from os.path import expandvars, basename, dirname, join
from numpy import linspace

def mkdir(dir):
    if not os.path.exists(dir):  os.makedirs(dir)
    return dir

def parse_sleep(sleep):
    MINUTE=60
    HOUR=60*MINUTE
    DAY=24*HOUR
    WEEK=7*DAY
    if isinstance(sleep,float) or isinstance(sleep,int):
        return sleep
    elif isinstance(sleep,str):
        try: return float(sleep)
        except ValueError: pass

        if sleep.endswith('s'): return float(sleep.strip('s'))
        elif sleep.endswith('m'): return float(sleep.strip('m'))*MINUTE
        elif sleep.endswith('h'): return float(sleep.strip('h'))*HOUR
        elif sleep.endswith('d'): return float(sleep.strip('d'))*DAY
        elif sleep.endswith('w'): return float(sleep.strip('w'))*WEEK
        else: raise ValueError
    else:
        raise ValueError


def bsub(jobname, command, logfile=None, submit=True, sleep='1m', **kwargs):
    # Just one command...
    if 'q' in kwargs and kwargs['q']  == 'local':
        if isinstance(command,str):
            job = command
        else:
            raise Exception("Cannot run job array locally.")

    else:
        if isinstance(command,str):
            job = create_job(jobname, command, logfile)
        else:
            job = create_job_array(jobname,command, logfile, sleep)
        # allow 'time submission' -> handles startup of jobs more easily
        if 'q' in kwargs:
            qval = kwargs['q']
            time_submit = False
            if ":" in qval:
                time_submit = True            
            try:
                qval = int(kwargs['q'])
                time_submit = True
            except ValueError:
                pass # do nothing
            # replace -q with -W if time is given! 
            if time_submit:
                kwargs['W']=qval
                kwargs.pop('q') # remove queue
        # END FIX
        opts = parse_opts(**kwargs)
        job = "bsub "+ opts + job

    print job
    if submit: status = sub.call(shlex.split(job))
    else:      status = 0
    print
    return status

def parse_opts(**kwargs):
    if 'q' not in kwargs: kwargs['q'] = 'long'
    if 'bullet' in kwargs['q']: 
        kwargs['q'] = kwargs['q'].split('-')[-1]
        kwargs['R'] = 'rhel60' if 'R' not in kwargs else kwargs['R']+' && rhel60'

    for k,v in kwargs.items(): 
        if v is None: kwargs.pop(k)
    return ''.join('-%s "%s" '%(k,v) for k,v in kwargs.items())
    

def create_job(jobname, command, logfile):
    params = dict(name=jobname, 
                  cmnd=command, 
                  log=logfile)
    if logfile is None:
        job = """-J %(name)s %(cmnd)s"""%(params)
    else:
        job = """-oo %(log)s -J %(name)s %(cmnd)s"""%(params)
    return job

def create_job_array(jobname, commands, logfiles=None, sleep='1m'):
    subdir=mkdir("sub")
    outdir=mkdir("log")

    subbase=os.path.join(subdir,basename(jobname))
    outbase=os.path.join(outdir,basename(jobname))

    create_scripts(commands,subbase,sleep)
    if logfiles is not None:
        link_logfiles(logfiles,outbase)

    submit="sh "+subbase+".${LSB_JOBINDEX}"; 
    output=outbase+".%I"

    njobs = len(commands)
    params = dict(name=jobname,
                  cmnd=submit,
                  log=output,
                  njobs=njobs)
                  
    job = """-oo %(log)s -J %(name)s[1-%(njobs)i] %(cmnd)s"""%(params)
    return job

def create_scripts(commands, subbase="submit", sleep='1m'):
    # Some interesting things we do here:
    # 1) cat the script for the logfile
    # 2) sleep to prevent overload
    # 3) Return the exit value of the command
    sleeps = linspace(0,parse_sleep(sleep),len(commands))
    for i, (command,sleep) in enumerate(zip(commands,sleeps)):
        filename = subbase+".%i"%(i+1)
        f = open(filename,'w')
        f.write('%s'%(os.path.basename(filename)).center(35,'#'))
        f.write('\n\n')
        f.write("cat $0;\n")
        f.write("sleep %i;\n"%(sleep))
        f.write(command)
        f.write("\nexit $?;")
        f.write("\n\n")
        f.write("Output follows...".center(35,'#'))
        f.write("\n\n")
        f.close()

def link_logfiles(logfiles,outbase="output"):
    for i,log in enumerate(logfiles):
        log = os.path.expandvars(log)
        output = "%s.%i"%(outbase,i+1)
        if os.path.lexists(log):    
            if log == os.devnull: pass
            else: os.remove(log)
        if os.path.lexists(output): os.remove(output)
        os.symlink(log,output)
    

# For checking logfiles...
def check_log(logfile, string='Successfully', exists=True):
    """ Often logfile doesn't exist because the job hasn't begun
    to run. It is unclear what you want to do in that case...
    logfile : String with path to logfile
    exists  : Is the logfile required to exist
    string  : Value to check for in existing logfile
    """
    if not os.path.exists(logfile):
	return not exists
    return string in open(logfile).read()

def isComplete(logfile):
    return check_log(logfile, "Successfully complete", True)

def isExited(logfile):
    return check_log(logfile, "Exited", True)

def random_sleep(sleep):
    RANDOM=32767
    # convert to seconds
    sleep = parse_sleep(sleep)
    return int(RANDOM/sleep)

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()
    
    commands = [
        "echo $HOME",
        "ls $HOME",
        "ls -lh $HOME"
        ]

    logfiles = [
        "${HOME}/tmp/lsf/echo_log.txt",
        "${HOME}/tmp/lsf/ls_log.txt",
        "${HOME}/tmp/lsf/ls-lh_log.txt",
        ]

    bsub("myArray", commands, logfiles, q='short')
