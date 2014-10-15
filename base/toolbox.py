from socket import gethostname
import StringIO, time
import numpy as np
#import os.path
#from os import system, path, listdir, getenv, environ
from xmltools import ReadConfigXML
import sys, os, re, random
#import shutil
import xml.dom.minidom
#from xml.dom.minidom import Node
#from xml.dom.minidom import getDOMImplementation
from string import split
import ROOT as R
if not sys.flags.interactive: R.gROOT.SetBatch(True)

class whitness(object):
    # this class can be used to add a logging to work that's been
    # done outside the normal realms of logfiles, stdout etc.
    def __init__(self,filename,**kwargs):
        self.filename = filename
        self.__dict__.update(kwargs)
        self.msg = None

    def check_batch(self):
        kstring = None
        LSB_JOBNAME = os.getenv('LSB_JOBNAME')
        LSB_BATCH_JID = os.getenv('LSB_BATCH_JID')
        if (LSB_JOBNAME is None or LSB_BATCH_JID is None):
            kstring = 'type=local '
        else:
            kstring = 'type=batch id=%s task=%s '%(LSB_BATCH_JID,LSB_JOBNAME)
        return kstring
    
    def system_stamp(self):
        return '%s: hostname=%s '%(time.ctime(),gethostname())
        
    def log(self,strg,mode='a'):
        f = open(self.filename,mode)
        self.msg = ''
        if mode=='w':
            self.msg = self.system_stamp()+'*INFO* writing new whitness file %s\n'%self.filename
            print self.msg
        self.msg+= self.system_stamp()+self.check_batch()+' %s'%strg
        f.write(self.msg+'\n')
        f.close()
    
    def done(self):
        f = open(self.filename,'a')
        self.msg = self.system_stamp()+'*INFO* exiting logger, bye'
        f.write(self.msg+'\n')
        f.close()
        
class logger(object):
    def __init__(self,logLevel=2):
        # 0: error 1: warning 2: info 3: debug
        self.level = logLevel
        pass
    def set_level(self,level):
        self.level = level
    
    def _getTimeStamp(self):
        import time
        timeStamp = time.ctime()
        return timeStamp
    
    def info(self,msg):
        t = self._getTimeStamp()
        if self.level>=2: print '*INFO* %s: %s'%(t,msg)

    def debug(self,msg):
        t = self._getTimeStamp()
        if self.level>=3: print '*DEBUG* %s: %s'%(t,msg)

    def warning(self,msg):
        t = self._getTimeStamp()
        if self.level>=1: print '*WARNING* %s: %s'%(t,msg)
        
    def error(self,msg):
        t = self._getTimeStamp()
        ExcText = '*FATAL* %s: %s'%(t,msg)
        raise Exception(ExcText)
##### done ####

DEBUG = False;

def make_envfile(ffile):
    lines = ""
    if os.path.isfile(ffile):
        print '*INFO* overwriting existing file %s'%ffile
        lines = '#!/bin/bash\n'
    for key in os.environ.keys():
        line = 'export %s="%s"\n'%(key,str(os.getenv(key)))
        lines+=line
    lines+='echo "DONE WITH RUNTIME ENV SETUP" \n'
    f = open(ffile,'w')
    f.write(lines)
    f.close()
    os.chmod(ffile,0755)
    print '*INFO* wrote env file %s'%ffile

def mkdir(dir):
    if not os.path.exists(dir):  os.makedirs(dir)
    return dir

def mkscratch():
    if os.path.exists('/scratch/'):    
        return(mkdir('/scratch/%s/'%os.getenv('USER',"zimmer")))
    elif os.path.exists('/tmp/'):
        return(mkdir('/tmp/%s/'%os.getenv("USER","zimmer")))
    else:
        raise Exception('...')

def isBatch():
    job_id = os.getenv("LSB_JOBID","None")
    kret = False
    if job_id != "None": kret = True
    print '*isBatch {}*'.format(kret)
    stop_here()
    return kret

def construction_site():
    raise Exception("** Error: the feature you have requested is not there (yet) **");

def do_nothing():
    print '*INFO* **do nothing**'

def stop_here(exc=None):
    import traceback;
    print '*** traceback of previous errors ***'
    print traceback.print_exc();
    print '*** traceback of previous errors ***'
    exception = '';
    if not exc==None:
        exception+=exc;
    raise Exception("*** STOP HERE ***"+exception);
        
def setPF(d=0):
    #    from os import getenv, system, envrion
    # handle p-files (for multiple job-submission)
    if (d>1):
        print '*INFO* SETTING PFILES'
    #HOME=os.getenv('HOME');
    USER=os.getenv('USER');
    PFILES=os.getenv('PFILES');
    LSB_JOBID=os.getenv('LSB_JOBID','NONE')
    if (LSB_JOBID=='NONE'):
        if (d>0):
            print 'WARNING: you are not running in batch mode, we don\'t change anything! ***'
        os.system('mkdir -p /tmp/'+USER);
        os.system('mkdir -p /tmp/'+USER+'/');
        os.system('mkdir -p /tmp/'+USER+'/pfiles');
        os.environ['HOME']='/tmp/'+USER; #Jim's suggestion
        try:
            PFILESNEW=PFILES.replace('/u/gl/'+USER+'/pfiles;','/tmp/'+USER+'/pfiles;');
        except AttributeError:
            PFILESNEW = '/u/gl/'+USER+'/pfiles' 
        os.environ['PFILES']=PFILESNEW;
    else:
        os.system('mkdir -p /scratch/'+USER);
        os.system('mkdir -p /scratch/'+USER+'/'+LSB_JOBID);
        os.system('mkdir -p /scratch/'+USER+'/'+LSB_JOBID+'/pfiles');
        os.environ['HOME']='/scratch/'+USER+'/'+LSB_JOBID; #Jim's suggestion
        PFILESNEW=PFILES.replace('/u/gl/'+USER+'/pfiles;','/scratch/'+USER+'/'+LSB_JOBID+'/pfiles;');
        os.environ['PFILES']=PFILESNEW;
        # test whether it's okay
    if (d>1):
        print '*INFO* ENVIRONMENT SET UP: '+os.environ['PFILES']+'\n'
        os.system('echo $PFILES');
        os.system('echo $HOME');

def cleanPF(d=0):
    # handle p-files (for multiple job-submission)
    USER=os.getenv('USER');
    LSB_JOBID=os.getenv('LSB_JOBID','NONE')
    if (LSB_JOBID=='NONE'):
        if (d>0):
            print 'WARNING: you are not running in batch mode, we don\'t change anything! ***'
        os.system('rm -rf /tmp/'+USER);

    else:
        if (d>1):
            print '*INFO* CLEAN UP PFILES '
        os.system('rm -rf /scratch/'+USER+'/${LSB_JOBID}');

def ContainsExpr(string,expr):
    ''' returns True if the expression expr is contained in string '''
    try:
        parts = string.partition(expr);
        if (len(parts[1])>0):
            return True;
        else:
            return False;
    except IndexError:
        False;

def runExec(execname,arglist='',d=DEBUG):
    import subprocess, traceback;
    path ='' # leave empty
    execstring= path+execname+arglist
    #print 'run: '+execstring
    rc = os.system(execstring)
    # do this only on slac...
    # okay if the command is not found, abort - that should kill jobs that catch error messages when running on slac batch
    so = subprocess.Popen('domainname', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = so.communicate();
    try:
        stdout = out[0].splitlines();
        stderr = out[1].splitlines();
    except IndexError:
        sys.exit();
    if (len(stderr)>0):
        raise RuntimeError;
    if (ContainsExpr(stdout[0],'slac.stanford.edu')):
        if not (d):
            if not (rc==0):

                # this is an attempt to catch exceptions while running (and not only at the end)
                print '**FATAL: runExec could not run properly'
                print '\n* stacktrace:\n'
                traceback.print_exc()
                print '\n* command: ',execname
                print '\n* arguments: ',arglist;
                print '\n* return code: ',rc;
                print '\n* stderr output: \n',stderr;
                print '\n**ABNORMAL_TERMINATION Good Bye.' 
                os.abort();

def run(cmd,stdout=None,batchargs=None,debug=False):    
    runcmd = ''
    if batchargs is None:
        runcmd = cmd
    else:
        args = ''
        if type(batchargs)==dict:
            for k in batchargs.keys():
                args+='%s %s'%(k,batchargs[k])
            #need to expand the batch args
        else:
            args = batchargs
        if not stdout is None:
            logfile = stdout
            runcmd = 'bsub -R \"rhel60 && linux64\" %s -o %s %s'%(args,stdout,cmd)
        else:
            runcmd = 'bsub -R \"rhel60 && linux64\" %s %s'%(args,cmd)
    print '*EXE* %s'%runcmd
    # implement trapping of bsub calls
    if debug:
        if not batchargs is None:
            jobid = random.randrange(100,900)
            print 'Job <%i> is submitted to default queue <short>'%jobid
            return jobid
        else:
            return
    output = {'returncode':None,'output':None,'error':None}
    if batchargs is None:
        os.system(runcmd)
        return None
    else:
        from subprocess import Popen,PIPE
        proc = Popen(runcmd, stderr=PIPE, stdout=PIPE, shell=True)
        proc.wait() # wait for it to finish                                                                                                      
        if proc.returncode!=0:
            raise RuntimeError("%s\n*ERROR*: error accessing batchsystem, reported RC %i"%(proc.stderr.read(),proc.returncode))
        else:
            output = proc.stdout.read()
            #Job <1848> is submitted to default queue <short>
            try:
                return int(re.search(r"<\d*>",output).group(0).replace("<","").replace(">",""))
            except IndexError:
                raise RuntimeError("*ERROR*: error accessing batchsystem when querying for job-ID here's the output\n%s"%output)
                
        
        
def is_file(fname,d=0):
    '''this subroutine checks the fetch-data-file
    and copies stuff to the designated areas'''
    # check if exist in direction
    if os.path.isfile(fname):
        if (d>2):
            print '*INFO* open file:', fname
        return True
    else:
        return False

def is_DirNotEmpty(fname):
    ''' subroutine checks whether dir exists and is not empty'''
    d = os.path.dirname(fname)
    if not os.path.exists(d):
        print 'FATAL: tag dir does not exist, exiting'
        sys.exit();
    elif (os.path.exists(d)):
        # check if non-empty
        if (os.listdir(fname)==[]):
            print 'FATAL: tag dir exists, but is empty, exiting'
            sys.exit();
        

# main program
def ReadConfigFile(fname):
    if (is_file(fname)):
        pars = [];
        # loop over all lines in file
        File  = open(fname,'r')
        lines = File.readlines()
        for line in lines:
            if not line.startswith("#"):
                linestr = line.rstrip()
                # now get rid of 'comments'
                source = linestr.partition("#");
                psource = source[0];
                par = psource.partition("="); 
                # get rid of line breaks, spaces and tabs
                spar = par[2]
                spar1 = spar.replace(" ","");
                spar2 = spar1.replace("\n","");
                spar3 = spar2.replace("\t","");
                spar4 = spar3.replace("'","");
                pars.append(spar4);
        # check length
        # if this number is not equal the last one, send exception
        NumberOfPars = pars[-1];
        if not (int(NumberOfPars)==len(pars)):
            print 'FATAL: error, check config file'
            print 'FATAL: parameters supplied: '+str(len(pars));
            print 'FATAL: parameters expected: '+NumberOfPars
            sys.exit();
        File.close();
        return pars;

    else:
        print 'FATAL: error, file does not exist'
        sys.exit();
        
def ReadDataFile(fname):
    name = [];
    ra = [];
    dec= [];
    J  = [];
    Junc= [];
    ID = [];

    if (is_file(fname)):
        # loop over all lines in file
        File  = open(fname,'r')
        lines = File.readlines()
        n = 0;
        for line in lines:
            if not line.startswith("#"):
                thisline = line.split();
                name.append(thisline[0]);
                dec.append(thisline[1]);
                ra.append(thisline[2]);
                J.append(thisline[3]);
                Junc.append(thisline[4]);
                ID.append(thisline[5]);
                n=n+1;
    else:
        print 'FATAL: error, file does not exist'
        sys.exit();
    output = [name, dec, ra, J, Junc, ID];
    File.close();
    return output;
    
def GetMassBounds():
    fname = 'dat/yield.dat'
    if (is_file(fname)):
        mass = [];
        # loop over all lines in file
        File  = open(fname,'r')
        lines = File.readlines()
        for line in lines:
            if not line.startswith("#"):
                thisline = line.split();
                mass.append(thisline[0]);
    else:
        print 'error, file does not exist'
        sys.exit();
    mMin = float(mass[0])*0.95; 
    mMax = float(mass[-1])*1.05;
    output = [mMin,mMax];
    File.close();
    return output;

def ReadFile(fname):
    pars = [];
    if (is_file(fname)):
        # loop over all lines in file
        File  = open(fname,'r')
        lines = File.readlines()
        for line in lines:
            if not line.startswith("#"):
                
                linestr = line.rstrip()
                # now get rid of 'comments'
                source = linestr.partition("#");
                psource = source[0];
                par = psource.partition("="); 
                # get rid of line breaks, spaces and tabs
                spar = par[2]
                spar1 = spar.replace(" ","");
                spar2 = spar1.replace("\n","");
                spar3 = spar2.replace("\t","");
                spar4 = spar3.replace("'","");
                pars.append(spar4);
        File.close();
    else:
        print 'FATAL: error, file does not exist'
        sys.exit();

    return pars;

def CreateEmptyTree(configFile):
    from ROOT import TTree, TFile, AddressOf, gROOT
    # read configurations file
    pars = ReadConfigXML(configFile);
    tag = pars['Tag'];
    path_output=pars['HOME']+'/results/'+tag+'/';
    # Make a tree
    f = TFile(path_output+'results.root','RECREATE');
    t = TTree('CompLike','Composite Likelihood Results');

    # Create a struct
    gROOT.ProcessLine(\
    "struct MyStruct{\
      Char_t ClusterName[20];\
      Char_t IRF[20];\
      Char_t tag[20];\
      Double_t mass;\
      Double_t MLE;\
      Double_t PosMinos;\
      Double_t NegMinos;\
      Double_t sigmav;\
      Double_t flux;\
      Double_t GalDiff;\
      Double_t IsoDiff;\
    };")
    from ROOT import MyStruct

    # Create branches in the tree
    s = MyStruct()
    t.Branch('ClusterName',AddressOf(s,'ClusterName'),'ClusterName/C');
    t.Branch('IRF',AddressOf(s,'IRF'),'IRF/C');
    t.Branch('tag',AddressOf(s,'tag'),'tag/C');
    t.Branch('mass',AddressOf(s,'mass'),'mass/D')
    t.Branch('MLE',AddressOf(s,'MLE'),'MLE/D')
    t.Branch('PosMinos',AddressOf(s,'PosMinos'),'PosMinos/D')
    t.Branch('NegMinos',AddressOf(s,'NegMinos'),'NegMinos/D')
    t.Branch('sigmav',AddressOf(s,'sigmav'),'sigmav/D')
    t.Branch('flux',AddressOf(s,'flux'),'flux/D')
    t.Branch('GalDiff',AddressOf(s,'GalDiff'),'GalDiff/D')
    t.Branch('IsoDiff',AddressOf(s,'IsoDiff'),'IsoDiff/D')
    
    f.Write()
    f.Close()
    print '*INFO* ** ROOT File with tree structure successfully created **';

def dumpConfig(configFile):
    ''' dumps the content from the config-file and associated into text file for parsing in C++'''
    pars=ReadConfigXML(configFile);
    newpars=[];
    # what do i need...
    storeDir=pars['HOME']+'/results/'+pars['Tag']+'/';
    HOME=os.getenv('HOME')+'/';
    DumpClusters=HOME+'dump-cluster.dat';
    datafile=DumpClusters; #gives the dumped cluster names
    if (pars['Ind_FluxFile']==''):
        pars['Ind_FluxFile']=HOME+'yield.dat'
    cdate=pars['cdate'];
    if (cdate=='11m'):
        dataSet='11'
    elif (cdate=='18m'):
        dataSet='18'
    elif (cdate=='24m'):
        dataSet='24'
    if bool(pars['C2_TieGal']):
        GalTie='yes';
    elif not bool(pars['C2_TieGal']):
        GalTie='no';
    if bool(pars['C2_TieEgal']):
        IsoTie='yes';
    elif not bool(pars['C2_TieEgal']):
        IsoTie='no';
    
    newpars.append(storeDir);
    newpars.append(datafile);
    newpars.append(pars['Ind_FluxFile']);
    newpars.append(pars['IRF']);
    newpars.append(dataSet);
    newpars.append(GalTie);
    newpars.append(IsoTie);
    # dump in rich text file
    HOME=os.getenv('HOME')+'/';
    target=HOME+'dump-parameters.dat';
    if (int(pars['chatter'])>2):
        print '*INFO* *** dumping parameters to file: '+target+' ***'
    File=open(target,'w');
    i=0;
    for i in range(0,len(newpars)):
        if sbool(pars['Debug']):
            print(newpars[i]);
        else:
            File.write(newpars[i]+'\n');
    File.close();
    # if (os.path.exists(target)):
    #     print '*INFO* config file exists'
    #     os.system('cat '+target);
    # print '*INFO* DumpConfig finished'
    
def DumpClusterNames(configFile):
    ''' dump just cluster-names '''
    pars = ReadConfigXML(configFile);
    fname= pars['ROIFile'];
    cluster = [];
    if (is_file(fname,int(pars['chatter']))):
        # loop over all lines in file
        File  = open(fname,'r')
        lines = File.readlines()
        for line in lines:
            if not line.startswith("#"):
                thisline = line.split();
                cluster.append(thisline[0]);
    if bool(pars['C2_DoMINOS']):
        cluster.append('compLike');
    File.close();
    # now dump the content...
    HOME=os.getenv('HOME')+'/';
    target=HOME+'dump-cluster.dat';
    if (int(pars['chatter'])>2):
        print '*INFO* *** dumping cluster names to file: '+target+' ***'
    File=open(target,'w');
    i=0;
    for i in range(0,len(cluster)):
        File.write(cluster[i]+'\n');
    # check if we do composite and add if necessary
    file.close();
#    print '*INFO* DumpClusterNames finished'

def sbool(string):
    ''' helper that converts python bool in stuff that scienceTools can read'''
    if (string=='True') or (string=='yes'):
        return True;
    else:
        return False;

def hbool(b):
    ''' returns a yes or no for true boolean'''
    if (b):
        return 'yes'
    else:
        return 'no';

def convBool(string):
    if string=='True':
        return True
    else:
        return False

class ROI(object):
    def __init__(self,d): 
        self.radius = None;
        self.__dict__.update(d);
        self.list_of_objects = [];

    def get_content(self):
        return self.__dict__;
    def add_object(self,astroobj):
        # add an astroobject to list of objects in roi
        self.list_of_objects.append(astroobj);
    def rm_object(self,obj_name):
        kReturn = False;
        for obj in self.list_of_objects:
            if obj.name==obj_name:
                try:
                    self.list_of_objects.remove(obj);
                    kReturn = True;
                except ValueError:
                    raise Exception("Found name in object list but cannot remove from list, check your data");
        if not kReturn:
            raise Exception("Could not find any matching object,",obj_name," abort.");

class astroobj(object):
    def __init__(self,d): self.__dict__.update(d);
    def get_content(self):
        return self.__dict__;
    def update(self,dkey,dvalue):
        # to update single fields with one method - if a key is not there, it will be appended to the class
        keys = self.__dict__.keys();
        for k in keys:
            if k==dkey:
                self.__dict__[k]=dvalue;
            else:
                self.__dict__[dkey]=dvalue;
        
class DataSet(object):
    def __init__(self,fname):
        # this one takes the datafile and provides the interface for the sub scripts
        self.list_of_rois = [];
        self.xmlfile = fname;
        import xml.dom.minidom
        if not (is_file(self.xmlfile)):
            raise Exception('FATAL: cannot open',fname);
        xmldoc = xml.dom.minidom.parse(self.xmlfile);
        for element in xmldoc.getElementsByTagName("ROI"):
            d = {};
            Keys = element.attributes.keys();
            #print Keys;
            for key in Keys:
                thisKey = str(key);
                value = str(element.attributes[key].value);
                d[thisKey]=value;
                #print d;
            roi = ROI(d);
            roi.radius = element.getAttribute("radius");
            for obj in element.getElementsByTagName("object"):
                # loop over all objects in each ROI that we want to include
                od = {};
                Keys = obj.attributes.keys();
                for key in Keys:
                    #print key
                    od[str(key)]=str(obj.attributes[key].value)
                astro = astroobj(od);
                roi.add_object(astro);
                del astro;
            self.list_of_rois.append(roi);
            del roi;
                
    def get_names(self):
        kReturn = [];
        for roi in self.list_of_rois:
            kReturn.append(roi.name);
        return kReturn;
    
    def get_positions(self):
        RA = [];
        DEC = [];
        for roi in self.list_of_rois:
            RA.append(roi.RA);
            DEC.append(roi.DEC);
        kReturn = [RA,DEC];
        return kReturn;

    def update(self,outputxml=None):
        # this method is used to write changes to the output xml
        import xml.dom.minidom
        xmlstring = '';
        if not (is_file(self.xmlfile)):
            raise Exception('FATAL: cannot open',fname);
        xmldoc = xml.dom.minidom.parse(self.xmlfile);
        for element in xmldoc.getElementsByTagName("ROI"):
            for roi in self.list_of_rois:
                if roi.name == element.getAttribute("name"):
                    # we have a match
                    d = roi.__dict__;
                    for key in d.keys():
                        element.setAttribute(key,str(d[key]));
                        #print 'key:',key,' value:',d[key];
                        if key=='list_of_rois':
                            pass;
                            #element.setAttribute(key,str(self.get_names()));
                    xmlstring+=element.toxml()+'\n';
        firstLine='<?xml version="1.0" ?>\n<list_of_roi title="my rois">'
        lastLine='</list_of_roi>'
        if not outputxml==None:
            f = open(outputxml,'w');
        else:
            f = open(self.xmlfile,'w');
        f.write(firstLine+'\n');
        f.write(xmlstring+'\n');
        f.write(lastLine);
        f.close();

    def selectRois(self,rois,sep=","):
        tlist = None
        tlist = rois.split(",")
        print tlist
        new_list = []
        for i,obj in enumerate(self.list_of_rois):
            if obj.name in tlist:
                # okay --- now need to check that the names are identical
                for item in tlist:
                    if item == obj.name and len(item)==len(obj.name):
                        new_list.append(obj)
        self.list_of_rois = new_list
        #print [roi.name for roi in self.list_of_rois]
        #sys.exit(1)
        return self.list_of_rois
        
def GetCoordinatesFromFitsFile(fitsfile):
    ''' uses pyfits to extract information about the roi center from a fitsfile '''
    try:
        import pyfits
    except IOError:
        raise Exception('FATAL: cannot import pyfits, cannot use this method');
    #print 'Use Fitsfile for coordinate extraction',fitsfile;
    hdu = pyfits.open(fitsfile);
    hduhdr = hdu[0].header;
    try:
        RA = hduhdr['CRVAL1'];
        DEC = hduhdr['CRVAL2'];
    except KeyError:
        raise Exception('FATAL: Cannot find the keywords for RA and DEC, check the fitsfile',fitsfile);
    return RA,DEC;

class InterpPowerLaw(object):
    """ Class for log interpolation. """
    def __init__(self, x, y):
        self.xarr, self.yarr = x, y
        self.log_xarr = np.log10(self.xarr)
        self.log_yarr = np.log10(self.yarr)
        
    def __call__(self, x):
        log_x = np.log10(x)
        i = np.where(self.log_xarr <= log_x)[0][-1]
        low , hi = self.log_xarr[i], self.log_xarr[i+1]
        ldata, hdata = self.log_yarr[i], self.log_yarr[i+1]
        return 10**(ldata + (hdata-ldata)*(log_x-low)/(hi-low))
    


def integrate(f,emin=None,emax=None):
    from scipy.integrate import quad
    
    # just does some integration
    energies, dnde = np.loadtxt(f,unpack=True)
    interp = InterpPowerLaw(energies,dnde)
    elow = energies[0]
    ehigh = energies[-1]
    if not emin is None:
        elow = float(emin)
    if not emax is None:
        ehigh = float(emax)

    integral = quad(interp, elow, ehigh, full_output=True)[0]
    return integral

def resolveFilename(fstr):
    is_envvar = False
    if fstr.startswith("$"):
        is_envvar = True
    if not is_envvar:
        if os.path.isfile(fstr):
            return fstr
        else:
            raise Exception("*ERROR* cannot find %s"%str)
    else:
        env_var = fstr.replace("$(","").replace(")","")
        if os.getenv(env_var)==None:
            raise Exception("*ERROR* filename appears to be defined as Env VAR but cannot find %s"%env_var)
        else:
            if os.path.isfile(os.getenv(env_var)):
                return os.getenv(env_var)
            else:
                raise Exception("ERROR* cannot find %s"%os.getenv(env_var))

def parse_options(optparse=None,sep=";"):
    # with a valid optparse entry, return a proper dict
    if optparse is None:
        return None
    else:
        config_args = {}
        pairs = optparse.split(sep)
        for p in pairs:
            the_pair = p.split("=")
            if len(the_pair)==2:
                config_args[the_pair[0]]=the_pair[1]
            else:
                print "*error* could not translate arg %s"%p
        return config_args

def make_covariance_plot(covar_mat):
    # given a covariance matrix, returns a TH2D object containing the matrix plotted
    int_mat = np.array(covar_mat)
    x_dim,y_dim = np.shape(int_mat)
    h2 = R.TH2D("h2","covariance matrix",x_dim,0,x_dim,y_dim,0,y_dim)
    R.gROOT.SetStyle("Plain")
    R.gStyle.SetOptStat(00000)
    #R.gStyle.SetPalette(51)
    print R.TColor.__dict__
    sys.exit()
    h2.SetContour(25)
    for index,value in np.ndenumerate(int_mat):
        #print index,value
        h2.SetBinContent(index[0]+1,index[1]+1,value)
        x_label = str(index[0]+1)
        y_label = str(index[1]+1)
        if x_dim>=10:
            x_label = ""
            if not index[0]%2==0:
                x_label = str(index[0]+1)
        if y_dim>=10:
            y_label = ""
            if not index[1]%2==0:
                y_label = str(index[1]+1)
        h2.GetXaxis().SetBinLabel(index[0]+1,x_label)
        h2.GetYaxis().SetBinLabel(index[1]+1,y_label)
    h2.GetXaxis().SetTitle("MINUIT FIT ID")
    h2.GetXaxis().CenterTitle()
    h2.GetYaxis().SetTitle("MINUIT FIT ID")
    h2.GetYaxis().CenterTitle()
    return h2

def checkRunningLSFJobs(joblist=None):
    if joblist is None:
        joblist=[]
    if not isinstance(joblist,list):
        raise Exception("must be a python list!")
    # use matching to extract the list of running jobs                                                                                       
    all_done = True
    from subprocess import Popen, PIPE
    proc = Popen("bjobs", stderr=PIPE, stdout=PIPE, shell=True)
    proc.wait() # wait for it to finish                                                                                                     \
                                                                                                                                             
    if proc.returncode!=0:
        raise RuntimeError("%s\n*ERROR*: error accessing batchsystem, reported RC %i"%(proc.stderr.read(),proc.returncode))
    output = proc.stdout.readlines()
    jobsInLSF = []
    for line in output:
        thisLine = line.split()
        if len(thisLine)!=0:
            try:
                jobsInLSF+= [int(thisLine[0])]
            except ValueError:
                pass
    # match against existing list:                                                                                                           
    for job in joblist:
        if job in jobsInLSF:
            all_done = False
    return all_done

def mkdir(dir):
    if not os.path.exists(dir):  os.makedirs(dir)
    return dir

def rmdir(dir):
    import shutil
    shutil.rmtree(dir)

