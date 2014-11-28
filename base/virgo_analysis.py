from base.likelihood import *
import base.toolbox
from base.xmltools import GetParameterIDForLikelihood
import xml.dom.minidom as xdom
from tempfile import NamedTemporaryFile
from skymaps import SkyDir
import numpy as np
import copy, sys, os, pickle, time, yaml
import lsf

def effectiveJ(mass,finalState="ee"):
    ratio_ee = {80.0: 6870.8334694828354, 100.0: 7857.3235653271013, 5.0: 1.0, 1000.0: 8796.6830234728004, 
                300.0: 9536.9635653948353, 10.0: 1.0, 12.0: 1.0001583731779042, 500.0: 9358.6246594617005, 
                15.0: 1.1489987912435411, 2000.0: 8120.2511392272336, 200.0: 9372.6704750819463, 
                50.0: 3922.6097800051202, 20.0: 24.791874384337806, 40.0: 2344.9138681816676, 
                25.0: 229.70062960201443, 30.0: 751.01037167073275}
    ratio_mu = {80.0: 2061.1646238830808, 100.0: 3075.946719746662, 5.0: 1.0, 1000.0: 10655.988532241725, 
                300.0: 8209.7408534267888, 10.0: 1.0, 12.0: 1.0, 500.0: 9672.8832784911829, 
                15.0: 1.000999814296067, 2000.0: 10768.481332069472, 200.0: 6502.5240241749752, 
                50.0: 551.14539383992917, 20.0: 1.3422163985180275, 40.0: 209.95731162025731, 
                25.0: 6.9740449704918319, 30.0: 33.672956505908054}
    return 1

eeIC = {20:'ee_IC_mx20GeV.dat',15:'ee_IC_mx15GeV.dat',25:'ee_IC_mx25GeV.dat',300:'ee_IC_mx300GeV.dat',
        500:'ee_IC_mx500GeV.dat',100:'ee_IC_mx100GeV.dat',5:'ee_IC_mx5GeV.dat',50:'ee_IC_mx50GeV.dat',
        12:'ee_IC_mx12GeV.dat',80:'ee_IC_mx80GeV.dat',10:'ee_IC_mx10GeV.dat',2000:'ee_IC_mx2000GeV.dat',
        200:'ee_IC_mx200GeV.dat',1000:'ee_IC_mx1000GeV.dat',40:'ee_IC_mx40GeV.dat',30:'ee_IC_mx30GeV.dat'}

mumuIC = {20:'mumu_IC_mx20GeV.dat',15:'mumu_IC_mx15GeV.dat',25:'mumu_IC_mx25GeV.dat',300:'mumu_IC_mx300GeV.dat',
          500:'mumu_IC_mx500GeV.dat',100:'mumu_IC_mx100GeV.dat',5:'mumu_IC_mx5GeV.dat',50:'mumu_IC_mx50GeV.dat',
          12:'mumu_IC_mx12GeV.dat',80:'mumu_IC_mx80GeV.dat',10:'mumu_IC_mx10GeV.dat',2000:'mumu_IC_mx2000GeV.dat',
          200:'mumu_IC_mx200GeV.dat',1000:'mumu_IC_mx1000GeV.dat',40:'mumu_IC_mx40GeV.dat',30:'mumu_IC_mx30GeV.dat'}

channel_mapping = {
    1  :  ["e+e-","ee"]                  ,
    2  :  ["mu+mu-","mumu","musrc"]      ,
    3  :  ["tau+tau-","tautau","tausrc"] ,
    4  :  ["bb-bar","bb","bbbar","bbar","bbsrc"],
    5  :  ["tt-bar","tt"]                ,
    6  :  ["gluons"]                     ,
    7  :  ["W+W-","w+w-","ww","wwsrc"]   ,
    8  :  ["ZZ","zz"]                    ,
    9  :  ["cc-bar","cc"]                ,
    10 :  ["uu-bar","uu"]                ,
    11 :  ["dd-bar","dd"]                ,
    12 :  ["ss-bar","ss"]                ,
    }

def setPF(dir):
    right   = os.environ['INST_DIR']+'/syspfiles'
    newPF = dir+";"+right
    os.environ["PFILES"]=newPF
    
def resetPF():
    right   = os.environ['INST_DIR']+'/syspfiles'
    left = "/u/gl/zimmer/pfiles"
    newPF = left+";"+right
    os.environ["PFILES"]=newPF
    
def channel2int(s):
    for k,v in channel_mapping.items():
        if s in v: return k
    else:  raise ValueError("Can't find value %s"%s)

def make_virgoSource(J=1e17,profile="NFW",outdir="/tmp/",template=None):
    #print '*virgo %s'%profile
    from uw.like.SpatialModels import Disk, SpatialMap
    from uw.darkmatter.spatial import NFW
    from uw.darkmatter.spectral import DMFitFunction
    from uw.like.roi_extended import ExtendedSource
    # returns the xmlsource for virgo with correct J
    if template is None:
        skydir = SkyDir(187.7006531,12.40596008,SkyDir.EQUATORIAL)
        rvir = 6.2 # use 6.2 degrees for r200
        if profile == "NFW":
            spatial_model = NFW(sigma=rvir,center=skydir)
        elif profile == "flat":
            spatial_model = Disk(sigma=rvir,center=skydir)
        template = os.path.join(outdir,"Virgo_%s.fits"%profile)
        spatial_model.save_template(template)
    
    # need to write that out at some point
    model = DMFitFunction(norm=J, sigmav=1e-26, channel0=4,channel1=1,mass=100,
                          bratio=1.0,file="$(INST_DIR)/Likelihood/data/gammamc_dif.dat")
    ps1 = ExtendedSource(
        name = 'Virgo',
        model = model,
        spatial_model=SpatialMap(file=template)
        )
    from uw.utilities.xml_parsers import unparse_diffuse_sources
    # only interested in the blurb 
    xml_blurb = unparse_diffuse_sources([ps1],strict=False,
                                        convert_extended=True,
                                        #extended_dir_name=outdir,
                                        expand_env_vars=False,filename=None)
    x = xdom.parseString(xml_blurb[-1]).firstChild
    [par.setAttribute("free","0") for par in x.getElementsByTagName("parameter") if not par.getAttribute("name")=="sigmav"]
    Jval = ("%1.5e"%float(J)).partition("e")
    for p in x.getElementsByTagName("parameter"):
        if p.getAttribute("name")=="norm":
            p.setAttribute("value","%1.5f"%float(Jval[0]))
            p.setAttribute("scale","1e%s"%Jval[-1])
            p.setAttribute("min","0.0001")
            p.setAttribute("max","10000.0")
        if p.getAttribute("name") in ["channel0","channel1"]:
            p.setAttribute("min","1")
            p.setAttribute("max","13")
        if p.getAttribute("name")=="bratio":
            p.setAttribute("min","0.")
            p.setAttribute("max","1.")
    return x

def run_prepare(roi,dry=True):
    print 'GAL: %s\nEGAL: %s'%(os.getenv("GAL"),os.getenv("EGAL"))
    #only does the srcmaps part
    from GtApp import GtApp
    srcMap = GtApp('gtsrcmaps', 'Likelihood')
    srcMap['irfs']=roi.configuration.IRF
    srcMap['scfile']=roi.files['ft2file']
    srcMap['expcube']=roi.files['expcube']
    srcMap['srcmdl']=roi.modelxml
    srcMap['bexpmap']=roi.files['binnedExpMap']
    srcMap['outfile']=roi.files['srcmap']
    srcMap['emapbnds']='no'
    #srcMap['npix']=1000 # see https://confluence.slac.stanford.edu/x/8xYpC
    srcMap['minbinsz']=0.05
    srcMap['cmap']=roi.files['ccube']
    srcMap.run(dry_run=dry)
    return srcMap.command()

def _setup(templates = None,rundir = "/afs/slac/g/glast/users/zimmer/FermiData/Clusters/Virgo_v3/"):
    
    #run_prepare(None)
    # input files
    jogler_path = "/nfs/farm/g/glast/u/jogler/data_for_stephan/"
    files = {
        'expcube':os.path.join(jogler_path,"ltcube_239557414_335188802.fits"),
        'binnedExpMap':os.path.join(jogler_path,"expcube_239557414_335188802_100-100000_P7SOURCE_V6.fits"),
        'ccube':os.path.join(jogler_path,"ccube_M87_GAL_0.1dpb_AIT_100-100000MeV_P7SOURCE_V6.fits"),
        'cmap':os.path.join(jogler_path,"cmap_M87_GAL_0.1dpb_AIT_100-100000MeV_P7SOURCE_V6.fits"),
        'ft2file':os.path.join(jogler_path,"virgo_merged-ft2-30s.fits")
        }

    ra = 187.7006531
    dec= 12.40596008
    
    # baseXml
#    baseXml = os.path.join(jogler_path,"M87_add_ps1_ps4_ps5.xml")
    baseXml = os.path.join(jogler_path,"double_disk_template.xml") # that's the double disk patch
    #baseXml = os.path.join(jogler_path,"M87_add_ps1_ps4_ps5.xml")
    my_run_dir = rundir 
    os.system("mkdir -p %s/log"%my_run_dir)
    config_dict = {}
    config_dict["IRF"]="P7SOURCE_V6"
    config_dict["EMin"]=100
    config_dict["EMax"]=1e5
    config_dict["optimizer"]="Minuit"
    config_dict["my_run_dir"]=my_run_dir
    config_dict["skydir"]=SkyDir(ra,dec,SkyDir.EQUATORIAL)
    config_dict["baseXml"]=baseXml
    config_dict["files"]=files
    config_dict["nullhypothesis"]=False
    config_dict["resultsfile"]=None
    config_dict["PhysicsMode"]="DM_A" # annhihilating DM
    config_dict["FinalState"]="bbar"
    config = configuration()
    config.update(config_dict)
    rois = []
    
    if templates is None:
        raise Exception("must be a not NONE list")
    for template in templates:
        roi = ROI(\
            name = template.name,
            RA = config.skydir.ra(),
            DEC= config.skydir.dec(),
            configuration = config,
            modelxml = template.Xml,
            roifile = None,
            eps = 1e-5,
            gal = template.getGalactic(),
            egal= template.getIsotropic(),
            extended = False)
        virgo_source = source(name="Virgo",parameter="sigmav",is_tied=0,
                              norm=template.getJ(),norm_scale=1,norm_error=0,bratio=1.0,
                              mass=100.,channel0=channel2int("bb"),channel1=1.0)
        virgo_source.DMFitSource()
        roi.list_of_sources+=[virgo_source]
        roi.sourceNames+=["Virgo"]
        srcmap = os.path.join(config.my_run_dir,"BinnedSrcMap_%s.fits"%roi.name)
        files = copy.deepcopy(config.files)
        files["srcmap"]=srcmap
        roi.files=copy.deepcopy(files)
        rois.append(roi)
    return rois
    
def make_prepare(roi,create_virgo=True,dry=True,template_dict=None):    
    if template_dict is None:
        raise NotImplementedError
    #if not template_dict.getIsotropic() is None:
    os.environ["EGAL"]=template_dict.getIsotropic()
    os.environ["GAL"]=template_dict.getGalactic()
    print os.getenv("GAL"), os.getenv("EGAL")
    #if not template_dict.getGalactic() is None:
    diff_names = ["iso_p7v6source","gal_2yearp7v6_v0"]        
    if create_virgo:
        x = xdom.parse(roi.configuration.baseXml)
        for src in x.getElementsByTagName("source"):
            if not src.getAttribute("name") in diff_names:
                #for p in src.getElementsByTagName("parameter"):
                #   p.setAttribute("free","0")
                pass
            else:
                spec = src.getElementsByTagName("spectrum")[0]
                spat = src.getElementsByTagName("spatialModel")[0]
                if src.getAttribute("name") == "iso_p7v6source":
                    src.setAttribute("name","EGAL")
                    spec.setAttribute("file","$(EGAL)")
                else:
                    src.setAttribute("name","GAL")
                    spat.setAttribute("file","$(GAL)")
        virgo = make_virgoSource(J=template_dict.getJ(),profile="NFW",outdir=roi.configuration.my_run_dir,template=template_dict.getFits())
        x.childNodes[-1].appendChild(virgo)
        f = open(roi.modelxml,'w')
        f.write(x.toxml())
        f.close()
    cmd = run_prepare(roi)
    cmds = "from base.toolbox import *\nimport os,sys\nsetPF()\nprint \"GAL:\",os.getenv('GAL'), \"EGAL:\",os.getenv('EGAL') \nprint '**** DEBUG XML *****'"
    cmds += "\nos.system('cat %s')\nprint '****** running gtsrcmaps ******'\nos.system('%s')\ncleanPF()"%(roi.modelxml,cmd)
    cmd_file = NamedTemporaryFile(prefix=os.getenv("MYTMPDIR")+"/").name
    f = open(cmd_file,"w")
    f.write(cmds)
    job = toolbox.run("python %s"%cmd_file,stdout=os.path.join("log/","SrcMaps_%s.log"%roi.name),batchargs="-q xlong",debug=dry)
    return job

norm = lambda mass,j : j/(4*np.pi*2*np.power(mass,2))

def prepareIC(roi,j,mass,final_state,tmpdir,specDir):
    if not final_state in ["ee","mumu"]:
        raise Exception("IC only implemented for ee, mumu")
    ### assembles the necessary file spectrum for IC emission! ###
    # first bit: build sources
    from uw.utilities.xml_parsers import parse_sources, write_sources
    from uw.like.scalemodels import ScaleFactorFileFunction
    ps,ds=parse_sources(roi.modelxml) # requires to pass the sources into the pointlike format!
    # get rid of virgo in the source model
    # is in the ds section
    new_ds = [s for s in ds if not s.name == "Virgo"]
    virgo = [s for s in ds if s.name == "Virgo"][0]
    # now modify the virgo source
    # that way the scale factor can be identified as <sigmav>
    filename = None
    if final_state in ["ee","mumu"]:
        if final_state == "ee":
            # e+/e- final state
            if not mass in eeIC.keys():
                raise Exception("Could not find mass point")
            filename = os.path.join(specDir,eeIC[mass])
        elif final_state == "mumu":
            # mu+/mu- final state
            if not mass in mumuIC.keys():
                raise Exception("Could not find mass point")
            filename = os.path.join(specDir,mumuIC[mass])
        j*=effectiveJ(mass,finalState=final_state) # include the IC correction
    if not os.path.isfile(filename):
        raise IOError("could not find file spectrum %s"%filename)
    if filename is None:
        raise Exception("could not find associated file")
    model = ScaleFactorFileFunction(ScaleFactor=1,normalization=norm(mass,j),file=filename)
    model.set_limits("Normalization",0,1e30)
    model.set_limits("ScaleFactor",0,1e8)
    model.setp(0,1)
    model.setp(1,norm(mass,j))
    model.free = np.array([True,False])
    virgo.model = model # replace DMFit model with file spectrum
    ds = new_ds 
    ds.append(virgo)
    ofile = os.path.join(tmpdir,"%s_%s_%sGeV.xml"%(roi.name,final_state,"%1.0f"%mass))
    write_sources(ps,ds, ofile, strict=True,
                  convert_extended=False,
                  extended_dir_name=None,
                  expand_env_vars=True)
    print '*INF* wrote temporary xml with IC: %s'%ofile
    return ofile

def process_likelihood(roi,configuration,mass_point,j=1.3e18,cl=2.71,minos=True,scan=False,scan_min=0,update_yaml=False,
                       scan_max=None,scan_npts=20,force_srcmap = True,std_diffuse=True,yamlfile=None,IC=True):
    if roi.__dict__.has_key("egal"):
        os.environ["EGAL"]=roi.egal
    if roi.__dict__.has_key("gal"):
        os.environ["GAL"]=roi.gal
    print '*INFO* entering process_likelihood at %s'%str(time.ctime())
    #print configuration()
    #sys.exit()
    finalstate = configuration.FinalState
    llh_scan = None
    par = "sigmav"
    tmpdir = configuration.getTempDir()
    if os.getenv("LSB_JOBID") is None:
        tmpdir = tmpdir.replace("scratch","tmp") # running locally
    if finalstate in ['ee','mumu']:
        # switch on demand to IC!
        specDir = os.path.join(os.getenv("PWD"),"spectra_v2")
        roi.modelxml = prepareIC(roi,j,mass_point,finalstate,tmpdir,specDir)
        par = "ScaleFactor"
        # USE local model
        tmpName = os.path.basename(roi.modelxml).replace(".xml",".fits")
        roi.files['srcmap']=os.path.join(tmpdir,tmpName)
        run_prepare(roi,dry=False) # force re-running gtsrcmaps # MUST re-run w/ IC!
    # enforce recreation of gtsrcmaps    
    if force_srcmap:
        tmpName = os.path.basename(roi.modelxml).replace(".xml",".fits")
        roi.files['srcmap']=os.path.join(tmpdir,tmpName)
    if not os.path.isfile(roi.files['srcmap']):
        run_prepare(roi,dry=False) # force re-running gtsrcmaps
    roi._prepare()
    print '*INFO* setup parameters for channel %s'%finalstate
    print '*INFO* resultsfile: %s'%configuration.resultsfile
    print '*INFO* GAL %s'%os.getenv("GAL")
    print '*INFO* EGAL %s'%os.getenv("EGAL")
    #sys.exit()
    src = roi.list_of_sources[-1]
    src.parameter = par # override the ROI setup
    # compile dictionary that we need to set
    d = dict(zip(["channel0","mass"],[channel2int(finalstate),float(mass_point)]))
    Id = roi.likelihood_fcn.par_index(src.name,par)
    roi.likelihood_fcn[Id].parameter.setBounds(0,1e12)
    roi.likelihood_fcn[Id].parameter.setScale(1e-26)
    if finalstate in ['ee','mumu']:
        # just set the norm!
        Id = roi.likelihood_fcn.par_index(src.name,"Normalization")
        jnorm = norm(mass_point,j)
        jnorm*= effectiveJ(mass_point,finalState=finalstate) # here we set the J-factor
        jval = ("%1.6e"%jnorm).split("e")
        val = float(jval[0])
        scale = float("1e%s"%jval[1])
        roi.likelihood_fcn[Id].parameter.setScale(scale)
        roi.likelihood_fcn[Id].parameter.setBounds(0,1e40) # let's be sure...
        roi.likelihood_fcn[Id].parameter.setValue(val)
        #roi.minos_id = 19 # hack it in... # or 23 for the point source model
    else:
        # not IC, change stuff to be within bounds!
        #roi.minos_id = 23 # hack it in... # or 19 for the double disk
        for key in d:
            Id = roi.likelihood_fcn.par_index(src.name,key)
            if key == "mass":
                mass_min,mass_max = roi.likelihood_fcn[Id].parameter.getBounds()
                if d[key] <= mass_min:
                    mass_min = 0.1*d[key]
                if d[key]>= mass_max:
                    mass_max = 1.9*d[key]
                roi.likelihood_fcn[Id].parameter.setBounds(mass_min,mass_max)
            roi.likelihood_fcn[Id].parameter.setValue(d[key])

            Id = roi.likelihood_fcn.par_index(src.name,"norm")
            val = str("%1.3e"%j).partition("e")

            print roi.likelihood_fcn[Id].parameter.getBounds()

            roi.likelihood_fcn[Id].parameter.setValue(1.)
            roi.likelihood_fcn[Id].parameter.setBounds(0.01,100)
            roi.likelihood_fcn[Id].parameter.setValue(float(val[0]))
            roi.likelihood_fcn[Id].parameter.setScale(float("1e%s"%val[-1]))
    print '*INFO* running with these parameters'
    roi.verifyByEye()
    print '*INFO* entering fit at %s'%str(time.ctime())
    roi.fit(mysource=src)
    minNeg = None
    minPos = None
    llh_scan = None
    # now need minos
    if scan:
        minos = False # disable minos!
        print '*INFO* requested SCAN of LLH rather than MINOS, skip MINOS'
    if minos:
        minos_id = None
        freePars = {i:p.getName() for i,p in enumerate(roi.likelihood_fcn.params()) if p.isFree()}
        for i,p in enumerate(sorted(freePars)):
            if freePars[p]==par:
                minos_id = i
        if minos_id is None:
            raise RuntimeError("could not find MINOS ID")
        print '*INFO* minos ID: ',minos_id
        try:
            minNeg,minPos = roi.minuit_object.Minos(minos_id,cl)
        except RuntimeError:
            print '*ERROR* could not determine MINOS, but we bravely carry on!'
        src.expand({"minos":(minPos,minNeg),"mass":mass_point})
    if scan:
        sc_min = scan_min
        sc_max = scan_max

        if not minNeg is None:
            sc_min = minNeg

        if not minPos is None:
            sc_max = minPos
        # do the scan
        print '*INFO* about to do SCAN'
        llh_scan = roi.likelihood_fcn.scan(src.name,par, xmin=sc_min, xmax=sc_max, npts=scan_npts,tol=1e-10, optimizer=configuration.optimizer, optObject=roi.minuit_object)

    # last not least, need store stuff
    out = {}
    out["fitResultXml"]=roi.exportFitResultToDict() # store stuff as dict instead of xml!
    out['mass']=mass_point
    Id = roi.likelihood_fcn.par_index(src.name,par)
    out['sigmav']={'mle':roi.likelihood_fcn[Id].parameter.getValue(),
                   'scale':roi.likelihood_fcn[Id].parameter.getScale(),
                   'Ts':roi.likelihood_fcn.Ts(src.name,reoptimize=True),
                   'Npred':roi.likelihood_fcn.logLike.NpredValue(src.name)}

    if minos:
        out['sigmav'].update({'pos':minPos,'neg':minNeg})
    if scan:
        out['sigmav'].update({'scan':llh_scan})
    print '*** done with fitting, calculating Npred for diffuse sources ***'

    if std_diffuse:
        # calculate TS values for diffuse components
        Id = roi.likelihood_fcn.par_index("GAL","Prefactor")
        out['GAL']={'Prefactor':roi.likelihood_fcn[Id].parameter.getValue()*roi.likelihood_fcn[Id].parameter.getScale(),
                    #'Ts':roi.likelihood_fcn.Ts("GAL",reoptimize=True),
                    'Npred':roi.likelihood_fcn.logLike.NpredValue("GAL")}
        Id = roi.likelihood_fcn.par_index("GAL","Index")
        out['GAL']["Index"]=roi.likelihood_fcn[Id].parameter.getValue()*roi.likelihood_fcn[Id].parameter.getScale()
    
    Id = roi.likelihood_fcn.par_index("EGAL","Normalization")
    out['EGAL']={'Normalization':roi.likelihood_fcn[Id].parameter.getValue()*roi.likelihood_fcn[Id].parameter.getScale(),
                #'Ts':roi.likelihood_fcn.Ts("EGAL",reoptimize=True),
                'Npred':roi.likelihood_fcn.logLike.NpredValue("EGAL")}
    out['CL']=cl
    print '*** calculating current val of LLH ***'
    out['llh']=roi.likelihood_fcn()
    print '*** done with fits, storing results ***'
    if not llh_scan is None:
        #if len(llh_scan)==2:
        #    llh_scan[1]*=out['llh']
        out['llh_scan']=llh_scan
    out['FinalState']=finalstate
    out['ROI']=roi.name
    out['STOOLS']="ST-%s"%os.getenv("INST_DIR").split("/")[-1]
    
    out["llh0"],out["xmlNull"]=roi.fitNull(export_fit=True) # now returns a tuple!
    
    out["xmlNull"]=roi.exportFitResultToDict()
    print '*INFO* done with likelihood at %s'%str(time.ctime())
    sleeptime = 15
    print '*INFO* sleeping for %i secs to avoid stressing disk'%sleeptime
    time.sleep(sleeptime)
    f = open(configuration.resultsfile,'rb')
    dumpstring = f.read()
    d = pickle.loads(dumpstring)
    d.update({mass_point:out}) # update old dict
    # could also try to write yaml
    if yamlfile is None:
        yamlfile = configuration.resultsfile.replace(".pkl",".yaml")
    dout = {}
    if os.path.isfile(yamlfile):
        print '*INFO* results yaml file exists already. Appending results - use force=True to start anew.'
        if update_yaml:
            dout = yaml.load(open(yamlfile,'rb'))
            print '*INFO* found following datapoints in yamlfile: {}'.format(dout.keys())
    dout.update(d)
    yaml.dump(dout,open(yamlfile,"wb"))
    # now write back
    d_pick = pickle.dumps(d,-1)    
    fo = open(configuration.resultsfile.replace(".pkl",".out"),'w')
    fo.write(str(d))
    fo.close()
    f = open(configuration.resultsfile,"wb")
    f.write(d_pick)
    f.close()
    del d_pick
    print '*INFO* done at %s'%str(time.ctime())

def process_likelihoodCR(roi,configuration,model="file",j=1.3e18,cl=2.71,minos=True,scan=False,scan_min=0,update_yaml=False,
                         scan_max=None,scan_npts=20,force_srcmap = False,std_diffuse=True,yamlfile=None):
    source = "virgo"
    par = "Normalization"
    if model=="powerlaw2":
        raise Exception("not implemented!")
    if roi.__dict__.has_key("egal"):
        os.environ["EGAL"]=roi.egal
    if roi.__dict__.has_key("gal"):
        os.environ["GAL"]=roi.gal
    print '*INFO* entering process_likelihood at %s'%str(time.ctime())
    #print configuration()
    #sys.exit()
    finalstate = "CR model"
    llh_scan = None
    tmpdir = configuration.getTempDir()
    if os.getenv("LSB_JOBID") is None:
        tmpdir = tmpdir.replace("scratch","tmp") # running locally

    # enforce recreation of gtsrcmaps    
    if force_srcmap:
        tmpName = os.path.basename(roi.modelxml).replace(".xml",".fits")
        roi.files['srcmap']=os.path.join(tmpdir,tmpName)
    if not os.path.isfile(roi.files['srcmap']):
        run_prepare(roi,dry=False) # force re-running gtsrcmaps
    roi._prepare()
    print '*INFO* setup parameters for channel %s'%finalstate
    print '*INFO* resultsfile: %s'%configuration.resultsfile
    print '*INFO* GAL %s'%os.getenv("GAL")
    print '*INFO* EGAL %s'%os.getenv("EGAL")
    #sys.exit()
    src = roi.list_of_sources[-1]
    src.set("parameter",par) # override the ROI setup
    # compile dictionary that we need to set
    print '*INFO* running with these parameters'
    roi.verifyByEye()
    print '*INFO* entering fit at %s'%str(time.ctime())
    roi.fit(mysource=src)
    minNeg = None
    minPos = None
    llh_scan = None
    # now need minos
    if minos:
        minos_id = None
        freePars = {i:p.getName() for i,p in enumerate(roi.likelihood_fcn.params()) if p.isFree()}
        for i,p in enumerate(sorted(freePars)):
            if freePars[p]==par:
                minos_id = i
        if minos_id is None:
            raise RuntimeError("could not find MINOS ID")
        print '*INFO* minos ID: ',minos_id
        try:
            minNeg,minPos = roi.minuit_object.Minos(minos_id,cl)
        except RuntimeError:
            print '*ERROR* could not determine MINOS, but we bravely carry on!'
        src.expand({"minos":(minPos,minNeg)})
    if scan:
        sc_min = scan_min
        sc_max = scan_max

        if not minNeg is None:
            sc_min = minNeg

        if not minPos is None:
            sc_max = minPos
        # do the scan
        llh_scan = roi.likelihood_fcn.scan(src.name,par, xmin=sc_min, xmax=sc_max, npts=scan_npts,tol=1e-10, optimizer=configuration.optimizer, optObject=roi.minuit_object)

    # last not least, need store stuff
    #roi.print_summary(src)
    # for now, store dict...
    out = {}
    out["fitResultXml"]=roi.exportFitResultToDict() # store stuff as dict instead of xml!
    out['model']="CR"
    Id = roi.likelihood_fcn.par_index(src.name,par)
    out[par]={'mle':roi.likelihood_fcn[Id].parameter.getValue(),
                   'scale':roi.likelihood_fcn[Id].parameter.getScale(),
                   'Ts':roi.likelihood_fcn.Ts(src.name,reoptimize=True),
                   'Npred':roi.likelihood_fcn.logLike.NpredValue(src.name)}

    if minos:
        out[par].update({'pos':minPos,'neg':minNeg})
    if scan:
        out[par].update({'scan':llh_scan})
    if std_diffuse:
        print '*running with Standard Diffuse model*'
        Id = roi.likelihood_fcn.par_index("GAL","Prefactor")
        out['GAL']={'Prefactor':roi.likelihood_fcn[Id].parameter.getValue()*roi.likelihood_fcn[Id].parameter.getScale(),
                    'Ts':roi.likelihood_fcn.Ts("GAL",reoptimize=True),
                    'Npred':roi.likelihood_fcn.logLike.NpredValue("GAL")}
        Id = roi.likelihood_fcn.par_index("GAL","Index")
        out['GAL']["Index"]=roi.likelihood_fcn[Id].parameter.getValue()*roi.likelihood_fcn[Id].parameter.getScale()
    
    Id = roi.likelihood_fcn.par_index("EGAL","Normalization")
    out['EGAL']={'Normalization':roi.likelihood_fcn[Id].parameter.getValue()*roi.likelihood_fcn[Id].parameter.getScale(),
                'Ts':roi.likelihood_fcn.Ts("EGAL",reoptimize=True),
                'Npred':roi.likelihood_fcn.logLike.NpredValue("EGAL")}
    out['CL']=cl
    out['llh']=roi.likelihood_fcn()
    if not llh_scan is None:
        #if len(llh_scan)==2:
        #    llh_scan[1]*=out['llh']
        out['llh_scan']=llh_scan
    out['FinalState']=finalstate
    out['ROI']=roi.name
    out['STOOLS']="ST-%s"%os.getenv("INST_DIR").split("/")[-1]
    out["llh0"]=roi.fitNull()                                               
    print '*INFO* done with likelihood at %s'%str(time.ctime())
    sleeptime = 15
    print '*INFO* sleeping for %is to avoid stressing disk'%int(sleeptime)
    time.sleep(sleeptime)    
    f = open(configuration.resultsfile,'rb')
    dumpstring = f.read()
    d = pickle.loads(dumpstring)
    d.update({"CR":out}) # update old dict
    # could also try to write yaml
    print '*INFO* sleeping for %is to avoid stressing disk'%int(sleeptime)
    time.sleep(sleeptime)    
    if yamlfile is None:
        yamlfile = configuration.resultsfile.replace(".pkl",".yaml")
    dout = {}
    if os.path.isfile(yamlfile):
        print '*INFO* results yaml file exists already. Appending results - use force=True to start anew.'
        if update_yaml:
            dout = yaml.load(open(yamlfile,'rb'))
            print '*INFO* found following datapoints in yamlfile: {}'.format(dout.keys())
    dout.update(d)
    yaml.dump(dout,open(yamlfile,"wb"))
    # now write back
    d_pick = pickle.dumps(d,-1)    
    fo = open(configuration.resultsfile.replace(".pkl",".out"),'w')
    fo.write(str(d))
    fo.close()
    f = open(configuration.resultsfile,"wb")
    f.write(d_pick)
    f.close()
    del d_pick
    print '*INFO* done at %s'%str(time.ctime())

def load(chunk):
    # returns ROI, configuration
    print '*INFO* reading pickle from %s'%chunk
    f = open(chunk,'rb')
    dumpstring = f.read()
    inputblock = pickle.loads(dumpstring)
    ROI = inputblock['ROI']
    configuration = inputblock['configuration']
    ROI.minos_id = inputblock['minos_id']
    ROI.set_configuration(configuration)
    return ROI,configuration

def pack(ROI,configuration,minos_id):
    # this code is used to pack!
    tempfile = NamedTemporaryFile(prefix=os.getenv('MYTMPDIR')+'/')
    dumpfile = tempfile.name
    tempfile.close()
    dump_dict = {'ROI':ROI,'configuration':configuration,'minos_id':minos_id}
    dumpstring = pickle.dumps(dump_dict,-1) # use highest protocol
    fdumpfile = open(dumpfile,'wb')
    fdumpfile.write(dumpstring)
    fdumpfile.close()
    return dumpfile

    


def likelihood(roi,mass_points,dry=True,ignore_batch=False,
               final_states=None,j=1.3e18,scan=False,scan_min=None,
               scan_max=None,scan_npts=20,debug=False,minos_id=19,IC=False,
               model="PS",sleep='1m',lsf_queue="bullet-xml",force_srcmap=True,
               std_diffuse=True,CR=False,update_yaml=False):
    # scan_max = 200
    # that's the drive routine
    print '*INFO* entering likelihood at %s'%str(time.ctime())
    if roi.__dict__.has_key("egal"):
        os.environ["EGAL"]=roi.egal
    if roi.__dict__.has_key("gal"):
        os.environ["GAL"]=roi.gal
    print 'GAL: %s\nEGAL: %s'%(os.getenv("GAL"),os.getenv("EGAL"))
    # pickle a dict
    fstates = ['bbar','mumu','ww','tautau']
    if not final_states is None:
        fstates = final_states
    batch_jobs = []
    if not CR:
        for f in fstates:
            tmp_masses = []
            masses = mass_points
            if IC:
                if f == "mumu": masses = mumuIC.keys()
                elif f == "ee": masses = eeIC.keys()
                #else:
                #    print '*INFO*: ignoring final state %s'%f
                #    continue
            if f == "ww":
                for m in mass_points:
                    if m>=80:
                        tmp_masses.append(m)
                masses = tmp_masses
            print '*INFO*: mass points for %s %s'%(f,str(masses))
            roi.configuration.FinalState=f
            roi.configuration.resultsfile = os.path.join(roi.configuration.my_run_dir,"outputLLH_%s_%s.pkl"%(roi.name,f))
            d = {}
            d_pick = pickle.dumps(d,-1)    
            f = open(roi.configuration.resultsfile,"wb")
            f.write(d_pick)
            f.close()
            del d_pick
            if not IC: chunk = pack(roi,roi.configuration,minos_id)
            #print '*INFO* packing chunk %s'%chunk
            cmds = []
            for m in masses:
                if IC: chunk = pack(roi,roi.configuration,minos_id)
                cmd = 'python fermi-virgo/base/virgo_analysis.py "runLikelihood" %s %1.8f %1.8e'%(chunk,m,j)
                if IC:                   cmd+=" --IC" # activate IC
                if scan:                 cmd+=" --scan"
                if not scan_min is None: cmd+=" --scan_min=%1.4e"%float(scan_min)
                if not scan_max is None: cmd+=" --scan_max=%1.4e"%float(scan_max)
                if not scan_npts is None:cmd+=" --scan_npts=%i"%int(scan_npts)
                if not std_diffuse: cmd+=" --no-std-diffuse"
                if update_yaml: cmd+=" --update-yaml"
                cmds.append(cmd)
            # build the job-array:
            jobname = "Likelihood_%s_%s_%s"%(model,roi.name,roi.configuration.FinalState)
            lsf.bsub(jobname, cmds, logfiles=None,submit=not(dry), sleep=sleep,q=lsf_queue)#"bullet-xxl")
    else:
        cmds = []
        roi.configuration.FinalState= "CR"
        roi.configuration.resultsfile = os.path.join(roi.configuration.my_run_dir,"outputLLH_%s_%s.pkl"%(roi.name,"CR"))
        d = {}
        d_pick = pickle.dumps(d,-1)    
        f = open(roi.configuration.resultsfile,"wb")
        f.write(d_pick)
        f.close()
        del d_pick
        chunk = pack(roi,roi.configuration,minos_id)
        print '*INFO* packing chunk %s'%chunk
        cmd = 'python fermi-virgo/base/virgo_analysis.py "runLikelihoodCR" %s %1.8f %1.8e'%(chunk,100.,j)
        if scan:                 cmd+=" --scan"
        if not scan_min is None: cmd+=" --scan_min=%1.4e"%float(scan_min)
        if not scan_max is None: cmd+=" --scan_max=%1.4e"%float(scan_max)
        if not scan_npts is None:cmd+=" --scan_npts=%i"%int(scan_npts)
        if not std_diffuse: cmd+=" --no-std-diffuse"
        if update_yaml: cmd+=" --update-yaml"

        print '*INFO* cmd to add to array: %s'%cmd
        cmds.append(cmd)
        # build the job-array:
        jobname = "Likelihood_%s_%s_%s"%(model,roi.name,roi.configuration.FinalState)
        lsf.bsub(jobname, cmds, logfiles=None,submit=not(dry), sleep=sleep,q=lsf_queue)#"bullet-xxl")
        
if __name__ == "__main__":
    # that's what is invoked once code is called... from outside', called with 1 attribute *chunk to unpack*
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    parser.add_option("--yaml",dest='yaml_out', type=str, default = None, help="")
    parser.add_option("--scan_min",dest='scan_min', type=float, default = 0, help="")
    parser.add_option("--scan_max",dest='scan_max', type=float, default = None, help="")
    parser.add_option("--scan_npts",dest='scan_npts', type=int, default = 20, help="")
    #parser.add_option("--minosID",dest='minosID', type=int, default = 19, help="")
    parser.add_option("--scan",dest='scan', action='store_true', default = False,
                      help = "do scan")
    parser.add_option("--update-yaml",dest='update_yaml', action='store_true', default = False,
                      help = "update already existing yamls")
    parser.add_option("--no-std-diffuse",dest='diffuse', action='store_false', default = True,
                      help = "use non-std diffuse model, e.g. split models")
    parser.add_option("--IC",dest='ic',action="store_true",default=False,help="enable IC (needs pre-determined masses etc.)")
    (opts, args) = parser.parse_args()
    args = sys.argv
    if len(args)<2:
        raise RuntimeError("Must call with mode")
    if args[1] == "runLikelihood":
        if len(args)<4:
            raise RuntimeError("Must call this code with a chunk argument followed by mass point")
        else:
            chunk = args[2]
            if not toolbox.is_file(chunk):
                raise RuntimeError("cannot process chunk, not a valid file")
            ROI,configuration=load(chunk)
            tempDir = configuration.getTempDir() # builds the temp dir!
            setPF(tempDir) # re-route pfiles
            make_srcmap = False
            if opts.ic:
                make_srcmap = True
            process_likelihood(ROI,configuration,float(args[3]),j=float(args[4]),scan=opts.scan,
                               scan_min=opts.scan_min,scan_max=opts.scan_max,scan_npts=opts.scan_npts,
                               force_srcmap=make_srcmap,std_diffuse=opts.diffuse,
                               yamlfile=opts.yaml_out, update_yam=opts.update_yaml)
            print '*INFO* done with running, attempting to remove chunk %s'%chunk
            configuration.removeTempDir() # remove the temp dir!
            os.remove(chunk)
            resetPF()
    elif args[1] == "runLikelihoodCR":
        if len(args)<4:
            raise RuntimeError("Must call this code with a chunk argument followed by mass point")
        else:
            chunk = args[2]
            if not toolbox.is_file(chunk):
                raise RuntimeError("cannot process chunk, not a valid file")
            ROI,configuration=load(chunk)
            tempDir = configuration.getTempDir() # builds the temp dir!
            setPF(tempDir) # re-route pfiles
            make_srcmap = False
            if opts.ic:
                make_srcmap = True
            print '*Use Std Diffuse: %s*'%opts.diffuse
            process_likelihoodCR(ROI,configuration,model="file",j=float(args[4]),scan=opts.scan,scan_min=opts.scan_min,
                                 scan_max=opts.scan_max,scan_npts=opts.scan_npts,force_srcmap=make_srcmap,
                                 std_diffuse=opts.diffuse,yamlfile=opts.yaml_out,update_yaml=opts.update_yaml)
            print '*INFO* done with running, attempting to remove chunk %s'%chunk
            configuration.removeTempDir() # remove the temp dir!
            os.remove(chunk)
            resetPF()            
    else:
        print '*do nothing*'
