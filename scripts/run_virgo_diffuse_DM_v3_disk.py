import base.virgo_analysis as VA
import sys, os, time
from tempfile import NamedTemporaryFile
from base.xmltools import expandEnvVarsInXml
from base.common import *
  
from optparse import OptionParser
usage = "Usage: %prog  [options] input"
description = "python script"
parser = OptionParser(usage=usage,description=description)
parser.add_option("--delay",dest='delay', type=int, default = 10, help="delay in seconds between submissions from this script")
parser.add_option("--dry",dest='dry', action='store_true', default = False,
                  help = "using this flag will create a bash script run_env.sh that contains the environment vars defined at runtime")
parser.add_option("--debug",dest='debug', action='store_true', default = False,
                  help = "")
parser.add_option("--createTemplate",dest='makeTemplate', action='store_true', default = False,
                  help = "if used, re-create template files")
parser.add_option("--nullfit",dest='null', action='store_true', default = False,
                 help = "using this flag to run nullhypothesis")
parser.add_option("--expand-vars",dest='evars', action='store_true', default = False,
                 help = "using this flag to run with temporary xml files to expand all env vars at runtime")
parser.add_option("--split",dest='split', action='store_true', default = False,
                 help = "only do split models")
parser.add_option("--config-overwrite",dest='cfg', type='str', default = None,
                 help = "config override of the form name=value, separate by using ;")
parser.add_option("--sleep",dest='sleep', type=float, default = 15,
                 help = "sleep period in seconds")
parser.add_option("--setenv",dest='setenv', type='str', default = None,
                 help = "set env vars, separate by using ;")
parser.add_option("--queue",dest='queue', type='str', default = "bullet-xxl",
                 help = "LSF queue")
parser.add_option("--ignore-logs",dest='ignore_logs', action='store_true', default = False,
                 help = "use this to not change log-dir content (expert option)")
parser.add_option("--fstate",dest="fstate",default=None,help="if set, specify the final state")
parser.add_option("--mass",dest="mass",default=None,help="if set, specify the masses")
parser.add_option("--model",dest="model",default=None,help="instead of running zillion models, choose which one to pick, std is Std_Gal_diffuse")
parser.add_option("--target",dest="target",default=None,help="instead of running over all targets: nfw, gao and miguel")
parser.add_option("--update-yaml",dest='update_yaml', action='store_true', default = False,help = "update already existing yamls")

(opts, args) = parser.parse_args()

step = sys.argv[1]

labels = ["nfw","miguel","gao"]
jvals = [4.41000000e+17,   1.18600000e+19,   4.08100000e+20]
targets = dict(zip(labels,jvals))
lorimerModels = ["Lorimer_z10_Ts100000","Lorimer_z4_Ts100000","Lorimer_z10_Ts150","Lorimer_z4_Ts150"]
snrModels = ["SNR_z10_Ts100000","SNR_z4_Ts100000","SNR_z10_Ts150","SNR_z4_Ts150"]
templates = []

if not opts.target is None:
    opts.target = opts.target.split(",")
    targets = {t for t in targets if t in opts.target}

for i,jlabel in enumerate(targets.keys()):
    for j,variant in enumerate(["fixed","split"]):
        psplit = False
        if variant == "split": psplit = True
        label = "%s.%s"%(jlabel,variant)
        modelFile = "/nfs/farm/g/glast/u55/zimmer/FermiData/Clusters/Virgo_v3/DoubleDisk/%s/Virgo_%s.xml"%(variant,jlabel.upper())
        extFits = "/nfs/farm/g/glast/u55/zimmer/FermiData/Clusters/Virgo_v3/Templates/coadd.%s.fits"%jlabel.upper()
        models = lorimerModels+snrModels
        if variant == "fixed":
            models+=["Std_Gal_diffuse"]
        for model in models:
            diffTag = None
            if model in lorimerModels+snrModels: diffTag = model
            templates.append(VirgoContainer(name="Virgo_%s.%s"%(model,label),J=targets[jlabel],fits=extFits,Xml=modelFile,diffuseTag=diffTag,split=psplit))

_templates = templates
if not opts.model is None:
    opts.model = opts.model.split(",")
    templates = [t for t in templates if t.name in opts.model]
    print 'working on models: {}'.format(templates)

for t in templates:
    if not t.diffuseTag is None:
        if opts.evars: t.dumpXml()

if opts.split: # test ONLY for split models
    print "test for split models"
    templates = [t for t in templates if t.split]
    print 'found %i split analyses'%len(templates)


    

rdir = "/afs/slac/g/glast/users/zimmer/FermiData/Clusters/Virgo_Diffuse_v3_DD/"
if not os.path.isdir(rdir): os.system("mkdir -p %s"%rdir)
rois = VA._setup(templates,rundir=rdir)
# uncomment
final_states = ['bbar','tautau','ww','mumu','ee']
jobs = []
masses = [5,10,20,50,100,200,500,1000.,2000.]

if not opts.fstate is None:
    if ',' in opts.fstate:
        final_states = opts.fstate.split(",")
    else:
        final_states = [opts.fstate]
if not opts.mass is None:
    if ',' in opts.mass:
        masses = [float(m) for m in opts.fstate.split(",")]
    else:
        masses = [float(opts.mass)]

if opts.debug:
    final_states = ['bbar']
    masses = [20]
    rois = [rois[-1]]


do_roi = True

if step == "PREPARE":
    for i,roi in enumerate(rois):
        jobs += [VA.make_prepare(roi,create_virgo=opts.makeTemplate,dry=opts.dry,template_dict=templates[i])]
        time.sleep(opts.sleep)
    print jobs

elif step == "LIKELIHOOD":
    for i,roi in enumerate(rois):
        std_diffuse=True
        if templates[i].split: 
            print '*working on %s'%templates[i].name
            std_diffuse=False
        VA.likelihood(roi,masses,dry=opts.dry,ignore_batch=True,final_states=final_states,
                      j=templates[i].getJ(),scan=False,scan_npts=15,model="DD",
                      lsf_queue=opts.queue,std_diffuse=std_diffuse,update_yaml=opts.update_yaml)#,debug
        time.sleep(opts.sleep)
else:
    raise Exception("Not supported!")
