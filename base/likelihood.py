# ScienceTools modules
#  note the order!
importError = False
try:
    from GtApp import GtApp
    from BinnedAnalysis import *
    from Composite2 import Composite2
    from tempfile import NamedTemporaryFile
    import LikelihoodState
except ImportError:
    print('*ERROR* could not import ScienceTools modules')
    importError = True
# standard python modules
import os, time, sys, math, copy
#sys.argv.append('-b') # forces to load root in batch mode
import ROOT as R

if not sys.flags.interactive:
    R.gROOT.SetBatch(True)
import pyfits, yaml
import xml.parsers.expat
import xml.dom.minidom as xdom
import StringIO
import os,gc, traceback, glob, time, random, shutil
import pickle
# numpy
import numpy as np
# matplotlib
import matplotlib

matplotlib.use('Agg') # if we don't have X11, use Agg as backend
import matplotlib.pyplot as plt
# my code
import toolbox, xmltools
bsub_opts = "-R \"rhel60 && linux64\""

LSD = 6
global_queues = {'express': 0, 'short': 1, 'medium': 2, 'long': 3, 'xlong': 4, 'xxl': 5}

def makeXmlElementFromDict(d,tagname="parameter"):
    """
    given a dictionary 'd' create an xml string.
    :param d:
    :return:
    """
    pp = xdom.parseString('<?xml version="1.0" ?><parameter_list title="jobTitle"/>')
    par = pp.createElement(tagname)
    if not isinstance(d,dict):
        raise Exception("can only be called with dict object")
    for key in d.keys():
        value = d[key]
        if type(value)==str:
            par.setAttribute(key,value)
        elif type(value)==str or type(value)==float:
            par.setAttribute(key,"%1.5e"%float(value))
        else:
            pass # do nothing.
    return par

def queue(i):
    for key in global_queues:
        if global_queues[key] == i:
            return key


def SearchForScaleFactors(xmlfile, srcname):
    print('*INFO* search for scalefactors in fitsfile')
    print('xmlfile %s source %s' % (xmlfile, srcname))
    # let's extract the FITS location
    scalefactor = 1.0
    fil = None
    extended = False
    xmlf = xdom.parse(xmlfile)
    xmlsrc = {}
    for s in xmlf.getElementsByTagName("source"):
        sname = s.getAttribute("name")
        xmlsrc[sname] = s
    if not srcname in xmlsrc.keys():
        raise Exception("*FATAL* cannot find source ", srcname, " in xmlfile, abort")
    else:
        # next check for spatialModel
        theSource = xmlsrc[srcname]
        if theSource.getAttribute("type") == "DiffuseSource":
            spatialModel = theSource.getElementsByTagName("spatialModel")
            if not len(spatialModel) > 0:
                raise Exception("*FATAL keyword extended in ROI file found, but xmlfile does not contain spatialModel, check both files")
            else:
                fil = spatialModel[0].getAttribute("file")
                if not toolbox.ContainsExpr(fil, ".fits"):
                    raise Exception("*FATAL no FITS file found")
                print('       extracting info from %s'%fil)
                hdu = pyfits.open(fil)[0]
                if hdu.header.has_key("sfactor"):
                    scalefactor = float(hdu.header['sfactor'])
                    extended = True
                else:
                    print('*WARNING* keyword extended in ROI file found, but scalefactor not present in FITS file, check header and ROI file')
    print('scalefactor %e' % scalefactor)
    return scalefactor, fil, extended


def format_output_str(inp):
    s = '%.3e' % float(inp)
    return s


def convBool(string):
    if string == 'True':
        return True
    else:
        return False


def conv2Days(TIME):
    # to convert seconds to days
    return TIME / (3600 * 24.)


def prettifyXML(xmlfile):
    # to fix a bug in xdom.toprettyxml()
    # loads the xml file and removes blank lines
    new_lines = []
    f = open(xmlfile, 'r')
    lines = f.readlines()
    for line in lines:
        # test this line:
        # let's try to reduce the line length...
        c_string = copy.deepcopy(line)
        c_string.replace("\n", "")
        for i in range(0, 10):
            new_str = c_string.partition("\t")
            if len(new_str[1]) > 0:
                c_string = new_str[2]
                #print 'step i',i,' string length',len(c_string)
        if len(c_string) > 1:
            new_lines.append(line)
    f.close()
    f = open(xmlfile, 'w')
    for l in new_lines:
        f.write(l)
    f.close()


class triple(dict):
    # a helper class, done!
    def __init__(self, name, val, plus, minus):
        self.name = name
        self.tuple = (float(minus), float(val), float(plus))
        self.value = float(val)
        self.limits = (-float(minus) + float(val), float(val) + float(plus))
        self.pos = float(plus)
        self.neg = float(minus)
        self.min = round(self.neg, 2)
        self.max = round(self.pos, 2)
        self.scale = 1.0
        self.error = 0.
        self.setboundsmanual = False
        super(triple,self).update(self.__dict__)
    
    def get(self,key):
        if key in self.__dict__:
            return self.__dict__[key]
    
    def __getstate__(self):
        # allow pickling
        d = copy.copy(self.__dict__)
        return d

    def __setstate__(self, state):
        self.__dict__ = state

    def get_limits(self):
        return self.limits[0] * self.scale, self.limits[1] * self.scale

    def get_tuple(self):
        return self.tuple

    def set_minmax(self, mmin, mmax):
        #print '*DEBUG* set_minmax()'
        self.setboundsmanual = True
        self.min = float(mmin)
        self.max = float(mmax)

    def make_bounds(self, factor=2.):
        if self.setboundsmanual:
            print('*INFO* parameter %s has already associated min/max values' % self.name)
            pass
        else:
            #print 'making bounds for ',self.name
            #print 'old values ',self.min,self.max
            #MIN = float(1./factor)*float(self.value)
            MAX = float(1. * factor) * float(self.value)
            self.min = 1e-2 # capping too low...
            self.max = MAX
            #print 'new values: ',self.min,self.max

    def setValue(self, value):
        self.value = value

    def setScale(self, scale):
        self.scale = scale


class my_value(object):
    def __init__(self, name, value, **kwargs):
        self.name = name
        self.scale = 1.0
        self.error = 0.
        self.value = float(value)
        self.__dict__.update(kwargs)

    def get(self,key):
        if key in self.__dict__:
            return self.__dict__[key]

    def __getstate__(self):
        # allow pickling
        d = copy.copy(self.__dict__)
        return d

    def __setstate__(self, state):
        self.__dict__ = state

    def setValue(self, value):
        self.value = value

    def setScale(self, scale):
        self.scale = scale


class source(object):
    # decorator class.../ have to think about it.
    def __init__(self, **kwargs):
        self.Type = 'Generic'
        self.name = None
        self.parameter = None
        self.scalefactorfromFITS = 1.0
        self.TSvalue = None  #each source should have a TS value...
        self.JvalName = None
        self.is_tied = "1"
        self.is_free = "1"
        self.__dict__.update(kwargs)
        self.xmlParameters = {}
    
    def getXmlParameters(self):
        return self.xmlParameters
    
    def __getstate__(self):
        # allow pickling
        d = copy.copy(self.__dict__)
        return d

    def __setstate__(self, state):
        self.__dict__ = state

    def __repr__(self):
        return "%s\t%s" % (self.name, str(self.__dict__))

    def has_key(self, key):
        if self.__dict__.has_key(key):
            return True
        else:
            return False

    def makeProfilePlot(self, fitsfile=None, output=None):
        from toolbox import GetCoordinatesFromFitsFile
        # new method to plot a profile
        # stolen from the ExtendedTemplateMaker
        ###
        if fitsfile is None or output is None:
            print('*INFO* makeProfilePlot needs to be called with input and output file')
            return 0
        else:
            print('*INFO* create plot of extended profile')
            gc.enable()
            im1 = plt.figure(figsize=(6, 6), dpi=80)
            matplotlib.rcParams['font.family'] = 'serif'
            #matplotlib.rcParams['verbose.level'] = 'silent' # to silence the ANNOYING outputs....
            #matplotlib.rcParams['verbose.fileo'] = 'matplotlib_warning.log'
            im1.set_label("Control Plot")
            #import plt as plt
            gc.enable()
            #################################
            # Now to do a little sanity check
            # and see if all the profiles are consistent
            # And the fits file
            fits = pyfits.open(fitsfile)
            header = fits[0].header
            data = fits[0].data

            npix, x0 = header['NAXIS1'], header['CRPIX1']
            isodd = npix % 2 # Deal with even and odd numbers of pixeheader['CDELT1'])
            delta = abs(header['CDELT1'])
            x_fits = np.arange(0.5 * (not isodd), x0) * delta
            pdf_fits = data[0, x0 - 1 * isodd:, x0 - 1 * isodd]
            RA, DEC = GetCoordinatesFromFitsFile(fitsfile)
            # Make a plot...
            identifier = 'fitsfile: %s \nCF: %2.4e RA: %f DEC: %f' % (fitsfile, self.scalefactorfromFITS, float(RA), float(DEC))
            plt.figtext(0.5, 0.935,identifier, ha='center', color='black')#, weight='bold', size='large')
            ax = im1.add_subplot(111)
            ax.plot(x_fits, pdf_fits, '-k', label="Fits Profile")
            ax.set_yscale('log')
            ax.legend(ncol=1, loc='best')
            ax.set_ylabel('intensity [1/sr]')
            ax.set_xlabel('r [deg]')
            # PLOT EVERYTHING
            im1.canvas.draw()
            cont_name = output + '/' + self.name + '.extendedProfile.png'
            im1.savefig(cont_name)
            plt.close()
            print('*INFO* file %s written' % cont_name)
            del im1
            del ax
            # keep pickle to allow easy access for later plots.
            d = {"radius [deg]":np.array(x_fits),"intensity [1/sr]":np.array(pdf_fits),"meta":identifier}
            f = open(cont_name.replace(".png",".pkl"),"wb")
            f.write(pickle.dumps(d,-1))
            f.close()
            gc.collect() # invoke the garbage collector                 


    def expand(self, my_dict):
        self.__dict__.update(my_dict)


    def make_scalefactor(self):
        Min = None
        Max = None
        if 'ScaleFactor_min' in self.__dict__.keys():
            #print '*found Scalefactor min*'
            Min = float(self.ScaleFactor_min)
        if 'ScaleFactor_max' in self.__dict__.keys():
            #print '*found Scalefactor max*'
            Max = float(self.ScaleFactor_max)
        if 'ScaleFactor' in self.__dict__.keys():
            scalefactor = triple('ScaleFactor', float(self.ScaleFactor), 0., 1.)
        else:
            scalefactor = triple('ScaleFactor', 0.5, 0., 1.)
            scalefactor.scale = 1.0
        if 'ScaleFactor_scale' in self.__dict__.keys():
            scalefactor.scale = float(self.ScaleFactor_scale)
        if not (Min is None and Max is None):
            scalefactor.set_minmax(Min, Max)
        elif Min is None and not Max is None:
            scalefactor.set_minmax(1e-8, Max)
        elif Max is None and not Min is None:
            scalefactor.set_minmax(Min, 1e8)
        if 'ScaleFactor_value' in self.__dict__.keys():
            scalefactor.value = float(self.ScaleFactor_value)
        self.xmlParameters["ScaleFactor"]=scalefactor
        return scalefactor

    def GenericSource(self):
        if self.has_key('xmlpars'):
            if len(self.xmlpars) > 0:
                raise Exception("** Source association duplicated, self.xmlpars exists already: ", self.xmlpars, "**")
        else:
            theParameters = {}
            d = {'xmlpars': theParameters}
            self.expand(d)
            self.Type = 'Generic'
        return

    # for triple values we use the triple class, otherwise my_value
    def DMFitSource(self, ScaleFactor=False):
        # that associates a number of parameters with the 
        # source in form of a dictionary containing the xml parameters
        if self.has_key('xmlpars'):
            if len(self.xmlpars) > 0:
                raise Exception("** Source association duplicated, self.xmlpars exists already: ", self.xmlpars, "**")
        else:
            mass = triple('mass', 150.0, 0, 2500) # set to 150 GeV
            #print mass
            #toolbox.stop_here()
            # not existent - nice.
            if self.has_key('bratio'):
                bratio = my_value('bratio', float(self.bratio))
            else:
                bratio = my_value('bratio', 1.0)
            if self.has_key('norm'):
                norm = triple('norm', float(self.norm), 0, 1e3)
                if self.has_key('norm_scale'):
                    norm.scale = float(self.norm_scale)
                if self.has_key('norm_error'):
                    norm.error = float(self.norm_error)
            else:
                norm = triple('norm', 1.0, 0, 1e3)
                norm.scale = 1e17
            if self.has_key('channel0'):
                channel0 = my_value("channel0", round(float(self.channel0), 0))
            else:
                channel0 = my_value('channel0', 4.0)
            if self.has_key('channel1'):
                channel1 = my_value('channel1', int(self.channel1))
            else:
                channel1 = my_value('channel1', 1.0)
                # the rest we don't want to touch
            norm.value *= self.scalefactorfromFITS
            norm.min = 1e-2
            norm.max = 1e2
            #norm.make_bounds()
            # use a triple instead
            if self.has_key('mass'):
                mass = triple('mass', float(self.mass), 0, 2500)
            if self.has_key('mass_scale'):
                mass.setScale(float(self.mass_scale))
            mass.set_minmax(0, 2500)
            sigmav = triple('sigmav', 1.1, 1.5e12, 1.5e-12)
            if self.has_key('sigmav_scale'):
                sigmav.scale = float(self.sigmav_scale)
            else:
                sigmav.scale = 1e-25
            theParameters = {'mass': mass, 'bratio': bratio,
                             'sigmav': sigmav, 'norm': norm, 'channel0': channel0, 'channel1': channel1}
            self.Type = 'DMFit'
            if ScaleFactor:
                theParameters['ScaleFactor'] = self.make_scalefactor()
                self.Type += "::ScaleFactor"
            d = {'xmlpars': theParameters}
            self.expand(d)
            for key in theParameters:
                if not key in self.xmlParameters:
                    self.xmlParameters[key]=theParameters[key]
            return

    def FileSource(self, ScaleFactor=False):
        if self.has_key('xmlpars'):
            if len(self.xmlpars) > 0:
                raise Exception("** Source association duplicated, self.xmlpars exists already: ", self.xmlpars, "**")
        else:
            # not existent - nice.
            if self.has_key('Normalization'):
                normalization = triple('Normalization', float(self.Normalization), 0., 1.)
            else:
                normalization = triple('Normalization', 1., 0., 1.)
            if self.has_key('Normalization_scale'):
                normalization.scale = float(self.Normalization_scale)
            else:
                normalization.scale = 1.0
            normalization.value *= self.scalefactorfromFITS
            if self.has_key('Normalization_min') and self.has_key('Normalization_max'):
                normalization.set_minmax(self.Normalization_min, self.Normalization_max)
            theParameters = {'Normalization': normalization}

            self.Type = 'FileSpectrum'
            if ScaleFactor:
                theParameters['ScaleFactor'] = self.make_scalefactor()
                self.Type += "::ScaleFactor"
            d = {'xmlpars': theParameters}
            self.expand(d)
            for key in theParameters:
                if not key in self.xmlParameters:
                    self.xmlParameters[key]=theParameters[key]
            return

    def PowerLaw2(self, ScaleFactor=False):
        # a dirty fix for testing
        if self.has_key('xmlpars'):
            if len(self.xmlpars) > 0:
                raise Exception("** Source association duplicated, self.xmlpars exists already: ", self.xmlpars, "**")
        else:
            # not existent - nice.
            if self.has_key('Normalization'):
                normalization = triple('Integral', float(self.Normalization), 0., 1.)
            else:
                normalization = triple('Integral', 1., 0., 1.)
            if self.has_key('Normalization_scale'):
                normalization.scale = float(self.Normalization_scale)
            else:
                normalization.scale = 1.0
            normalization.value *= self.scalefactorfromFITS
            if self.has_key('Normalization_min') and self.has_key('Normalization_max'):
                normalization.set_minmax(self.Normalization_min, self.Normalization_max)
            theParameters = {'Integral': normalization}

            if self.has_key('Index'):
                # if it's in the keys, then self.Index is allowed...
                index = triple('Index', self.Index, -5, -0.1)
                index.set_minmax(-5, -0.1)
                if 'Index_min' in self.__dict__.keys() and 'Index_max' in self.__dict__.keys():
                    index.set_minmax(self.Index_min, self.Index_max)
                theParameters['Index'] = index
            self.Type = 'PowerLaw2'
            if ScaleFactor:
                theParameters["ScaleFactor"] = self.make_scalefactor()
                self.Type += "::ScaleFactor"
            d = {'xmlpars': theParameters}
            self.expand(d)
            for key in theParameters:
                if not key in self.xmlParameters:
                    self.xmlParameters[key]=theParameters[key]
            return


class configuration(object):
    # done!
    def __init__(self, xmlfile=None, cfgdefault="fermi-virgo/base/defaults_p7v6.xml"):#None):
        self.cfgdefault = cfgdefault
        self.loghome = os.getenv("LOGDIR")
        if self.loghome is None:
            raise Exception("$LOGDIR not defined, source setup.sh first!")
        self.config_xml = xmlfile
        self.LightCube = "None"
        self.Ind_DoFit = "True"
        self.C2_DoFit = "True"
        self.scan_likelihood = "False"
        self.seed = None
        self.C2_scan_likelihood = None
        self.C2_refit = "True"
        self.CL = 2.71
        self.C2_cummulative_fit = "False"
        self.C2_NumberOfTiedRois = None
        self.C2_freeze_sources = "False" # new expert option
        self.resultxml = None
        self.nullhypothesis = "False"
        self.mass = '0'
        self.mass_points = '12,20,50,100,500,1000,2000'
        self.whitness = None
        self.FixDiffuse = "True"
        self.event_class = "source" # hidden
        self.make_covariance = "False"
        self.make_residuals = "True"
        self.calculateSourceTs = "True" # expert option
        self.env_vars = {}
        self.load_defaults()
        self.chatter = 2
        self.tempdir = None
        self.set_env(copy.copy(self.env_vars))
        self.makeRandomSeed()
        self.seed = 0

    def makeRandomSeed(self):
        from random import randint
        return randint(1,1000) # make random seed
        
    def __getstate__(self):
        # allow pickling
        d = copy.copy(self.__dict__)
        return d

    def __setstate__(self, state):
        self.__dict__ = state


    def __call__(self):
        # if called directly, returns the old dict form used in other code
        pars = {}
        for k in self.__dict__.keys():
            pars[str(k)] = str(self.__dict__[k])
        return pars

    def set(self, key, value):
        self.__dict__.update({key: value})
    
    def get(self,key,callable=str):
        if key in self.__dict__:
            if type(self.__dict__[key])!=callable:
                print '*INFO* explicit type conversion from %s to %s'%(type(self.__dict__[key]),callable)
            return callable(self.__dict__[key])
    
    def readin_envvars(self):
        for key in self.env_vars:
            if self.chatter>2:
                print '*INFO* setting %s to %s'%(key,self.env_vars[key])
            os.environ[key] = self.env_vars[key]
            
    def load_defaults(self):
        # build standard vars
        print('*INFO* loading defaults from %s' % self.cfgdefault)

        xmlfile = xdom.parse(self.cfgdefault)
        os.environ['GAL'] = "/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/ring_2year_P76_v0.fits"
        os.environ['EGAL'] = "/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/isotrop_2year_P76_source_v0.txt"
        self.env_vars['GAL'] = os.getenv("GAL")
        self.env_vars['EGAL'] = os.getenv("EGAL")

        my_dict = {}
        for par in xmlfile.getElementsByTagName("parameter"):
            my_dict[str(par.getAttribute("name"))] = str(par.getAttribute("value"))
        self.__dict__.update(my_dict)
        if int(self.chatter) > 4:
            print(self.__dict__)
            #toolbox.stop_here()

    def read_xml(self):
        my_dict = {}
        xmlfile = xdom.parse(self.config_xml)
        # need to incorporate env block
        for var in xmlfile.getElementsByTagName("var"):
            varname = str(var.getAttribute("name"))
            value = str(var.firstChild.data)
            os.environ[varname] = value
            self.env_vars[varname] = value
        for par in xmlfile.getElementsByTagName("parameter"):
            my_dict[str(par.getAttribute("name"))] = str(par.getAttribute("value"))
        self.__dict__.update(my_dict)
        self.__read__()

    def read_yaml(self,yamlfile):
        print '*INFO* reading configuration as yaml file'
        my_dict = yaml.load(open(yamlfile).read())
        if 'var' in my_dict:
            for key in my_dict['var']:
                os.environ[key]=my_dict['var'][key] 
            self.env_vars.update(my_dict['var'])
        self.__dict__.update(my_dict)
        self.__read__()
        
    def __read__(self):
        try:
            self.prepdir = self.HOME + '/Prepared/'
            if len(self.PrepTag) > 0:
                self.prepdir += self.PrepTag + '/'
            else:
                self.prepdir += self.Tag + '/'
            self.outdir = self.HOME + '/results/' + self.Tag + '/'
            self.likeModelVersion = '_' + self.InclusionRadius + 'C_' + self.catalog
            if not os.getenv("INST_DIR") is None:
                self.stoolsversion = os.getenv('INST_DIR').partition("ScienceTools/")[2]
            else:
                self.stoolsversion = "N/A"
            lifetime_s = (float(self.TStop) - float(self.TStart))
            lifetime_d = conv2Days(lifetime_s)
            if int(lifetime_d) == 0:
                self.lifetime = str(lifetime_s) + 's'
            else:
                self.lifetime = str(int(lifetime_d)) + 'd'
                # if C2_scan_likelihood is not explicitly set, i.e. in config xml,
            # override None with the value of scan_likelihood
            if self.C2_scan_likelihood is None:
                self.C2_scan_likelihood = self.scan_likelihood
                #sys.exit()


        except KeyError:
            raise Exception('**ERROR: config file malformed, cannot parse directories')

    def write_xml(self, xmlfile):
        kstring = '<?xml version="1.0" ?>\n<parameter_list title="jobTitle"/>'
        xmlf = xdom.parse(StringIO.StringIO(kstring))
        env_block = xmlf.createElement("environmentVars")
        for key in self.env_vars:
            # ugly parsing... to avoid minidom's broken output
            var = xdom.parseString('<var name="%s">%s</var>' % (key, self.env_vars[key]))
            env_block.appendChild(var.firstChild)
        xmlf.lastChild.appendChild(env_block)
        for key in sorted(self.__dict__):
            if isinstance(self.__dict__[key], str):
                p = xmlf.createElement("parameter")
                p.setAttribute("name", key)
                p.setAttribute("value", self.__dict__[key])
                xmlf.lastChild.appendChild(p)
        f = open(xmlfile, 'w')
        #f.write(xmlf.toprettyxml())
        f.write(xmlf.toxml())
        # FIXME! toprettyxml creates bullshit xml text node (BUG!), normal xml does not have line breaks
        f.close()
        print('*INFO* wrote configuration file %s' % xmlfile)
    
    def write_yaml(self,yamlfile):
        yaml.dump(self(),open(yamlfile,'w'))
        print('*INFO* wrote configuration yaml file %s'%yamlfile)

    def set_tag(self, tag=None):
        # that should get me plots for each mass point
        if self.PhysicsMode.startswith("DM"):
            self.specialTag += self.mass + self.FinalState
        if not tag is None:
            self.specialTag += tag

    def set_env(self, mydict):
        for key in mydict.keys():
            os.environ[key] = mydict[key]

    def setResultsFile(self, resultsfile):
        self.resultxml = resultsfile
        if self.chatter > 4:
            print('*INFO* resultsfile %s' % os.path.abspath(self.resultxml))

    def expand(self, my_dict):
        # okay we have derived variables (see in read_xml) but also could add arbitrary stuff to the configuration (not written to XML)
        self.__dict__.update(my_dict)

    def update(self, my_dict):
        if my_dict.has_key("cfgefault"):
            self.cfgdefault = str(my_dict['cfgdefault'])
            self.load_defaults()
        self.__dict__.update(my_dict)
        if convBool(self.nullhypothesis):
            print('*** TESTING FOR NULLHYPOTHESIS ***')
            self.likeModelVersion += '_null'

    def MakeNewResultsXML(self):
        individual = convBool(self.Ind_DoFit)
        combined = convBool(self.C2_DoFit)
        if combined and individual:
            return True
        if combined:
            if not individual:
                return False
            else:
                return True
        if individual:
            if not combined:
                return True
        if not combined and not individual:
            print(self.__dict__)
            raise Exception("Neither individual nor combined fit requested, giving up")

    def CheckXML_Integrity(self, xmlfile, rois):
        self.setResultsFile(xmlfile)
        # this one performs a bunch of tests
        xmlf = None
        if not os.path.exists(self.resultxml):
            return False
        else:
            if not os.path.getsize(self.resultxml) > 0:
                print('*INFO* file %s exists but size zero' % self.resultxml)
                return False
            else:
                # now need to parse it
                print('*INFO* file %s exists and is nonzero, perform integrity check' % self.resultxml)
                try:
                    xmlf = xdom.parse(self.resultxml)
                except xml.parsers.expat.ExpatError:
                    print('*WARNING* parsing failed, assume corrupt file')
                    return False
                if xmlf is None:
                    return False
                else:
                    if not convBool(self.Ind_DoFit) and convBool(self.C2_DoFit):
                        # now perform deep check
                        rois_from_xml = xmlf.getElementsByTagName("ROI")
                        if len(rois_from_xml) != len(rois):
                            return False
                        else:
                            # now check if names are the same
                            roi_names_from_xml = []
                            roi_names = []
                            for r in rois: roi_names.append(str(r.name))
                            for r in rois_from_xml: roi_names_from_xml.append(str(r.getAttribute("name")))
                            # sort it 
                            sorted_roi_xml = sorted(roi_names_from_xml)
                            sorted_rois = sorted(roi_names)
                            if sorted_roi_xml == sorted_rois:
                                cl = xmlf.getElementsByTagName("CompLike")
                                if len(cl) > 0:
                                    # now need to check if complike object exists
                                    if len(cl.getElementsByTagName("parameter")) > 0:
                                        return False
                                    else:
                                        print('*INFO* stored roi information fits expectation, appending content')
                                        return True
                                else:
                                    print('*INFO* stored roi information fits expectation, appending content')
                                    return True
                            else:
                                return False

    def buildResultsFile(self, resultxml=None):
        if resultxml is None:
            resultxml = self.resultxml
            # now if we have NOT IndividualFit but Combined AND file exists, then just return 0
        if os.path.isfile(resultxml):
            print("*INFO* result xml %s exists already, overwriting content."%resultxml)
            return
        # builds the basic structure of the resultsfile
        mass = self.mass
        mode = self.PhysicsMode
        first_line = "<?xml version=\"1.0\" ?>\n"
        xmlstr = first_line
        #  file's not existent, so let's build one

        irf = str(self.IRF)
        path = str(self.outdir)
        runtime = str(time.ctime())

        second_line = "<result runtime=\"%s\" dataperiod=\"%s\" " \
                      "sciencetools=\"%s\" irf=\"%s\" path=\"%s\" " \
                      "physicsmode=\"%s\" mass=\"%s\">\n" % (str(runtime), str(self.lifetime),
                                                             str(self.stoolsversion), str(irf),
                                                             str(path), str(mode), str(mass))
        last_line = "</result>"
        f = open(resultxml, 'w')
        f.write(first_line)
        f.write(second_line)
        f.write(last_line)
        f.close()
        xmlfile = xdom.parse(resultxml)
        for r in xmlfile.getElementsByTagName("result"):
            r.setAttribute("GAL", str(os.getenv('GAL')))
            r.setAttribute("EGAL", str(os.getenv('EGAL')))
            cobj = xmlfile.createElement("CompLike")
            iobj = xmlfile.createElement("IndividualFit")
            r.appendChild(iobj)
            r.appendChild(cobj)
            xmlstr += r.toxml()

        f = open(resultxml, 'w')
        f.write(xmlstr)
        f.close()
        # done

    def getTempDir(self):
        base = "/tmp"
        BATCH = toolbox.isBatch()
        if BATCH:
            base = "/scratch"
        job_id = str(np.random.randint(1,999999999))
        usr = os.getenv("USER","zimmer")
        o = os.path.join(base,usr,job_id)
        print '*INFO* using temp path: %s'%o
        toolbox.mkdir(o)
        self.tempdir = o
        return o

    def removeTempDir(self):
        o = self.tempdir
        print '*INFO* removing %s'%o
        toolbox.rmdir(o)
        return 

if importError:
    pass
else:

    class CompositeLikelihood2(Composite2):
        # my reminder for super:
        # from python website: delegates a proxy to the parent class, hence needs to reference to itself!
        def __init__(self, optimizer='Minuit'):
            """
            constructor that included listOfTargets, shadows Composite2 constructor
            :param optimizer:
            """
            self.name = None
            self.listOfTargets = []
            self.listOfRois = []
            self.configuration = None
            super(CompositeLikelihood2,self).__init__(optimizer='Minuit')
        def setConfiguration(self,config):
            if not isinstance(config,configuration):
                raise Exception("Must be a configuration instance")
            self.configuration = config
        def getListOfRoisAsString(self,sep=","):
            roi_list = ""
            for _roi in self.listOfRois:
                roi_list+=_roi.name
                if _roi.name!=self.listOfRois[-1]:
                    roi_list+=sep
            # bugfix- get rid of trailing separator
            roi_list+="\n"
            return roi_list.replace("%s\n"%sep,"")

        def getListOfRoisByName(self):
            return [_roi.name for _roi in self.listOfRois]

        def setName(self,name):
            self.name = str(name)

        def getListOfTargets(self,suppressZero=False):
            """
            the list of targets in the whole sample; targets are defined in the ROI (objects in ROI xml)
            :param suppressZero: use this keyword to suppress sources that are zero'ed
            :return: listOfTargets (list)
            """
            retList = []
            if not suppressZero:
                retList = self.listOfTargets
            else:
                for target in self.listOfTargets:
                    if target.has_key("Normalization"):
                        if float(target.Normalization)!=0:
                            retList.append(target)
                    else:
                        retList.append(target)
            return retList

        def addComponent(self, theROI):
            """
            shadows Composite2::addComponent, but keeps track of added sources.
            :param theROI:
            :raise:
            """
            if not isinstance(theROI, ROI):
                raise Exception("must be an ROI class object")
            self.listOfTargets += theROI.list_of_sources
            self.listOfRois.append(theROI)
            super(CompositeLikelihood2,self).addComponent(theROI.likelihood_fcn)

        def __append_to_class(self, my_dict):
            # to artificially extend the class, assume this to be the minos stuff etc.
            # NOT intended to be used outside class
            # to avoid clashes, use a my-prefix for all variable names
            self.__dict__.update(my_dict)


        def store_result(self, resultxml, tagName='CompLike', suppressZeroSources=False):
            """
            to store the results of the composite likelihood stuff
            okay assume that all variables we want to store have a my_ in the prefix
            NEW XML structure:
            <?xml version="1.0" ?>
            <result EGAL="/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/isotrop_2year_P76_source_v0.txt" GAL="/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/ring_2year_P76_v0.fits" dataperiod="91d" irf="P7SOURCE_V6" mass="0" path="/afs/slac/g/glast/users/zimmer/FermiData/Clusters/results/DIST4/" physicsmode="CR" runtime="Sun Mar 31 12:40:06 2013" sciencetools="09-31-00">
              <IndividualFit>
                 <ROI>
                  <par/>
                 </ROI>
              </IndividualFit>
              <CompLike>
                 <sequence
                  <par/>
                 </sequence>
              </CompLike>
            </result>
            :param resultxml: the output xml file
            :param tagName: IndividualFit or CompLike
            :param suppressZeroSources: true of false (for accounting)
            :return: LLH value
            """
            if self.configuration.Queue!="None":
                sleeptime = random.randrange(10,100)
                print("*INFO* done collecting all results, sleeping for %i seconds to avoid too many disk access at once"%sleeptime)
                time.sleep(sleeptime)
            if os.path.isfile(resultxml):
                try:
                    open(resultxml,"read")
                except xml.parsers.expat.ExpatError:
                    sleeptime = random.randrange(10,100)
                    print("*INFO* could not read %s, sleeping for %i seconds and try again."%(resultxml,sleeptime))
                    time.sleep(int(sleeptime))
            else:
                raise Exception("** FATAL ** resultfile %s not existent, run configuration::buildResultFile() first"%resultxml)
            xmlfile = xdom.parse(resultxml)
            result = xmlfile.getElementsByTagName("result")[-1]
            ParentElement = result.getElementsByTagName(tagName)[-1]
            elements = {}
            searchName = None
            if tagName=="IndividualFit":
                elements = {el.getAttribute("name"):el for el in ParentElement.getElementsByTagName("ROI")}
                searchName = self.name
            else:
                elements = {el.getAttribute("ROIs"):el for el in ParentElement.getElementsByTagName("sequence")}
                searchName = self.getListOfRoisAsString()
            currentElement = None
            new = True
            if len(elements.keys())!=0:
                if searchName in elements.keys():
                    currentElement = elements[searchName]
                    new = False
            if new:
                if tagName == "IndividualFit":
                    currentElement = xmlfile.createElement("ROI")
                else:
                    currentElement = xmlfile.createElement("sequence")
                # okay everything is defined.
            currentElement.setAttribute("NumberOfTargets","%i"%len(self.getListOfTargets(suppressZero=suppressZeroSources)))
            if tagName == "IndividualFit":
                currentElement.setAttribute("name",searchName)
            else:
                currentElement.setAttribute("ROIs",searchName)
                currentElement.setAttribute("NumberOfTiedRois","%i"%len(self.listOfRois))
            if convBool(self.configuration.nullhypothesis):
                currentElement.setAttribute("LLHNull","%1.8e"%self())
            else:
                currentElement.setAttribute("LLH","%1.8e"%self())
                # done storing ROI information. Now onto the parameters.
            if currentElement.hasAttribute("LLH") and currentElement.hasAttribute("LLHNull"):
                Ts = -2*(float(currentElement.getAttribute("LLH"))-float(currentElement.getAttribute("LLHNull")))
                attr = "Ts"
                if tagName == "CompLike":
                    attr+="_Global"
                currentElement.setAttribute(attr,"%1.4e"%Ts)
            my_keys = []
            for k in self.__dict__.keys():
                if k.startswith("my_"):
                    my_keys.append(k)
            for k in my_keys:
                parIsNew = True
                curr_dict = self.__dict__[k]
                if not isinstance(curr_dict,dict):
                    raise Exception("my_key %s must be a dictionary to proceed"%k)
                par = makeXmlElementFromDict(curr_dict,tagname="parameter")
                par.setAttribute("CL","%1.3f"%self.configuration.get("CL",float))
                # does this parameter perhaps already exist, if so, overwrite, otherwise create!
                for p in currentElement.getElementsByTagName("parameter"):
                    if p.hasAttribute("name") and curr_dict.has_key("name"):
                        if p.getAttribute("name")==curr_dict["name"]:
                            if p.hasAttribute("CL") and float(p.getAttribute("CL"))==self.configuration.get("CL",float):
                                parIsNew = False
                                for key in curr_dict:
                                    value = curr_dict[key]
                                    if type(value)==str:
                                        p.setAttribute(key,value)
                                    elif type(value)==int or type(value)==float:
                                        p.setAttribute(key,"%1.5e"%value)
                                    else:
                                        pass # do nothing.
                if parIsNew:
                    currentElement.appendChild(par)
            if new:
                ParentElement.appendChild(currentElement)

            f = open(resultxml, 'w')
            f.write(xmlfile.toprettyxml())
            f.close()
            prettifyXML(resultxml)
            return self()

        def save_covariance_matrix(self, oplot):
            if self.covariance is None:
                print('*WARNING* requested covariance plot, but matrix is empty')
                return -1
            else:
                covar_mat = np.array(self.covariance)
                h2 = toolbox.make_covariance_plot(covar_mat)
                c = R.TCanvas("c", "covariance matrix", 400, 400)
                h2.Draw("colz")
                c.SaveAs(oplot)

        def my_minos(self, configuration, likelihood, name, parameter):
            # almost done...
            d = {}
            Id = likelihood.par_index(name, parameter)
            mle = likelihood[Id].parameter.getValue()
            scale = likelihood[Id].parameter.getScale()
            # double check
            print(self)
            print(likelihood)
            print(name)
            print(parameter)
            print(configuration.CL)
            try:
                minNeg, minPos = self.minosError(likelihood, name, parameter, float(configuration.CL))
            except TypeError:
                print('*WARNING* cannot find MINOS errors')
                # just a temporary fix!
                minNeg = 0.
                minPos = 0.
                #sys.exit()
            physics_par = triple(parameter, mle, minPos, minNeg)
            physics_par.scale = scale
            d['my_sigmav'] = physics_par
            if configuration.PhysicsMode == 'DM_D':
                my_scale = d['my_sigmav'].scale
                d['my_sigmav'].scale = 1. / my_scale
                try:
                    d['my_sigmav'].value = 1. / physics_par.value
                except ZeroDivisionError:
                    d['my_sigmav'].value = 0.
                try:
                    d['my_sigmav'].neg = 1. / minNeg
                except ZeroDivisionError:
                    d['my_sigmav'].neg = 0.
                try:
                    d['my_sigmav'].pos = 1. / minPos
                except ZeroDivisionError:
                    d['my_sigmav'].pos = 0.

                try:
                    if convBool(configuration.C2_TieEgal):
                        Id = likelihood.par_index(configuration.isotropic, 'Normalization')
                        norm = likelihood[Id].parameter.getValue()
                        isotropic = my_value('EGAL', norm)
                        d['my_egal'] = isotropic
                    if convBool(configuration.C2_TieGal):
                        Id = likelihood.par_index(configuration.galactic, 'Value')
                        value = likelihood[Id].parameter.getValue()
                        galactic = my_value('GAL', value)
                        d['my_gal'] = galactic
                except RuntimeError:
                    pass

            d['likelihood_fcn'] = likelihood
            # cool now we have all stuff...
            #makes sure that the composite2 object has a common likelihood fcn as well
            self.__append_to_class(d)

        def scan_likelihood(self, config, parameter, ROIS, src, label=None):
            # pickles the info that's presently available and calls external instance ()
            logfile = ROIS[0].configuration.loghome + "/" + ROIS[0].configuration.Log + '/' + 'scanLLH_%s.log' % ROIS[
                0].name
            if os.getenv("TEMPDIR") is None:
                raise Exception("$TEMPDIR not defined, source setup.sh first!")
            tempfile = NamedTemporaryFile(prefix=os.getenv('TEMPDIR') + '/')
            dumpfile = tempfile.name
            tempfile.close()
            rois = []
            for r in ROIS:
                rois.append(r)
                # assemble dumpfile
            pickle_comp = PickleCompLike(self)
            dump_dict = {'config': config, 'mode': label, 'rois': rois, 'parameter': parameter,
                         'compPickle': pickle_comp}
            dumpstring = pickle.dumps(dump_dict, -1)

            #dumpfile = 'dump.dat'
            fdumpfile = open(dumpfile, 'w')
            fdumpfile.write(dumpstring)
            fdumpfile.close()

            command_line = ''
            # now need to call outside world
            if not config.Queue == 'None':
                queue_id = global_queues[config.Queue] - 1
                if queue_id < 0:
                    set_queue = queue(0)
                else:
                    set_queue = queue(queue_id)
                command_line = 'bsub %s -q %s -o %s ' % (bsub_opts,set_queue, logfile)
                #command_line = 'bsub -q %s '%(set_queue) # email version
            command_line += 'time python src/scan_llh.py %s' % dumpfile#.name
            print('*INFO* spawning command %s' % command_line)
            if ROI.configuration.Queue != "None":
                sleeptime = random.randrange(60, 300) # sleep between 15 and 120 seconds
                time.sleep(sleeptime)
            toolbox.runExec(command_line)


        def CALL_scan_likelihood(self, config, parameter, ROI, src, label=None):
            print('*SCAN LLH BLOCK*\nself: %s\nconfig: %s\nparameter: %s\nROI: %s'%(str(self),str(config),str(parameter),str(ROI)))
            print('ROI.__dict__: %s\nsrc: %s\nsrc.__dict__: %s\nlabel: %s'%str(ROI.__dict__),str(src),str(src.__dict__),label)
            print('****************')
            key = 'my_' + parameter
            print('*CALL* %s'%key)
            xval = float(self.__dict__[key].value)
            #print '*DEBUG* value',xval
            # that should give me the LLH to the MINOS step
            xmax = 1.3 * xval
            xmin = xval * 0.7
            if float(self.__dict__[key].pos) > 0.:
                xmax = xval + float(self.__dict__[key].pos)
                #print '*DEBUG* set max',xmax
            if float(self.__dict__[key].neg) < 0.:
                #print '*DEBUG* NegMinos',float(self.__dict__[key].neg)
                xmin = xval + float(self.__dict__[key].neg)
                #print '*DEBUG* set min',xmin
            print('*INFO* scan range min=%e max%e' % (xmin, xmax))
            #!FIXME! change to 0 later!
            current_log_val = -ROI.likelihood_fcn()
            print('*INFO* LLH value at MLE', current_log_val)
            #if not label is None:
            #    ROI.verifyByEye()
            # let's change the parameter
            Id = ROI.likelihood_fcn.par_index(src.name, parameter)
            ROI.likelihood_fcn[Id].parameter.setBounds(1e-10, 1e10) # relax parameters
            if not label is None:
                ROI.verifyByEye()
            print('src %s par %s' % (src.name, src.parameter))
            xval, dloglx = ROI.likelihood_fcn.scan(src.name, src.parameter, xmin=xmin, xmax=xmax, npts=20, verbosity=2)
            # now we need to scale dlogLx with the actual logLikelihood Value
            dlogLx = np.array(dloglx)
            dlogLx += current_log_val
            fname = config.outdir + '/'
            if not label is None:
                fname += label
            else:
                fname += 'ROI_' + ROI.name
            fname += '_' + src.parameter + '_' + config.mass + '.scan'
            print('*INFO* scanfile written to %s' % fname)
            f = open(fname, 'w')
            for i in range(len(xval)):
                #print "%f \t%f"%(xval[i],dlogLx[i])
                f.write("%f \t%f\n" % (xval[i], dlogLx[i]))
            f.close()
            # now let's add the plot
            gc.enable()
            im1 = plt.figure(figsize=(6, 6), dpi=80)
            matplotlib.rcParams['font.family'] = 'serif'
            #matplotlib.rcParams['verbose.level'] = 'silent' # to silence the ANNOYING outputs....
            #matplotlib.rcParams['verbose.fileo'] = 'matplotlib_warning.log'
            im1.set_label("Control Plot")
            plt.figtext(0.5, 0.935, 'ROI: %s' % src.name, ha='center', color='black')#, weight='bold', size='large')
            gc.enable()
            ax = im1.add_subplot(111)
            ax.plot(xval, dlogLx, '-k', label=src.name)
            #ax.set_yscale('log')
            ax.set_ylabel('-2#Delta Log[L]')
            ax.set_xlabel(self.__dict__[key].name)
            # PLOT EVERYTHING
            im1.canvas.draw()
            fname += '.png'
            im1.savefig(fname)
            plt.close()
            print('*INFO* file %s written' % fname)
            del im1
            del ax
            gc.collect() # invoke the garbage collector

            #sys.exit()

            #toolbox.stop_here()

        def MakeStackedResidualMap(self, rois, outfile, **kwargs):
            # makes only sense in combined setup
            # can use the newly created AddFits from combinefits module
            print('*INFO* Making Stacked Residual Map')
            from combinefits import AddFits

            sum_cmap = outfile.replace('_spatialResiduals.fits', '_cmap_sum.fits')
            sum_model = outfile.replace('_spatialResiduals.fits', '_model_sum.fits')
            # count maps
            count_maps = []
            model_maps = []

            length = len(rois)
            for i in range(length):
                if rois[i].countsMap is None or rois[i].modelMap is None:
                    pass
                else:
                    count_maps.append(rois[i].countsMap)
                    model_maps.append(rois[i].modelMap)

            length = len(count_maps)
            if length == 0:
                print('*WARNING* cannot stack, all maps are None')
                return 1
            elif length == 1:
                print('*WARNING* no need to stack, only 1 ROI in stack')
                return 1
            else:
                print('*INFO* count maps: ', count_maps)
                print('*INFO* model maps: ', model_maps)
                for i in range(1, length):
                    AddFits(count_maps[0], count_maps[i], sum_cmap)
                    AddFits(model_maps[0], model_maps[i], sum_model)
                    #print '*INFO* sum files',sum_cmap,sum_model
                cmap_fits = pyfits.open(sum_cmap)[0]
                model_fits = pyfits.open(sum_model)[0]
                model_array = np.array(model_fits.data)
                cmap_array = np.array(cmap_fits.data)
                residual_array = (cmap_array - model_array) / np.sqrt(model_array)
                data = np.where(-np.isnan(residual_array), residual_array, 0)
                model_fits.data = data
                if os.path.isfile(outfile):
                    os.remove(outfile)
                model_fits.writeto(outfile)
                print('*INFO* stacked residual map written %s' % outfile)
                # now we can invoke the histogram creation directly
                dummy_roi = ROI(residualMap=outfile)
                dummy_roi.MakeResidualHistogram(fitsfile=outfile)
                # and now invoke the residual histogram
                return 0

    class ROI(object):
        # still work needed here...
        def __init__(self, **kwargs):
            self.debug = 0
            self.modelxml = None
            self.modelxml_backup = None
            self.list_of_sources = []
            self.files = {}
            self.extended = False
            self.likelihood_fcn = None
            self.individual = True
            self.minuit_object = None
            self.fitQuality = None
            self.MINOS_ID = None # should not be necessary...
            self.fit_result = None
            self.eps = 1.e-5
            self.configuration = None
            self.fittedSrcModel = None
            self.roifile = None
            self.name = None
            self.sourceNames = []
            self.modelMap = None
            self.countsMap = None
            self.residualMap = None
            self.isnull = False
            self.__dict__.update(kwargs)

        def setTempDir(self,tempdir):
            if tempdir is None:
                raise Exception("must ne not none!")
            if not os.path.isdir(tempdir):
                os.system("mkdir -p %s"%tempdir)
            self.tempdir = tempdir
            print '*INFO* set temp dir to %s'%tempdir

        def getListOfTargets(self,suppressZero=False):
            """
            the list of targets in the whole sample; targets are defined in the ROI (objects in ROI xml)
            :param suppressZero: use this keyword to suppress sources that are zero'ed
            :return: listOfTargets (list)
            """
            retList = []
            if not suppressZero:
                retList = self.list_of_sources
            else:
                for target in self.list_of_sources:
                    if target.has_key("Normalization"):
                        if float(target.Normalization)!=0:
                            retList.append(target)
                    else:
                        retList.append(target)
            return retList

        def set_configuration(self, config):
            if not isinstance(config, configuration):
                raise Exception("Can only be used with valid config object")
            else:
                self.configuration = config

        def __repr__(self):
            # add repr string to make ID easier.
            kstr = "\n---------------------------------------------------------------------------------"
            kstr += "\nROI: %s" % self.name
            kstr += "\nfiles: %s" % str(self.files)
            kstr += "\nsrcModel: %s" % str(self.modelxml)
            kstr += "\nsources: %s" % str(self.sourceNames)
            kstr += "\n---------------------------------------------------------------------------------"
            return kstr

        def set_fitQuality(self, qual):
            self.fitQuality = qual

        def __getstate__(self):
            # allow pickling                                                                                                       
            d = copy.copy(self.__dict__)
            if d.has_key("likelihood_fcn"):
                del d['likelihood_fcn']
            if d.has_key("minuit_object"):
                del d['minuit_object']
            return d


        def __setstate__(self, state):
            self.__dict__ = state

        def save_covariance_matrix(self, oplot):
            if self.likelihood_fcn.covariance is None:
                print('*WARNING* requested covariance plot, but matrix is empty')
                return -1
            else:
                covar_mat = np.array(self.likelihood_fcn.covariance)
                h2 = toolbox.make_covariance_plot(covar_mat)
                c = R.TCanvas("c", "covariance matrix", 400, 400)
                h2.Draw("colz")
                c.SaveAs(oplot)

        def calculateTSvalue(self, source):
            # for a given source calculate the TSvalue
            LLH = LikelihoodState.LikelihoodState(self.likelihood_fcn)
            try:
                source.TSvalue = self.likelihood_fcn.Ts(source.name, reoptimize=True)
            except RuntimeError:
                print('*CAUGHT RUNTIME ERROR, for TS value %s return -999' % source.name)
                source.TSvalue = -999
            LLH.restore()
            return source.TSvalue

        def exportFitResultToDict(self):
            gtlikeXml = NamedTemporaryFile(prefix=os.getenv("TEMPDIR") + "/").name
            self.likelihood_fcn.writeXml(gtlikeXml)
            d = xmltools.unmarshalGtlikeXml(gtlikeXml)
            os.remove(gtlikeXml) # cleaning up
            return d

        def export_TSvalue(self):
            # now this is a new interesting feature: we write one continuous xml file with all sources that we
            # have in our likelihood and their associated TS values.... who knows what that may be useful for
            sources = self.list_of_sources
            GAL = source(name=str(self.configuration.galactic), Type="diffuse_template")
            EGAL = source(name=str(self.configuration.isotropic), Type="diffuse_template")
            sources.append(GAL)
            sources.append(EGAL)
            kstring = '<?xml version="1.0" ?>\n<some_instance/>'
            xmlf = xdom.parse(StringIO.StringIO(kstring))
            roi_el = xmlf.createElement("ROI")
            LLH = self.likelihood_fcn()
            roi_el.setAttribute("name", str(self.name))
            roi_el.setAttribute("LLH", str(LLH))
            # and now we add one source line there...
            if not convBool(self.configuration.nullhypothesis):
                for s in sources:
                    src = xmlf.createElement("source")
                    src.setAttribute("name", s.name)
                    TS = self.calculateTSvalue(s)
                    src.setAttribute("TS", str(TS))
                    roi_el.appendChild(src)
                    # returns an xml-dom element!
            return roi_el

        def readSourceList(self,profile=True):
            # reads the roifile and looks for source parameters
            xmlfile = xdom.parse(self.roifile)
            for roi in xmlfile.getElementsByTagName("ROI"):
                srcmodel = os.path.expandvars(str(roi.getAttribute("srcModel")))
                if roi.getAttribute("name") == self.name:
                    for src in roi.getElementsByTagName("object"):
                        keys = src.attributes.items()
                        d = {}
                        # cool - now we just need to loop through this stuff of the format:
                        #[(key,value),(key,value)]
                        for k in keys:
                            d[str(k[0])] = str(k[1])
                        d['srcModel'] = xdom.parse(srcmodel)
                        thesource = source()
                        subtask_output = SearchForScaleFactors(self.modelxml, src.getAttribute("name"))
                        self.extended = subtask_output[2]
                        thesource.scalefactorfromFITS = subtask_output[0]
                        if self.extended:
                            # here we add a profile plot SZ 2012-02-14
                            thesource.name = src.getAttribute("name")
                            if profile:
                                thesource.makeProfilePlot(fitsfile=subtask_output[1], output=self.configuration.outdir)
                        thesource.expand(d)
                        if d['Type'] == 'DMFit':
                            thesource.DMFitSource()
                        elif d['Type'] == 'DMFit::ScaleFactor':
                            thesource.DMFitSource(ScaleFactor=True)
                        elif d['Type'] == 'FileSpectrum':
                            thesource.FileSource()
                        elif d['Type'] == 'FileSpectrum::ScaleFactor':
                            thesource.FileSource(ScaleFactor=True)
                            #print '*DEBUG* ',thesource.__dict__
                        elif d['Type'] == 'PowerLaw2':
                            thesource.PowerLaw2()
                        elif d['Type'] == 'PowerLaw2::ScaleFactor':
                            thesource.PowerLaw2(ScaleFactor=True)
                            #print '*DEBUG* ',thesource.__dict__

                        elif d['Type'] == 'Generic':
                            thesource.GenericSource()
                            # a really generic source, no fancy stuff when setting the parameters
                        #print thesource.__dict__
                        self.list_of_sources.append(thesource)
                        self.sourceNames.append(thesource.name)

        def fitNull(self,export_fit=False,cleanup=True):
            bestModel = os.path.join(self.tempdir,"model.xml")
            nullModel = os.path.join(self.tempdir,"nullModel.xml")
            oldModel = self.modelxml
            shutil.copy(oldModel,bestModel)
            print '*** DOING NULLFIT ***'
            self.make_nullfit()
            shutil.copy(self.modelxml,nullModel)
            self.modelxml = nullModel
            self._prepare()
            self.configuration.nullhypothesis = "True"
            self.fit(None)
            LLHNull = self.likelihood_fcn()
            _dict = self.exportFitResultToDict()
            # got back to best fit setup
            self.modelxml = bestModel
            self._prepare()
            self.configuration.nullhypothesis = "False"
            if cleanup:
                self.cleanup()
            if export_fit:
                return (LLHNull,_dict)
            else:
                return LLHNull,None

        def fit(self, mysource):
            print('*INFO* now entering FIT ')
            #print self.minuit_object
            #print self
            # okay here's the clue - for individual fit - put J-values to 0 for all stuff...
            # done....
            d = {}
            #print self.likelihood_fcn
            #sys.exit()
            # assume that we don't want individual minos stuff..., fit the roi with all sources but get parameters...
            print('** FIT OF ALL SOURCES IN ONE ROUND **')
            self.likelihood_fcn.fit(covar=True, tol=self.eps, optimizer=self.configuration.optimizer,
                                    optObject=self.minuit_object)
            if convBool(self.configuration.nullhypothesis):
                return 0
            else:
                for src in self.list_of_sources:
                    Id = self.likelihood_fcn.par_index(src.name, src.parameter)
                    mle = self.likelihood_fcn[Id].parameter.getValue()
                    sigmav = triple(src.parameter, mle, 0, 0)
                    sigmav.scale = self.likelihood_fcn[Id].parameter.getScale() # store the scale...
                    myKey = 'my_sigmav_%s' % src.name
                    d[myKey] = sigmav
                    if src.Type == 'FileSpectrum::ScaleFactor' or src.Type == 'DMFit':
                        if self.configuration.PhysicsMode == 'DM_D':
                            sca = d[myKey].scale
                            d[myKey].scale = 1. / sca
                            try:
                                d[myKey].value = 1. / sigmav.value
                            except ZeroDivisionError:
                                d[myKey].value = 0.

                try:
                    # isotropic diffuse
                    Id = self.likelihood_fcn.par_index(self.configuration.isotropic, 'Normalization')
                    norm = self.likelihood_fcn[Id].parameter.getValue()
                    # galactic diffuse
                    Id = self.likelihood_fcn.par_index(self.configuration.galactic, 'Value')
                    value = self.likelihood_fcn[Id].parameter.getValue()
                    gal = my_value('gal', value)
                    egal = my_value('egal', norm)
                    d['my_gal'] = gal
                    d['my_egal'] = egal

                except RuntimeError:
                    # no diffuse sources included... no problems...
                    pass
                mysource.expand(d) # that associates the values from this section with source parameters
                # done
            return self.likelihood_fcn() 
        
        def _prepare(self):
            print '*INFO* source model %s'%self.modelxml
            likeObs = BinnedObs(srcMaps=self.files["srcmap"], expCube=self.files["expcube"],
                                binnedExpMap=self.files["binnedExpMap"], irfs=self.configuration.IRF)
            self.likelihood_fcn = BinnedAnalysis(likeObs, srcModel=self.modelxml)
            old_verbosity = self.likelihood_fcn.verbosity
            self.likelihood_fcn.verbosity = 3
            print ('*INFO* Energy Range for Analysis %1.1f -- %1.1f MeV'%(float(self.configuration.EMin), float(self.configuration.EMax)))
            self.likelihood_fcn.setEnergyRange(float(self.configuration.EMin), float(self.configuration.EMax))
            print('*INFO* Energy Range for Analysis %1.1f -- %1.1f MeV' %(float(self.configuration.EMin), float(self.configuration.EMax)))
            self.likelihood_fcn.verbosity = old_verbosity
            self.likelihood_fcn.optimizer = self.configuration.optimizer
            #print '*have setup a Binned Analysis object*',self.likelihood_fcn
            self.minuit_object = eval("pyLike.%s(self.likelihood_fcn.logLike)" % self.configuration.optimizer)
            print(self.likelihood_fcn)
            #print self.minuit_object
            #sys.exit()
            self.minuit_object.setStrategy(2)
            self.likelihood_fcn.setFitTolType(1)
            #print '*did the eval of pyLike*'

        def make_nullfit(self):
            if self.modelxml is None:
                raise Exception("Cannot be called without a source model")
            print('*** REQUESTING NULL FIT ***')
            self.isnull = True
            xf = xdom.parse(self.modelxml)
            for src in xf.getElementsByTagName("source"):
                srcname = str(src.getAttribute("name"))
                for source in self.list_of_sources:
                    #print '*working on %s*'%source()
                    if srcname == source.name:
                        #print '*match source* %s'%source.name
                        for par in src.getElementsByTagName("parameter"):
                            parname = str(par.getAttribute("name"))
                            if parname == source.parameter:
                                par.setAttribute("value", str(0))
                                par.setAttribute("min", str(-5))
                                par.setAttribute("max", str(+5))
                                par.setAttribute("scale", str(1))
                                par.setAttribute("free", str(0))
                                #print src.toxml()
            self.modelxml_backup = self.modelxml
            if os.getenv("TEMPDIR") is None:
                raise Exception("$TEMPDIR not defined, source setup.sh first!")
            self.modelxml = NamedTemporaryFile(prefix=os.getenv("TEMPDIR") + "/").name
            print('*INFO* using temporary source model %s' % self.modelxml)
            f = open(self.modelxml, "w")
            f.write(xf.toxml())
            f.close()

        def cleanup(self):
            if self.isnull:
                print('*INFO* removing temporary xml model %s' % self.modelxml)
                os.remove(self.modelxml)
                self.modelxml = self.modelxml_backup
            del self.likelihood_fcn 
            del self.minuit_object           
            return


        def load_fitresult(self, freeze=False):
            # this one loads the saved xml and interprets it as input xml
            saveVersionXML = self.configuration.likeModelVersion
            saveVersionXML += '_' + str(self.name) + '_m' + self.configuration.mass + '_' + self.configuration.mode
            fittedSrcModel = self.configuration.outdir + saveVersionXML + '_fittedSrc.xml'

            if os.path.exists(fittedSrcModel) and os.path.getsize(fittedSrcModel) > 0:
                sourceNames = []
                # need a bit of code here
                if freeze:
                    shutil.copy(fittedSrcModel, fittedSrcModel + ".bak") # make a backup
                    sourceNames = self.sourceNames + ["GAL", "EGAL"]
                    print('*INFO* created backup of fitted xml model %s.bak' % fittedSrcModel)
                    # finally load fit result
                xm = xdom.parse(fittedSrcModel)
                if freeze:
                    for src in xm.getElementsByTagName("source"):
                        if not src.getAttribute("name") in sourceNames:
                            for p in src.getElementsByTagName("parameter"):
                                p.setAttribute("free", "0")
                    f = open(fittedSrcModel, "w")
                    f.write(xm.toxml())
                    f.close()
                    print('*INFO* all sources except GAL,EGAL & targets frozen')

                self.modelxml = fittedSrcModel
                self._prepare()

                fitqual = str(xm.childNodes[0].getAttribute("fitQuality"))
                if fitqual != 'None':
                    int_fitqual = int(fitqual)
                    self.set_fitQuality(int_fitqual)
            else:
                raise Exception("*FATAL* cannot find fit result in xml format, check content of %s" % fittedSrcModel)

        def prepare_and_set_fit(self, mass, mode, newDM=False):
            # to combine and refactor prepare and set
            print('*INFO* new prepare_and_set_fit method')
            if convBool(self.configuration.nullhypothesis):
                self.make_nullfit()
            self._prepare()
            if convBool(self.configuration.nullhypothesis):
                return 0
            #for src in self.list_of_sources:
            #    print '*DEBUG* *MINOS_ID: * src %s Id %i'%(src.name,src.MINOS_ID)
            # can modify the likelihood values for each source and their mass....
            #print '*did the eval of pyLike*'
            else:
                for src in self.list_of_sources:
                    if src.is_free == "1":
                        src.MINOS_ID = xmltools.GetParameterIDForLikelihood(self.modelxml, src.name)

                for src in self.list_of_sources:
                    if src.is_free == "1":
                        src.MINOS_ID = xmltools.GetParameterIDForLikelihood(self.modelxml, src.name)
                    # can modify the likelihood values for each source and their mass....
                #self.verifyByEye()
                if newDM:
                    print('*INFO* requested new DM/CR implementation from FileFunction')
                    #toolbox.construction_site()
                else:
                    #                 print '* use old implementation*'
                    #                 print '**** DEBUG ****'
                    #                 self.verifyByEye()
                    #                 print '**** DEBUG ****'
                    for src in self.list_of_sources:
                        #print 'src name',src.name
                        #print '**** DEBUG ****',src.__dict__
                        if "ScaleFactor" in src.Type or src.Type == 'DMFit': # enable scalefactor class to be used much more widely
                            if src.Type == 'DMFit':
                                # set final state accordingly
                                finalStates = {'ee': 1, 'mumu': 2, 'tautau': 3, 'bbbar': 4, 'ww': 7, 'zz': 8}
                                try:
                                    FSTATE = round(float(finalStates[self.configuration.FinalState]), 0)
                                    src.xmlpars['channel0'].value = FSTATE
                                except KeyError:
                                    raise Exception(
                                        "Fatal: FinalState %s parameter in configuration xml not understood, supported types %s" % (
                                            self.configuration.FinalState, str(finalStates)))
                            elif src.Type == 'FileSpectrum::ScaleFactor':

                                # todays feature: normalization valid above low energy limit in
                                # filefunction, needs to be rescaled if energy is changed
                                for s in src.srcModel.getElementsByTagName("source"):
                                    if s.getAttribute("name") == src.name:
                                        xml_spectrum = s.getElementsByTagName("spectrum")[-1]
                                        if not xml_spectrum.hasAttribute("file"):
                                            raise Exception(
                                                "*FATAL* source %s is designated FileSpectrum but has no attribute file" % src.name)
                                        else:
                                            fname = toolbox.resolveFilename(xml_spectrum.getAttribute("file"))
                                            # candidate for removal - does some iffy stuff...
                                            #old_integral = toolbox.integrate(fname) # should be 1
                                            #new_integral = toolbox.integrate(fname,emin=self.configuration.EMin,emax=self.configuration.EMax)
                                            #scalefactor_spectrum = new_integral/old_integral
                                            #if scalefactor_spectrum!=1:
                                            #    print '*INFO* energy range has changed, need to rescale %s by %1.1f'%(src.name,scalefactor_spectrum)
                                            #src.xmlpars['Normalization'].value*=scalefactor_spectrum

                            if mode.startswith("DM"):
                                if mode == 'DM_D':
                                    # double check!
                                    src.xmlpars['mass'].setValue(0.5 * float(mass))
                                    #src.xmlpars['mass'].make_bounds()
                                    j_value = src.xmlpars['norm'].value
                                    src.xmlpars['norm'].setValue(0.5 * float(mass) * j_value * 2.)
                                    src.xmlpars['norm'].make_bounds()
                                else:
                                    if src.xmlpars.has_key('mass'):
                                        src.xmlpars['mass'].setValue(mass)
                                    else:
                                        print('*WARNING* caught KeyError with in src.xmlpars[mass], likelihood_class.py:1103')
                                        #src.xmlpars['mass'].make_bounds()
                                        #j_value = src.xmlpars['norm'].value
                                        #src.xmlpars['norm'].value=j_value#/2.#(8.*pi) # is that *REALLY* correct?
                                        #src.xmlpars['norm'].make_bounds()
                        for par in src.xmlpars.keys():
                            Id = self.likelihood_fcn.par_index(src.name, par)
                            if src.is_free != "1":
                                self.likelihood_fcn[Id].parameter.setFree(False)
                            if int(self.configuration.chatter) > 4:
                                print('*DEBUG: * ID %i name: %s parameter: %s' % (Id, src.name, par))
                                print('processing parameter: %s'%str(src.xmlpars[par].__dict__))
                            if src.xmlpars[par].name == 'channel0':
                                self.likelihood_fcn[Id].parameter.setScale(1)
                                self.likelihood_fcn[Id].parameter.setBounds(0, 12)
                                self.likelihood_fcn[Id].parameter.setValue(src.xmlpars['channel0'].value)

                            if isinstance(src.xmlpars[par], triple):
                                if int(self.configuration.chatter) > 5:
                                    print('***** DEBUG ****** %s'%str(src.__dict__))
                                    print('***** DEBUG ****** %s'%str(src.xmlpars[par].__dict__))
                                    print('***** DEBUG ****** current val %s'%self.likelihood_fcn[Id].parameter.getValue())
                                    print('***** DEBUG ****** current scale %s'%self.likelihood_fcn[Id].parameter.getScale())
                                if not src.xmlpars[par].name == 'sigmav':
                                    src.xmlpars[par].make_bounds()
                                    #print 'new bounds ',src.xmlpars[par].min,src.xmlpars[par].max
                                    #try:
                                    #self.likelihood_fcn[Id].parameter.setBounds(src.xmlpars[par].min,src.xmlpars[par].max)
                                self.likelihood_fcn[Id].parameter.setScale(src.xmlpars[par].scale)
                                #print 'attempting to set this triple',src.xmlpars[par].__dict__ 
                                try:
                                    self.likelihood_fcn[Id].parameter.setValue(float(src.xmlpars[par].value))
                                    self.likelihood_fcn[Id].parameter.setBounds(src.xmlpars[par].min,
                                                                                src.xmlpars[par].max)
                                except RuntimeError:
                                    #print ' change order to prevent crash '
                                    print(src.xmlpars[par].min, src.xmlpars[par].value, src.xmlpars[par].max)
                                    print(self.likelihood_fcn[Id].parameter.getBounds())
                                    print(self.likelihood_fcn[Id].parameter.getValue())

                                    self.likelihood_fcn[Id].parameter.setBounds(src.xmlpars[par].min,
                                                                                src.xmlpars[par].max)
                                    self.likelihood_fcn[Id].parameter.setValue(float(src.xmlpars[par].value))
                            if src.xmlpars[par].name == 'norm' and "DMFit" in src.Type and src.xmlpars[par].error != 0:
                                self.likelihood_fcn[Id].addPrior("LogNormalLog")
                                self.likelihood_fcn[Id].setPriorParams(Log10_Mean=np.log10(src.xmlpars['norm'].value),
                                                                       Log10_Sigma=src.xmlpars[par].error)
                                self.likelihood_fcn[Id].setFree(1)
                                print('set J prior with Log10_Mean=%1.4e'%np.log10(src.xmlpars['norm'].value))
                                print('and Log10_Sigma=%1.4e'%src.xmlpars[par].error)

            print(str(self.likelihood_fcn))
            #print self.minuit_object
            self.verifyByEye()

        def verifyByEye(self, par=None):
            # use to check the values of the fit
            print(str(self.likelihood_fcn))
            print('***** LIKELIHOOD VALUES FOR THE FIT NOW *******')
            continue_print = True
            ii = 0
            while continue_print:
                try:
                    name = self.likelihood_fcn[ii].parameter.getName()
                    if par == name:
                        print('id\t%i\t free %s\t name %s\t\t value %f\t scale %1.1e\t bounds %s' % (
                            ii, str(self.likelihood_fcn[ii].parameter.isFree()),
                            str(self.likelihood_fcn[ii].parameter.getName()),
                            self.likelihood_fcn[ii].parameter.getValue(),
                            float(self.likelihood_fcn[ii].parameter.getScale()),
                            str(self.likelihood_fcn[ii].parameter.getBounds())))

                    elif par is None:
                        print('id\t%i\t free %s\t name %s\t\t value %f\t scale %1.1e\t bounds %s' % (
                            ii, str(self.likelihood_fcn[ii].parameter.isFree()),
                            str(self.likelihood_fcn[ii].parameter.getName()),
                            self.likelihood_fcn[ii].parameter.getValue(),
                            float(self.likelihood_fcn[ii].parameter.getScale()),
                            str(self.likelihood_fcn[ii].parameter.getBounds())))
                        ii += 1
                except AttributeError:
                    print('index failure: %i'%ii)
                    continue_print = False


        def plot(self, srcname, save=True):
            counts_Spectra = "%s/%s_counts_spectra.fits" % (self.configuration.outdir, srcname)
            if self.configuration.mass != 0:
                if self.configuration.__dict__.has_key("FinalState"):
                    counts_Spectra = "%s/%s_counts_spectra_m%s_%s.fits" % (
                        self.configuration.outdir, srcname, self.configuration.mass, self.configuration.FinalState)
                else:
                    counts_Spectra = "%s/%s_counts_spectra_m%s.fits" % (
                        self.configuration.outdir, srcname, str(self.configuration.mass))
            self.likelihood_fcn.writeCountsSpectra(counts_Spectra)

        def MakeCombinedSpectralPlot(self, mysource):
            spectrum = self.configuration.outdir + mysource.name + '_spectral_' + self.configuration.specialTag + '.root'
            residual = self.configuration.outdir + mysource.name + '_residual_' + self.configuration.specialTag + '.root'
            out = self.configuration.outdir + mysource.name + '_comb_spec_' + self.configuration.specialTag
            if not (os.path.isfile(spectrum) or os.path.isfile(residual)):
                print('*WARNING* could not combine spectral plots, check rootfiles')
                return 0
            file1 = R.TFile(spectrum, "READ")
            file2 = R.TFile(residual, "READ")
            canvas = R.TCanvas("canvas", "a canvas", 500, 750)
            canvas.Divide(1, 2)
            # lets get the pad
            #print '*DEBUG*'
            #file1.GetListOfKeys().ls()
            #file2.GetListOfKeys().ls()
            #print '*DEBUG*'

            pad1 = file1.GetListOfKeys().Last().ReadObj()
            pad2 = file2.GetListOfKeys().Last().ReadObj()
            canvas.cd(1)
            pad1.DrawClonePad()
            canvas.cd(2)
            pad2.DrawClonePad()
            canvas.SaveAs(out + '.png')
            file1.Close()
            file2.Close()

        def MakeSpatialResidualPlot(self, run=True,sleeptime=None):
            # pickles the info that's presently available and calls external instance ()
            logfile = self.configuration.loghome + self.configuration.Log + '/' + 'runMakeSpatialResidual_%s.log' % self.name
            if os.getenv("TEMPDIR") is None:
                raise Exception("$TEMPDIR not defined, source setup.sh first!")
            tempfile = NamedTemporaryFile(prefix=os.getenv('TEMPDIR') + '/')
            dumpfile = tempfile.name
            tempfile.close()

            d = {"roi": self, "configuration": self.configuration}
            dumpstring = pickle.dumps(d, -1)

            #dumpfile = 'dump.dat'
            fdumpfile = open(dumpfile, 'w')
            fdumpfile.write(dumpstring)
            fdumpfile.close()

            command_line = ''
            # now need to call outside world
            if not self.configuration.Queue == 'None':


                if global_queues[self.configuration.Queue] >= 3:
                    queue_id = 3 # use long for spatial maps...
                else:
                    queue_id = 2
                if queue_id < 0:
                    set_queue = queue(0)
                else:
                    set_queue = queue(queue_id)
                command_line = 'bsub %s -q %s -o %s ' % (bsub_opts, set_queue, logfile)
                #command_line = 'bsub -q %s '%(set_queue) # email version
            command_line += 'time python src/MakeSpatialResiduals.py %s' % dumpfile#.name
            print('*INFO* spawning command %s' % command_line)

            if self.configuration.Queue != "None":
                if sleeptime is None:
                    sleeptime = random.randrange(60, 300) # sleep between 15 and 120 seconds
                time.sleep(sleeptime)
            toolbox.runExec(command_line)

        def CALL_MakeSpatialResidualPlot(self, run=True, **kwargs):
            #print "****D E B U G ****"
            #print self.__dict__
            # that code is called from elsewhere
            print '*SRC MODEL:',self.fittedSrcModel
            if not convBool(self.configuration.MakeMaps):
                return 0
                # uses a GtApp for gtmodelmaps and does fits magic
            #from GtApp import GtApp

            gc.enable()
            app = GtApp('gtmodel')
            self.modelMap = self.configuration.outdir + self.name + '_SpatialModel_' + self.configuration.specialTag + '.fits'
            self.countsMap = self.files['cmap']
            self.residualMap = self.configuration.outdir + self.name + '_SpatialResidual_' + self.configuration.specialTag + '.fits'
            if not run:
                return 0
            toolbox.setPF()
            app.run(srcmaps=self.files['srcmap'],
                    srcmdl=self.fittedSrcModel,
                    outfile=self.modelMap,
                    irfs=self.configuration.IRF,
                    expcube=self.files['expcube'],
                    bexpmap=self.files['binnedExpMap'],
                    chatter=4,
                    **kwargs)
            toolbox.cleanPF()
            # now we are almost done
            if os.path.isfile(self.residualMap):
                os.remove(self.residualMap)

            counts_fits = pyfits.open(self.countsMap)[0]
            model_fits = pyfits.open(self.modelMap)[0]

            counts_data = np.array(counts_fits.data)
            model_data = np.array(model_fits.data)

            diff_data = counts_data - model_data

            # now in this residual map we can interpret the stuff as sigma-maps!
            res_data = diff_data / np.sqrt(model_data)
            data = np.where(-np.isnan(res_data), res_data, 0)
            model_fits.data = data
            model_fits.writeto(self.residualMap)
            print('*INFO* spatial residual file written: %s' % self.residualMap)
            # now we can invoke the histogram creation directly
            self.MakeResidualHistogram(data=data)
            gc.collect()

        def MakeResidualHistogram(self, data=None, fitsfile=None, makeROOT=False):
            # set makeROOT to False if you don't wish to create ROOT files
            # in this design we could use it in external code too

            # digression:
            # from likelihood_class import ROI
            # roi = ROI(name="empty",residualMap="some_file.fits")
            # roi.MakeResidualHistogram(fitsfile=residualMap)

            print('*INFO* create histogram from residual map')
            Data = None
            histfile = self.residualMap.replace(".fits", "_hist.root")
            if (data is None) and (fitsfile is None):
                raise Exception("*ERROR* Need to call function with argument for either a np.array or fitsfile")
            elif data is None:
                hdu = pyfits.open(fitsfile)[0]
                Data = np.array(hdu.data)
            else:
                Data = np.array(data)
            xdim, ydim = Data.shape
            # the ROOT stuff
            c = R.TCanvas("c", "A Canvas", 400, 400)
            histogram = R.TH1D("Spatial Residuals", "", 25, -5, 5)
            histogram.GetXaxis().SetTitle("Residuals [#sigma]")
            histogram.SetLineStyle(1) # solid
            histogram.SetLineColor(1)
            func = R.TF1("func", "gaus", -5, 5)
            func.SetLineStyle(2) # dashed
            func.SetLineColor(1)

            for x in range(xdim):
                for y in range(ydim):
                    value = Data[x][y]
                    histogram.Fill(value)
                # now fit to Gaussian
            print('*INFO* fitting residuals to Gaussian')
            histogram.Fit(func, "R")
            # canvas
            c.cd()
            R.gStyle.SetOptFit(1111)
            histogram.DrawNormalized()
            # debug
            histogram.SetTitle(histfile)
            #func.Draw("lsame")
            c.SaveAs(histfile.replace(".root", ".png"))
            # now let's write the root file
            if makeROOT:
                rfile = R.TFile(histfile, "recreate")
                histogram.SetName("histogram")
                histogram.Write()
                func.Write()
                rfile.Write()
                rfile.Close()
                print('*INFO* residual histogram written: %s' % histfile)
            c.Delete()
            histogram.Delete()
            func.Delete()

        def print_summary(self, mysource):
            print('*INFO* SUMMARY OF OBSERVATIONS, SOURCE: %s' % mysource.name)
            try:
                print("Npred Isotropic Diffuse : %1.4e"%self.likelihood_fcn.logLike.NpredValue(self.configuration.isotropic))
                if convBool(self.configuration.calculateSourceTs):
                    print("TS value: %1.4e"%self.likelihood_fcn.Ts(self.configuration.isotropic, reoptimize=True))
                print("Npred Galactic Diffuse: %1.4e"%self.likelihood_fcn.logLike.NpredValue(self.configuration.galactic))
                if convBool(self.configuration.calculateSourceTs):
                    print("TS value: %1.4e"%self.likelihood_fcn.Ts(self.configuration.galactic, reoptimize=True))

            except RuntimeError:
                pass
            if not convBool(self.configuration.nullhypothesis):
                print("Npred DM                : %1.4e"%self.likelihood_fcn.logLike.NpredValue(mysource.name))
                if convBool(self.configuration.calculateSourceTs):
                    print("TS value                : %1.4e"%self.likelihood_fcn.Ts(mysource.name, reoptimize=True))
            print("-log(Like)              : %1.8e"%self.likelihood_fcn())
            print("Nobs                    : %s"%str(self.likelihood_fcn.nobs))
            return None

        def writeXML(self, srcmodel=None):
            # # that takes all the sources and employs the writeXML method in the BinnedAnalysis
            # saveVersionXML = self.configuration.likeModelVersion
            # if (mode.startswith("DM")):
            #     saveVersionXML+='_m'+self.configuration.mass+'_'+mode+
            # saveVersionXML = self.configuration.likeModelVersion+'_m'+x
            if not srcmodel is None:
                self.fittedSrcModel = srcmodel
            else:
                self.configuration.saveVersionXML = self.configuration.likeModelVersion
                self.configuration.saveVersionXML += '_' + str(self.name) + '_m' + self.configuration.mass + '_' + self.configuration.mode
                self.fittedSrcModel = self.configuration.outdir + self.configuration.saveVersionXML + '_fittedSrc.xml'
            LLH = self.likelihood_fcn()
            self.likelihood_fcn.writeXml(xmlFile=self.fittedSrcModel)
            # and afterwards store fit quality
            xm = xdom.parse(self.fittedSrcModel)
            xm.childNodes[0].setAttribute("fitQuality", str(self.fitQuality))
            xm.childNodes[0].setAttribute("LLH", str("%1.8e" % LLH)) # that exports the LLH value for the fit
            xm.childNodes[0].setAttribute("name",self.name)
            # if this option is set to False, no Ts calculations are performed.
            if convBool(self.configuration.calculateSourceTs):
                for src in xm.getElementsByTagName("source"):
                    srcname = str(src.getAttribute("name"))
                    src.setAttribute("Ts","%1.4e"%self.likelihood_fcn.Ts(srcname,reoptimize=True))
            f = open(self.fittedSrcModel, 'w')
            f.write(xm.toxml())
            print('*INFO* Save fitted model file as %s' % self.fittedSrcModel)

        def store_result(self, resultxml, suppressZeroSources=False):
            """
            to store the results of the composite likelihood stuff
            okay assume that all variables we want to store have a my_ in the prefix
            NEW XML structure:
            <?xml version="1.0" ?>
            <result EGAL="/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/isotrop_2year_P76_source_v0.txt" GAL="/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/ring_2year_P76_v0.fits" dataperiod="91d" irf="P7SOURCE_V6" mass="0" path="/afs/slac/g/glast/users/zimmer/FermiData/Clusters/results/DIST4/" physicsmode="CR" runtime="Sun Mar 31 12:40:06 2013" sciencetools="09-31-00">
              <IndividualFit>
                 <ROI>
                  <par/>
                 </ROI>
              </IndividualFit>
              <CompLike>
                 <sequence
                  <par/>
                 </sequence>
              </CompLike>
            </result>
            :param resultxml: the output xml file
            :param tagName: IndividualFit or CompLike
            :param suppressZeroSources: true of false (for accounting)
            :return: LLH value
            """
            if self.configuration.Queue!="None":
                sleeptime = random.randrange(10,100)
                print("*INFO* done collecting all results, sleeping for %i seconds to avoid too many disk access at once"%sleeptime)
                time.sleep(sleeptime)
            if os.path.isfile(resultxml):
                try:
                    open(resultxml,"read")
                except xml.parsers.expat.ExpatError:
                    sleeptime = random.randrange(10,100)
                    print("*INFO* could not read %s, sleeping for %i seconds and try again."%(resultxml,sleeptime))
                    time.sleep(int(sleeptime))
            else:
                raise Exception("** FATAL ** resultfile %s not existent, run configuration::buildResultFile() first"%resultxml)
            xmlfile = xdom.parse(resultxml)
            result = xmlfile.getElementsByTagName("result")[-1]
            ParentElement = result.getElementsByTagName("IndividualFit")[-1]
            elements = {el.getAttribute("name"):el for el in ParentElement.getElementsByTagName("ROI")}
            searchName = self.name
            currentElement = None
            new = True
            if len(elements.keys())!=0:
                if searchName in elements.keys():
                    currentElement = elements[searchName]
                    new = False
            if new:
                currentElement = xmlfile.createElement("ROI")
                currentElement.setAttribute("name",searchName)

            currentElement.setAttribute("NumberOfTargets","%i"%len(self.getListOfTargets(suppressZero=suppressZeroSources)))
            if convBool(self.configuration.nullhypothesis):
                currentElement.setAttribute("LLHNull","%1.8e"%self.likelihood_fcn())
            else:
                currentElement.setAttribute("LLH","%1.8e"%self.likelihood_fcn())
                # done storing ROI information. Now onto the parameters.
            if new:
                ParentElement.appendChild(currentElement)

            f = open(resultxml, 'w')
            # almost done... now just writing everything...
            f.write(xmlfile.toprettyxml())
            f.close()
            prettifyXML(resultxml)
            return self.likelihood_fcn()

    def read_rois(configuration, combined=False, read_sources=True,selected=None):
        # if combined=True -> return ROIs + ordering
        order = None
        xmlfile = xdom.parse(configuration.ROIFile)
        ROIS = []

        for roi in xmlfile.getElementsByTagName("ROI"):
            extended = False
            if roi.hasAttribute("extended"):
                if roi.getAttribute("extended") == "True":
                    extended = True
            else:
                extended = False
            keep = True
            name = str(roi.getAttribute("name"))
            if not selected is None:
                if name in selected:
                    keep = True
                else: 
                    keep = False
            if configuration.cdate == 'None':
                Name = name
            else:
                Name = name + '-' + configuration.cdate.partition("m")[0]
            files = {}
            # the files we need
            try:

                files['srcmap'] = glob.glob(configuration.prepdir + Name + '_BinnedSrcMap.fits*')[0]
                files['expcube'] = configuration.LightCube #glob.glob(configuration.prepdir+'expCube_'+Name+'.fits*')[0]
                files['binnedExpMap'] = glob.glob(configuration.prepdir + Name + '_expCube2.fits*')[0]
                # some additionals
                files['ccube'] = glob.glob(configuration.prepdir + Name + '_cCube.fits*')[0]
                files['cmap'] = glob.glob(configuration.prepdir + Name + '_cMap.fits*')[0]
            except IndexError:
                print('*COULD NOT FIND FILES FOR %s in %s' % (Name, configuration.prepdir))
                keep = False
                # ....
            newROI = None
            if keep:
                newROI = ROI(
                    name=name,
                    RA=str(roi.getAttribute("RA")),
                    DEC=str(roi.getAttribute("DEC")),
                    configuration=configuration,
                    modelxml=os.path.expandvars(str(roi.getAttribute("srcModel"))),
                    roifile=str(configuration.ROIFile),
                    eps=1e-5,
                    extended=extended,
                    files=files
                   )
                if read_sources:
                    newROI.readSourceList() # to build up the sources for each roi
                ROIS.append(newROI)
            #     print ROIS
            #     sys.exit()

        if combined:
            order = [None, None]
            list_of_roi = xmlfile.getElementsByTagName("list_of_roi")[-1]
            str_list = None
            if list_of_roi.hasAttribute("ordering"):
                str_list = str(list_of_roi.getAttribute("ordering")).split(",")
            order[0] = str_list
            roi_dict = {}
            for roi in ROIS:
                roi_dict[roi.name] = roi
            order[1] = roi_dict
            return ROIS, order
        else:
            return ROIS


