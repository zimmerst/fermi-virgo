import virgo_analysis as VA
import sys, os, yaml
from tempfile import NamedTemporaryFile
from xmltools import expandEnvVarsInXml
    
class VirgoContainer(object):
    def __init__(self,**kwargs):
        self.name = None
        self.J = None
        self.fits = None
        self.gal = "/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/ring_2year_P76_v0.fits"
        self.egal = "/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v2r0/isotrop_2year_P76_source_v0.txt"
        self.diffuseTag = None
        self.DiffDir = "/afs/slac/g/glast/users/zimmer/FermiData/Diffuse/SNRCatalog/"
        self.split = False
        self.Xml = None
        self.__dict__.update(kwargs)
        if not self.diffuseTag is None:
            self.gal = self.diffuseTag
            self.egal= "%s_v2_clean_isotropic.txt"%self.diffuseTag
            if self.split: 
                self.egal = os.path.join(self.DiffDir,self.egal)
                self.DiffDir = ""
                self.gal = self.diffuseTag
            else:          
                self.gal = "%s_v2.fits"%self.diffuseTag
    def makeDiffuseTag_v2(self):
        self.gal = "%s_v2.fits"%self.diffuseTag
    def dumpXml(self):
        # this creates a temporary source model to be used instead.
        os.environ['GAL']=self.getGalactic()
        os.environ['EGAL']=self.getIsotropic()
        print "$GAL: ",os.getenv("GAL")
        print "$EGAL: ",os.getenv("EGAL")
        tmp_xml = NamedTemporaryFile(prefix=os.getenv("TEMPDIR")+"/").name
        self.oldXml = self.Xml
        xo = expandEnvVarsInXml(self.Xml)
        f = open(tmp_xml,'w')
        f.write(xo.toprettyxml())
        self.Xml  = tmp_xml
        print '*INFO* use temporary xml %s'%self.Xml

    def cleanUp(self):
        if self.oldXml:
            os.remove(self.Xml)
        self.Xml = self.oldXml
        print '*INFO* restored xml %s'%self.Xml

    def __call__(self):

        return self.__dict__
    def getJ(self):
        return float(self.J)
    def getFits(self):
        return self.fits
    def getIsotropic(self):
        if self.egal is None:
            return None
        else:
            return os.path.join(self.DiffDir,self.egal)
    def getGalactic(self):
        if self.gal is None:
            return None
        else:
            #return self.gal #os.path.join(self.DiffDir,self.gal)
            return os.path.join(self.DiffDir,self.gal)
