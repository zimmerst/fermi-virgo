# here we collect all relevant xml tools
import os
import xml.dom.minidom
import StringIO
# string operations

# pretty much obsolete
class pickle_dict(object):
    def __init__(self,tmpfile):
        self.tmpfile = tmpfile
        self.xmlfile = None
        self.dicts = []

    def __call__(self):
        kstring = '<?xml version="1.0" ?>\n<pickle/>'
        # add environment block
        self.add_dict(os.environ,name='environ')
        self.xmlfile = xml.dom.minidom.parse(StringIO.StringIO(kstring))
        for d in self.dicts:
            self.xmlfile.childNodes[-1].appendChild(d)
        return self.xmlfile
            
    def add_dict(self,DICT,name=None):
        # adds a pickled version of the dict
        dict_name = str(name)
        if name is None:
            dict_name = str(DICT)
        kstring = '<?xml version="1.0" ?>\n<pickle/>'
        xmlf = xml.dom.minidom.parse(StringIO.StringIO(kstring))
        el = xmlf.createElement("dict")
        el.setAttribute("name",str(dict_name))
        for key in DICT.keys():
            dict_el = xmlf.createElement("var")
            dict_el.setAttribute("key",str(key))
            dict_el.setAttribute("value",str(DICT[key]))
            el.appendChild(dict_el)
        self.dicts.append(el)

def xcomments(xmlfile):
    # given an xmlfile, remove comments, replace by xcomments
    xmlf = xml.dom.minidom.parse(xmlfile)
    for p in xmlf.getElementsByTagName("parameter"):
        comment = str(p.getAttribute("comment"))
        if len(comment)!=0:
            p.removeAttribute("comment")
            p.setAttribute("xcomment",comment)
    f = open(xmlfile,'w')
    f.write(xmlf.toxml())
    f.close()

def str_contains(string,sign):
    ''' returns True if sign is contained within string, False otherwise'''
    ps = string.partition(sign);
    if (len(ps[1])>0):
        return True;
    else:
        return False;
# other readouts

def ReadConfigXML(fname):
    doc = xml.dom.minidom.parse(fname);
    parname = [];
    parval = [];
    for par in doc.getElementsByTagName("parameter"):
        parname.append(str(par.getAttribute("name")));
        parval.append(str(par.getAttribute("value")));
    pars = dict(zip(parname,parval));
    # and here we also set the variables from the envVar block
    for var in doc.getElementsByTagName("var"):
        varname = str(var.getAttribute("name"))
        value = str(var.firstChild.data)
        os.environ[varname]=value        
    return pars;

# prettify newly written xmls

def prettifyXML(xmlfile):
    # to fix a bug in xml.dom.minidom.toprettyxml();
    # loads the xml file and removes blank lines
    import copy;
    new_lines  = [];
    f = open(xmlfile,'r');
    lines = f.readlines();
    for line in lines:
        # test this line:
        # let's try to reduce the line length...
        c_string = copy.deepcopy(line);
        c_string.replace("\n","");
        for i in range(0,10):
            new_str = c_string.partition("\t");
            if len(new_str[1])>0:
                c_string = new_str[2];
                #print 'step i',i,' string length',len(c_string);
        if len(c_string)>1:
            new_lines.append(line);
    f.close();
    f = open(xmlfile,'w');
    for l in new_lines:
        f.write(l);
    f.close();




# source model xml tools (for use with science tools)

def GetCombinedXML(ListOfXMLs):
    ''' reads in all xmls that are provided in a list and merges them together, removing duplicates, returns the combined xml in proper markup '''
    firstLine='<?xml version="1.0" ?> \n'+'<source_library title="source library">';
    lastLine='</source_library>';   
    EndXML = firstLine+'\n';
    source_names = [];
    source_xmls = [];
    included = []
    for xmlf in ListOfXMLs:
        # loop over all xml files
        print 'xml file now',xmlf;
        doc = xml.dom.minidom.parse(xmlf);
        for source in doc.getElementsByTagName("source"):
            source_names.append(source.getAttribute("name"));
            source_xmls.append(source.toxml()+'\n')
            included.append(True);
    # have everything together
    # loop over it again to find duplicates
    i = 0; j = 0;
    for i in range(0,len(source_names)):
        for j in range(i+1,len(source_names)):
            if (source_names[i]==source_names[j]):
                included[j] = False;
    final_list = [];
    for i in range(0,len(source_names)):
        if (included[i]):
            final_list.append(source_xmls[i]);
    # loop over it to append it to the end xml
    for l in final_list:
        EndXML+=l+'\n';
    EndXML+= lastLine;
    return EndXML;

def GetSourceIDs(xmlfile):
    '''returns a dict containing the name of the sources in a given xml along with its ID to use for likelihood, it sorts the names alphabetically'''
    names = [];
    IDsUnsorted=[];
    doc = xml.dom.minidom.parse(xmlfile);
    i = 0;
    for source in doc.getElementsByTagName("source"):
        names.append(str(source.getAttribute("name")));
        IDsUnsorted.append(i);
        i+=1;
    d = dict(zip(names,IDsUnsorted));
    keys = d.keys();
    keys.sort();
    IDsSorted = [];
    i=0;
    for i in range(0,len(keys)):
        IDsSorted.append(i);
        #i+=1;
    D = dict(zip(keys,IDsSorted));
    return D;

def GetParameterIDForLikelihood(xmlfile,name,dbg=False):
    ''' can be used for the individual likelihood analysis, returns the ID based on the number of free parameters in the xml in alphabetical ordering '''
    sources = [];
    doc = xml.dom.minidom.parse(xmlfile);
    for source in doc.getElementsByTagName("source"):
        kFree = False;
        for par in source.getElementsByTagName("parameter"):
            if par.getAttribute("free")=="1":
                kFree = True;
        if kFree:
            sources.append(source.getAttribute("name"));
    # okay now we have the sources
    sources_sorted = sorted(sources);
    if dbg:
        print '*SOURCES SORTED*',sources_sorted
    parameters = [];
    # loop again through and count parameters and names:
    for src in sources_sorted:
        for source in doc.getElementsByTagName("source"):
            if source.getAttribute("name")==src:
                for par in source.getElementsByTagName("parameter"):
                    if par.getAttribute("free")=="1":
                        parameters.append(par.getAttribute("name"));

    d = dict(zip(sources_sorted,parameters));
    if dbg:
        print '*DICT*',d
    keys = d.keys();
    keys.sort();
    i = 0;
    Id = [];
    for k in keys:
        #print d[k];
        i+=1;
        Id.append(i);
    dd = dict(zip(keys,Id));
    try:
        return dd[name];
    except KeyError:
        print 'FATAL: Cannot find ID, return 0'
        return 0;

def GetSourceNameFromXML(xmlfile,ID):
    ''' returns the source name and the parameter name for the ID in the xmlfile, can be used to test GetParameterIDForLikelihood '''
    names = [];
    pars = [];
    doc = xml.dom.minidom.parse(xmlfile);
    for source in doc.getElementsByTagName("source"):
        names.append(str(source.getAttribute("name")));
        parSet = [];
        for parameter in source.getElementsByTagName("parameter"):
            if (parameter.getAttribute("free")=='1'):
                parname = parameter.getAttribute(str("name"));
                parSet.append(parname);
        pars.append(parSet);
    d = dict(zip(names,pars)); # contains now the names along with the name of the parameters
    keys = d.keys();
    keys.sort();
    parameterName = [];
    parameterID = []
    i = 0;
    sourceName = [];
    for k in keys:
        #        print k;
        #        print d[k];
        for el in d[k]:
            sourceName.append(k);
            parameterName.append(el);
            parameterID.append(i);
            i+=1;
    FoundMatch = False;
    for i in range(0,len(sourceName)):
        if (i == ID):
            print '*INFO* xml:',xmlfile,'ID: ',ID,'found match, source: ',sourceName[i],' parameter Name: ',parameterName[i];
            FoundMatch = True;
    if not (FoundMatch):
        print 'Cannot find matching ID'

    return FoundMatch;
    
    
    

def ReplaceSource(fname,sID,dpars):
    keys = dpars.keys();
    doc = xml.dom.minidom.parse(fname);
    source = doc.getElementsByTagName("source")[sID];
    source.setAttribute("name",str(dpars['name']));
    pars = source.getElementsByTagName("parameter");
    for p in pars:
        for d in keys:
            if (d==p.getAttribute("name")):
                p.setAttribute("value",str(dpars[d]));
            elif (str_contains(d,"#")):
                # disentangle:
                newvar = d.partition("#");
                newpar = newvar[0];
                newval = newvar[2];
                if (p.getAttribute("name")==newpar):
                    #                    print 'found variable with hash, requires special treatment', d
                    p.setAttribute(newval,str(dpars[d]));
    return source.toxml()+'\n';

def MakeXML(fname,pars,d=0):
    ''' requires the name of the input xml (i.e. template) together with a dict containing 'source_id' and parameters to be replaced, if source_id is found to match with the sources in the template, the content of this source is then replaced. see xmltools.dummy() for details '''
    ''' enter a variable where you do not want to change only the value but some other xml aspect as variable#name and then the associated value in the dict '''
    # first and last line:
    firstLine='<?xml version="1.0" ?> \n'+'<source_library title="source library">';
    lastLine='</source_library>';
    # 1st step, loop over all sources
    modified = [];
    newxml = '';
    newxml+=firstLine+'\n';
    if (d>=2):
        print '*INFO* fname: ',fname;
    doc = xml.dom.minidom.parse(fname);
    i = 0;
    for s in doc.getElementsByTagName("source"):
        modified.append(False);
        i+=1;
    # okay, next step, check which ones we'll replace:
    for i in range(0,len(modified)):
        for j in pars:
            #print 'par: ',j,' in pars: ',pars;
            if (j['source_id']==i):
                modified[i]=True;
    #print modified;
    # finally, loop over sources again
    i = 0;
    for s in doc.getElementsByTagName("source"):
        if (modified[i]):
            for jp in pars:
                if (i == jp['source_id']):
                    newxml+=ReplaceSource(fname,i,jp);
        else:
            newxml+=s.toxml()+'\n';
        i+=1;
    newxml+=lastLine;
    return newxml;

def MakeXMLFromSeedfile(seedxml,roi,radius,suffix):
    ''' little tool that takes a seed xml, changes all parameters to being fixed and extracts all sources around a given radius '''
    from os import getenv;
    from catalogSorter3 import CalculateRadius;
    from toolbox import GetCoordinatesFromFitsFile
    print 'MakeXMLFromSeedFile'
    REFRA = float(roi.RA);
    REFDEC = float(roi.DEC);
    #print 'REF',REFRA,REFDEC;
    name = getenv('PWD')+'/'+roi.name+'_seeded.xml';
    # okay take some of the stuff...
    f = open(seedxml,'r');
    flines = f.readlines()
    firsts = flines[0];
    firstLine = '<?xml version="1.0" ?> \n'+firsts;
    lastLine = flines[-1]
    f.close();
    xmlstring = '';
    seed = xml.dom.minidom.parse(seedxml)
    RA = 0.; DEC=0.;
    for source in seed.getElementsByTagName("source"): # to get each source
        try:
            RA = float(source.getAttribute("RA"));
            DEC = float(source.getAttribute("DEC"));
        except ValueError:
            for model in source.getElementsByTagName("spatialModel"):
                RA,DEC = GetCoordinatesFromFitsFile(model.getAttribute("file"));

        #print 'Coordinates',RA,DEC;
        if (CalculateRadius(RA,REFRA,DEC,REFDEC)<=radius):
            old_name = str(source.getAttribute("name"));
            source.setAttribute("name",old_name+suffix); # to make it compatible with old code
            for par in source.getElementsByTagName("parameter"):
                #print par.getAttribute("name");
                if (par.getAttribute("name")=="Integral"):
                    par.setAttribute("free","1");
                else:
                    par.setAttribute("free","0");
                    
            xmlstring+=source.toxml()+'\n';
    f = open(name,'w');
    #print firstLine,lastLine;
    f.write(firstLine);
    f.write(xmlstring+'\n');
    f.write(lastLine);
    f.close();
    return name;

def UpdateDataXML(fname):
    # opens the e.g. cluster.xml and extracts the info xml from it - then open the info.xml and extract all objects and append to cluster.xml as object
    from toolbox import GetCoordinatesFromFitsFile;
    from catalogSorter3 import CalculateRadius;
    # get first and last line
    f = open(fname,'r');
    lines = f.readlines();
    f.close();
    del f;
    firstLine = lines[0]+lines[1];
    lastLine = lines[-1];
    xmlfile = xml.dom.minidom.parse(fname);
    xmlstring = '';
    for roi in xmlfile.getElementsByTagName("ROI"):
        infoxml_str = roi.getAttribute("xml");
        if infoxml_str == 'None':
            raise Exception("Error, no xml field specified in ",str(roi.getAttribute("name"))," occurred in file: ",fname);
        infoxml = xml.dom.minidom.parse(infoxml_str);
        roi_RA = roi.getAttribute("RA");
        roi_DEC = roi.getAttribute("DEC");
        for source in infoxml.getElementsByTagName("source"):
            name = source.getAttribute("name");
            for model in source.getElementsByTagName("spatialModel"):
                #print 'check spatial model',model;
                #d = model.attributes.keys();
                #for key in d:
                #    print 'key: ',key,' value: ',model.attributes[key].value;
                RA = model.getAttribute("RA");
                DEC = model.getAttribute("DEC");
                if (len(RA)==0 or len(DEC)==0):
                    #print 'try to extract coordinates from fitsfile'
                    for model in source.getElementsByTagName("spatialModel"):
                        fitsfile = str(model.getAttribute("file"));
                        #print '*** use: ',fitsfile,' **';
                        RA,DEC = GetCoordinatesFromFitsFile(fitsfile);
                        #print '**** test*** ';
                        #print 'coordinates: ',RA,DEC;
            minos_ID = GetParameterIDForLikelihood(infoxml_str,name)
            for par in source.getElementsByTagName("parameter"):
                if (par.getAttribute("name")=='Integral' or par.getAttribute("name")=='norm'):
                    norm = par.getAttribute("value");
                    norm_error = par.getAttribute("error");
                    norm_scale = par.getAttribute("scale");
            # okay now we have everything together to create a object-xml
            obj = xmlfile.createElement("object");
            obj.setAttribute("name",name);
            obj.setAttribute("RA",str(RA));
            obj.setAttribute("DEC",str(DEC));
            obj.setAttribute("norm",norm);
            obj.setAttribute("norm_error",norm_error);
            obj.setAttribute("norm_scale",norm_scale);
            obj.setAttribute("minos_ID",str(minos_ID));
            obj.setAttribute("distance_to_roi",str(CalculateRadius(float(roi_RA),float(RA),float(roi_DEC),float(DEC))));
            roi.appendChild(obj);
        xmlstring+=roi.toprettyxml();
    print xmlstring;
    f = open(fname,'w');
    f.write(firstLine);
    f.write(xmlstring+'\n');
    f.write(lastLine);
    f.close();
    del f;
            
def dummy():
    # dummy test
    print 'these parameters are supplied in a dict (here we give 2 sources)'
    dname1 = ['source_id','name','norm','sigmav','RA','DEC','norm#error'];
    dname2 = ['source_id','name','Normalization'];
    dvalue1 = [0,'blub',0,0,0,0,0];
    dvaule2 = [1,'BLUB',0];
    d1 = dict(zip(dname1,dvalue1));
    d2 = dict(zip(dname2,dvaule2));
    print d1,d2;
    d = [d1,d2];
    print 'combine them into a list:\n',d
    print 'now print out the output of xmltools.MakeXML(): ',MakeXML("template_dmcluster_p7v4.xml",d);
    return 1;

def dummy2():
    DM_source = ['source_id','name','norm','norm#error','RA','DEC'];
    DM_values = [0,'AWM7','1.8','0.03','43.1','-23.3'];
    
    DM = dict(zip(DM_source,DM_values));
    d = [DM];
    print 'combine them into list',d;
    print 'output of xmltools.MakeXML',MakeXML("template_dmcluster_p6v3.xml",d);
    return 1;


def expandEnvVarsInXml(xfile,ddict=None):
    ''' tries to expand environment variables in xml '''
    from logging import Logger
    import re
    #Log = Logger()
    x = xml.dom.minidom.parse(xfile)
    for source in x.getElementsByTagName("source"):
        atts = dict(zip(source.attributes.keys(),[att.value for att in source.attributes.values()]))
        spec = source.getElementsByTagName("spectrum")[0]
        atts_spec = dict(zip(spec.attributes.keys(),[att.value for att in spec.attributes.values()]))
        spat = source.getElementsByTagName("spatialModel")[0]
        atts_spat = dict(zip(spat.attributes.keys(),[att.value for att in spat.attributes.values()]))
        # NOW CHECK FOR ATTRIBUTES
        all_atts = [atts,atts_spec,atts_spat]
        for att_block in all_atts:
            for key in att_block:
                if "$" in att_block[key]:
                    #print '*FOUND BLOCK*'
                    # now try to identify var
                    results = re.findall(r'\$\(.+\)',att_block[key])
                    #print results
                    var = None
                    if len(results)!=0:
                        for result in results:
                            var = (result.replace("$(","")).replace(")","")
                            if os.getenv(var):
                                att_block[key]=att_block[key].replace(result,os.getenv(var))
                                #print att_block[key]
                            else:
                                print "found env var %s but not expanded"%var
                    else:
                        print "Could not interpret %s"%att_block[key]
        # now need to re-write info
        for key in atts:
            source.setAttribute(key,atts[key])
        for key in atts_spec:
            spec.setAttribute(key,atts_spec[key])
        for key in atts_spat:
            spat.setAttribute(key,atts_spat[key])
        # and now parameters
        for par in source.getElementsByTagName("parameter"):
            atts = dict(zip(par.attributes.keys(),[att.value for att in par.attributes.values()]))
            for key in atts:
                if "$" in atts[key]:
                    # now try to identify var
                    results = re.findall(r'\$\(.+\)',atts[key])
                    var = None
                    if len(results)!=0:
                        for result in results:
                            var = (result.replace("$(","")).replace(")","")
                            if os.getenv(var):
                                atts[key]=atts[key].replace(result,os.getenv(var))
                            else:
                                print "found env var %s but not expanded"%var
                    else:
                        print "Could not interpret %s"%atts[key]
    # THIS CODE IS SLOW!
    return x
