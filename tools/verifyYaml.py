'''
Created on Nov 28, 2014

@author: zimmer
'''
import yaml, sys, glob, os, numpy as np

from optparse import OptionParser
usage = "Usage: %prog  [options] input"
description = "python script"
parser = OptionParser(usage=usage,description=description)
parser.add_option("--debug",dest='debug', action='store_true', default = False,
                  help = "")
parser.add_option("--mass",dest="mass",default=None,help="if set, specify the masses")
(opts, args) = parser.parse_args()

masses = [5,10,20,50,100,200,500,1000.,2000.]
if not opts.mass is None:
    if ',' in opts.fstate:
        masses = [float(m) for m in opts.fstate.split(",")]
    else:
        masses = [float(opts.mass)]
masses = sorted(masses)
print 'INFO: checking for masses {}'.format(masses)
ifile = sys.argv[1]
print 'INFO: processing path %s'%ifile
files = glob.glob(ifile)
print 'INFO: found %i files'%len(files)

for fi in files:
    IS_OK = True
    fo = open(os.path.abspath(fi),'rb')
    d = yaml.load(fo)
    mkeys = sorted(d.keys())
    not_in_mkeys = [m for m in masses if not m in mkeys]
    if len(not_in_mkeys):
        IS_OK= False
        print os.path.basename(fi), not_in_mkeys
