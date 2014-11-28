'''
Created on Nov 28, 2014

@author: zimmer
'''
import yaml, sys, glob, os

ifile = sys.argv[1]
print 'INFO: processing path %s'%ifile
files = glob.glob(ifile)
print 'INFO: found files: {}'.format(files)

for fi in files:
    fo = open(os.path.abspath(fi),'rb')
    d = yaml.load(fo)
    print os.path.basename(fi), d.keys()
