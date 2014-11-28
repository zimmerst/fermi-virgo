'''
Created on Nov 28, 2014

@author: zimmer
'''
import yaml, sys, glob, os

ifile = glob.glob(sys.argv[1])

for fi in ifile:
    fo = open(os.path.abspath(fi),'rb')
    d = yaml.load(fo)
    print os.path.basename(fi), d.keys()
