#!/usr/bin/env python
import sys,os

nargs = len(sys.argv)
print sys.argv[0], nargs
if nargs!=2:
    print "Forgot to specify a file name"
    exit(0)
else:
    filename = sys.argv[1]
    
os.system("root -l 'TMVAGui.C(\""+filename+"\")'")
    
                    

