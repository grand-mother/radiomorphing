#!/usr/bin/env python
import os
from os.path import split, join, realpath
import sys

#import time
import subprocess
import shlex


for j in range(5):
    if j>1:
     try:
        print p.communicate()
     except NameError:
        continue
                

                cmd1='ishell -c "get grand/sim/flat-150x67km2/taus/events.%s.json flat-150x67km2/taus/"' %(str(j) )
                try:
                    p=subprocess.Popen(shlex.split(cmd1))
                except OSError:
                    continue

                

                
try:                
    print p.communicate()
except NameError:
    pass   

print "Job done"