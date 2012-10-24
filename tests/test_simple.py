#!/usr/bin/env python
from subprocess import call

files = ['04-nbody-p1024_m1.in', \
        'nbody-p256.in',  \
        'nbody-p128.in',  \
        '05-nbody-p2048_m1.in']

for f in files:
    print f
    st_cpu = call("./gravidy -i input/"+f+" -t 1 -r cpu &> output_"+f, shell=True)
    if st_cpu != 0:
        print "Error with CPU version"
#    st_gpu = call("./gravidy -i input/"+f+" -t 1 -r gpu", shell=True)
#    if st_gpu != 0:
#        print "Error with CPU version"
