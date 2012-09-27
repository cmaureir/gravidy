#!/usr/bin/env python
from subprocess import call

files = ['04-nbody-p1024_m1.in', \
        '05-nbody-p2048_m1.in',  \
        '06-nbody-p4096_m1.in',  \
        '07-nbody-p8192_m1.in',  \
        '08-nbody-p16384_m1.in', \
        '09-nbody-p32768_m1.in']

for f in files:
    print f
    st_cpu = call("./gravidy -i input/"+f+" -t 1 -r cpu", shell=True)
    if st_cpu != 0:
        print "Error with CPU version"
#    st_gpu = call("./gravidy -i input/"+f+" -t 1 -r gpu", shell=True)
#    if st_gpu != 0:
#        print "Error with CPU version"
