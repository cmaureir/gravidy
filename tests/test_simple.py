#!/usr/bin/env python
from subprocess import call

files = []

files.append('input/01-nbody-p128_m1.in')
files.append('input/02-nbody-p256_m1.in')
files.append('input/03-nbody-p512_m1.in')
files.append('input/04-nbody-p1024_m1.in')
files.append('input/05-nbody-p2048_m1.in')
files.append('input/06-nbody-p4096_m1.in')
files.append('input/07-nbody-p8192_m1.in')
files.append('input/08-nbody-p16384_m1.in')

for f in files:
    print f
    st_cpu = call("time ./gravidy -i "+f+" -t 1", shell=True)
    if st_cpu != 0:
        print "Error with CPU version"
