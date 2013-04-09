#!/usr/bin/env python

import re

def clean(line):
    line = re.sub(r'[^e0-9-.\ ]+','',line)
    return line

f = open('1k.out')
for line in f.readlines():
    print(clean(line))
