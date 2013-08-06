#!/usr/bin/env python

from imports import *

def read_gravidy_output(fname,N,S):
    summary_array = []
    snapshot_array = []
    f = open(fname, 'r')

    line_counter = 0
    # Skipping first line
    f.readline()

    # The output file will contain a summary line with system properties
    # for each snapshop (S)
    # and then N lines with the information of all the particles of the ststem
    for line in f.readlines():
        line = map(float, filter(None, line.strip().split(' ')))

        if not line_counter%(N+1):
            summary_array.append(line)
        else:
            snapshot_array.append(line)
        line_counter += 1

    summary_array = np.array(summary_array)
    snapshot_array = np.array(snapshot_array)

    return summary_array, snapshot_array


