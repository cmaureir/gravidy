#!/usr/bin/env python
# Script to parse the output files from GraviDy

import argparse
import os
from numpy import ceil

class Particles:
    def __init__(self):
        self.time  = []
        self.index = []
        self.mass  = []
        self.rx    = []
        self.ry    = []
        self.rz    = []
        self.vx    = []
        self.vy    = []
        self.vz    = []
        self.ax    = []
        self.ay    = []
        self.az    = []
        self.jx    = []
        self.jy    = []
        self.jz    = []
        self.dt    = []
    def add_elements(self,attributes):

        self.time.append(attributes[0])
        self.index.append(attributes[1])
        self.mass.append(attributes[2])
        self.rx.append(attributes[3])
        self.ry.append(attributes[4])
        self.rz.append(attributes[5])
        self.vx.append(attributes[6])
        self.vy.append(attributes[7])
        self.vz.append(attributes[8])
        self.ax.append(attributes[9])
        self.ay.append(attributes[10])
        self.az.append(attributes[11])
        self.jx.append(attributes[12])
        self.jy.append(attributes[13])
        self.jz.append(attributes[14])
        self.dt.append(attributes[15])



class OutputFile:
    def __init__(self, fname):
        self.fname = fname
        self.n     = 0
        self.e2    = 0
        self.eta   = 0
        self.t     = 0
        self.t_rh  = 0
        self.t_cr  = 0

        self.info_lines       = []
        self.log_lines        = []
        self.radii_lines      = []
        self.all_lines        = []
        self.all_lines_split  = []
        self.particles        = []
        self.radii            = []
        self.radii_percentage = []
        self.times            = []
        self.iterations       = []
        self.nsteps           = []
        self.energy           = []
        self.rel_energy       = []
        self.cum_energy       = []
        self.clock_time       = []
        self.gflops           = []
        self.particles        = Particles()

    def read_output_file(self):
        with open(self.fname) as fn:
            for line in fn:
                if line.startswith('#'):
                    self.info_lines.append(line.strip())
                elif line.startswith('00'):
                    self.log_lines.append(line.strip())
                elif line.startswith('01'):
                    self.radii_lines.append(line.strip())
                else:
                    self.all_lines.append(line.strip())
        fn.close()

    def generate_directory(self):
        only_file = self.fname.split("/")[-1]
        os.mkdir(only_file)

        info = open(only_file+"/info", "w")
        for i in self.info_lines:
            info.write(i)
            info.write("\n")
        info.close()

        log = open(only_file+"/log", "w")
        for i in self.log_lines:
            log.write(i[3:])
            log.write("\n")
        log.close()

        radii = open(only_file+"/radii", "w")
        for i in self.radii_lines:
            radii.write(i[3:])
            radii.write("\n")
        radii.close()

        # All for every integration time
        for t in range(int(ceil(self.t))+1):
            lower = t*self.n
            higher = (t+1)*self.n
            self.all_lines_split.append(self.all_lines[lower:higher])
            all = open(only_file+"/all.t"+str(t), "w")
            for i in self.all_lines_split[t]:
                all.write(i)
                all.write("\n")
            all.close()



    def parse_lines(self):
        ########################################################
        # File information: self.info_lines
        ########################################################
        tmp_lines = self.info_lines[:-1] # removing header of the general log
        self.n          =   int(tmp_lines[0].split()[2])
        self.e2         = float(tmp_lines[1].split()[2])
        self.eta        = float(tmp_lines[2].split()[2])
        self.t          = float(tmp_lines[3].split()[2])
        self.t_rh       = float(tmp_lines[4].split()[2])
        self.t_cr       = float(tmp_lines[5].split()[2])

        ########################################################
        # General log: self.log_lines
        ########################################################
        for i in self.log_lines:
            i = i.split()[1:]
            self.times.append(i[0])
            self.iterations.append(i[1])
            self.nsteps.append(i[2])
            self.energy.append(i[3])
            self.rel_energy.append(i[4])
            self.cum_energy.append(i[5])
            self.clock_time.append(i[6])
            self.gflops.append(i[7])

        ########################################################
        # Lagrange radii: self.radii_lines
        ########################################################
        self.radii = [[] for i in range(len(self.radii_lines[0].split()) - 2)]
        ratio = 1.0/len(self.radii)
        self.radii_percentage = [i*ratio for i in range(1,len(self.radii))]
        for i in self.radii_lines:
            line = i.split()[2:]
            for j in range(len(line)):
                self.radii[j].append(line[j])

        ########################################################
        # Particles information: self.all_lines
        ########################################################
        for i in self.all_lines:
            self.particles.add_elements([float(j) for j in i.split()])

######################################## Parsing arguments
def parse_args():
    desc='Parsing GraviDy output file'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input', type=str, help='Input file')
    parser.add_argument('-g', '--generate', help='Create a directory with different\
    files for every type of logs and integration time')
    args = parser.parse_args()
    if args.input == None:
        parser.print_help()
        sys.exit(0)

    options = dict()
    options['input_file'] = args.input
    options['generate'] = args.generate
    return options

######################################## Main
if __name__ == "__main__":

    ops = parse_args()
    fname = ops['input_file']
    out = OutputFile(fname)
    out.read_output_file()
    out.parse_lines()
    if ops['generate'] != None:
        out.generate_directory()

