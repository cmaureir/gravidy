# Require the amuse code http://amusecode.org/
from amuse.support import data
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.plummer import new_plummer_sphere
import numpy as np
import argparse
import sys

# Parsing options
desc='Generate nbody initial conditions using the Plummer Model.'
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-p', '--particles', type=int, help='Number of particles')
args = parser.parse_args()
if args.particles == None:
    parser.print_help()
    sys.exit(0)

#1 nbu in m/s
msscale=65.58
mtot=1.0

# Using amuse to generate the system
convert_nbody = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
plummer = new_plummer_sphere(args.particles, convert_nbody)

# Opening the file
name = "nbody-p%d" % (args.particles)
f_out = open(name, "w")

# Shortcut for the variables
m = plummer.mass.value_in(units.MSun)*mtot
r = plummer.position.value_in(units.parsec)
v = plummer.velocity.value_in(units.ms)/msscale*np.sqrt(mtot)

# Writing on the file
for i in range(args.particles):
    s = "%0.15f %0.15f %0.15f %0.15f %0.15f %0.15f %0.15f\n" % (m[i],r[i][0],r[i][1],r[i][2],v[i][0],v[i][1],v[i][2])
    f_out.write(s)
