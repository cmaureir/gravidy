#!/usr/bin/env python

from subprocess import call

ini = 0.001
end = 0.1
step = (end-ini)/8
print([ini+step*i for i in range(0,9)])

# Modifying ETA_N value
#for i in etas:
#    print('Trying eta:', i)
#    cmd = "sed -i 's/ETA_N\ 0.*/ETA_N\ "+i+"/g' include/common.hpp"
#    response = call(cmd, shell=True)
#    if response != 0:
#        print('Error: sed ', i)
#    else:
#        print('[PASS] : sed')
#
#    cmd = "make distclean all"
#    response = call(cmd, shell=True)
#    if response != 0:
#        print('Error: make ', i)
#    else:
#        print('[PASS] : make')
#
#    cmd = "./gravidy -i input/04-nbody-p1024_m1.in -t 10 >> etas_e"
#    response = call(cmd, shell=True)
#    if response != 0:
#        print('Error: gravidy ', i)
#    else:
#        print('[PASS] : execution')
