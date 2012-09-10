#!/bin/bash
for i in $(/bin/ls -1 input/ | grep -v py);
do
    echo $i >> CPU_E
    echo $i >> GPU_E
     ./gravidy -i input/$i -t 1 -r cpu >> CPU_E
     ./gravidy -i input/$i -t 1 -r gpu >> GPU_E
done
