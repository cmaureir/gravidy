#!/bin/bash
for i in $(/bin/ls -1 input/ | grep -v py);
do
    echo $i >> $i.cpu.out
     ./gravidy -i input/$i -t 1 -r cpu &>> $i.cpu.out
    echo $i >> $i.gpu.out
     ./gravidy -i input/$i -t 1 -r gpu &>> $i.gpu.out
done
