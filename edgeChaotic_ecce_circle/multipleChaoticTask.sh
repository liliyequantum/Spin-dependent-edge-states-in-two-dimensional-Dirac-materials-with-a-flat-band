#!/bin/bash
for i in $(seq 23.0 0.5 27.5)
do
   matlab nohup -nodisplay -r "sigmaDiffChaotic($i);exit" > ./log/log_sigmaDiffmu_$i.txt &
done