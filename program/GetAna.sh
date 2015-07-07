#!/bin/bash

for (( i=100; i<=266; i++ ))
do
  if [ -f ../root/ana.$i.0309.0.root ]
  then
    result=`./getAna $i 1 0309 0 0 7`
    echo $i" "$result
  else
    echo $i 0 0 0 0 0 0  0 0 0 0 0 0  0 0 0 0 0 0  0 0 0 0 0 0  0 0 0 0 0 0  0 0 0 0 0 0  0 0 0 0 0 0  0 0 0 0 0 0 0 0
  fi
done
