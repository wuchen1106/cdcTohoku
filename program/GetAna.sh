#!/bin/bash

for (( i=0; i<270; i++ ))
do
  if [ -f ../root/ana.$i.0309.0.root ]
  then
    result=`./getAna $i 1 0309 0 100 7`
    echo $i" "$result
  fi
done
