#!/bin/bash
build(){
   runNo=$1
   startNo=$2
   endNo=$3
   echo "building $runNo $startNo $endNo"
   /data3/work/eventbuilder/eventbuilder $runNo $startNo $endNo > /data3/work/data-built/$runNo.log 2> /data3/work/data-built/$runNo.err
   echo "translating $runNo $startNo $endNo"
   /data3/work/binary2root/binary2rootC $runNo > /data3/work/root/$runNo.log 2> /data3/work/root/$runNo.err
   echo "  Finished!"
}

#build 104 0 12
#build 105 0 11
#build 106 0 6
#build 107 0 5
#build 108 0 5
#build 109 0 6
#build 116 0 3
for (( i=216; i<=226; i++ ))
do
   n=`ls -l /data3/work/data/run_000"$i"_*.dat|wc -l`
   ((n--))
   build $i 0 $n
done
