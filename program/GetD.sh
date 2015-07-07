#!/bin/bash

#for (( i=100; i<266; i++ ))
#do
#	s=`du -s ../root/run_000"$i"_built.root|sed 's/\(\w*\) *.*/\1/'`
#	echo "$i: $s KB"
#	if [ $s -gt 8 ]
#	then
#		echo "	larger than 8 KB"
#		for (( j=1; j>0; j++ ))
#		do
#			sleep 2
#			n=`ps -ef | grep getD | wc -l`
#			echo "		$j: $n"
#			if [ $n -lt 7 ]
#			then
#				echo "submit $i"
#				nohup ./getD $i &
#				break
#			fi
#		done
#	fi
#done

#nohup ./getD 234 &
#nohup ./getD 235 &
#nohup ./getD 165 &
#nohup ./getD 168 &
#nohup ./getD 169 &
#nohup ./getD 172 &
#nohup ./getD 204 &

#nohup ./getD 149 &
#nohup ./getD 168 &

nohup ./getD 250 &
