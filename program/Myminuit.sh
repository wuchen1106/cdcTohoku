#!/bin/bash

for (( i=100; i<266; i++ ))
do
	if [ ! -f ../root/d_"$i".root ] 
#	if [ ! -f ../good/"$i" ] 
	then
		continue
	fi
	s=`du -s ../root/d_"$i".root|sed 's/\(\w*\) *.*/\1/'`
	echo "$i: $s KB"
	if [ $s -gt 8 ]
	then
		echo "	larger than 8 KB"
		for (( j=1; j>0; j++ ))
		do
			sleep 2
			n=`ps -ef | grep getD | wc -l`
			if [ $n -lt 2 ]
			then
				n=`ps -ef | grep myminuit | wc -l`
				echo "		$j: $n"
				if [ $n -lt 6 ]
				then
					echo "submit $i"
					NMIN=50
					if [ $i -lt 215 -a $i -gt 209 ]; then NMIN=5;
					elif [ $i -lt 252 -a $i -gt 249 ]; then NMIN=100;
					elif [ $i -eq 236 ]; then NMIN=24;
					elif [ $i -eq 226 ]; then NMIN=251;
					elif [ $i -eq 221 ]; then NMIN=100;
					elif [ $i -eq 199 ]; then NMIN=80;
					elif [ $i -eq 198 ]; then NMIN=200;
					elif [ $i -lt 198 -a $i -gt 194 ]; then NMIN=5;
					elif [ $i -lt 179 -a $i -gt 172 ]; then NMIN=5;
					elif [ $i -lt 179 -a $i -gt 172 ]; then NMIN=5;
					elif [ $i -lt 157 -a $i -gt 154 ]; then NMIN=10;
					elif [ $i -lt 155 -a $i -gt 148 ]; then NMIN=5;
					elif [ $i -eq 146 ]; then NMIN=14;
					elif [ $i -eq 132 ]; then NMIN=14;
					elif [ $i -eq 131 ]; then NMIN=7;
					elif [ $i -eq 130 ]; then NMIN=3;
					elif [ $i -eq 129 ]; then NMIN=10;
					elif [ $i -eq 128 ]; then NMIN=10;
					elif [ $i -eq 127 ]; then NMIN=3;
					elif [ $i -eq 126 ]; then NMIN=3;
					elif [ $i -eq 125 ]; then NMIN=3;
					elif [ $i -eq 124 ]; then NMIN=85;
					elif [ $i -eq 123 ]; then NMIN=45;
					elif [ $i -lt 123 -a $i -gt 102 ]; then NMIN=5; fi
					echo ./myminuit $i 0309 0 0 $NMIN 7 ">>" logs/log.myminuit.$i 2">>" logs/err.myminuit.$i "&"
					nohup ./myminuit $i 0309 0 0 $NMIN 7 >> logs/log.myminuit.$i 2>> logs/err.myminuit.$i &
					break
				fi
			fi
		done
	fi
done

#i=237; nohup ./myminuit $i 0309 0 0 50 7 >> logs/log.myminuit.$i 2>> logs/err.myminuit.$i &
