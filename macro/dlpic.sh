#!/bin/bash
for (( i=220; i<=223; i++ ))
do
#	DLkl_daq /data3/work/user/chen/pdf "$i"_eff1.pdf "$i"_eff2.pdf "$i"_eff1.png "$i"_eff2.png 1000
#	mv *.png *.pdf pics
	DLkl_daq /data3/work/user/chen/root run_000"$i"_built.root 1000
done
