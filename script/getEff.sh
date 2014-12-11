#!/bin/bash
if [ ! $# == 2 ]
then
	echo "$0 [startNo] [endNo]"
	exit 1
fi
startNo=$1
endNo=$2
../program/getEff $startNo $endNo

prePWD=`pwd`
