#!/bin/bash
# $1 is the process name
for n in {20..26..2}
do
	for p in {0..7}
	do
		for b in {12..16}
		do
			`./$1 $n $n $p $b > out-$n-$p-$b.txt`
		done
	done
done
