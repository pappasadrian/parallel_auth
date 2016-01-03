#!/bin/bash
#$1 is how many nodes (power of two)
#$2 is how many processes for the actual program (power of two), 0-7
#ppn will be evaluated automatically
p=$2
q=20
c=20
b=12
nodes=$1
ppn=p-nodes
two=2
for q in {20..26..2}
do
	for c in {20..26..2}
	do
		for b in {12..16..2}
		do
			
			nodespowered=$((2**nodes))
			ppnpowered=$((2**ppn))
			`sed "5s/.*/#PBS ­-l nodes=$nodespowered:ppn=$ppnpowered/" run.sh > run-$nodes-$ppn.sh`
			`sed "13s/q c p b/$q $c $p $b/g" run-$nodes-$ppn.sh > run-$nodes-$ppn-$q-$c-$p-$b.sh`
			chmod +x run-$nodes-$ppn-$q-$c-$p-$b.sh
			#`./run-$nodes-$ppn-$q-$c-$p-$b.sh`
		done
	done
done

