#!/bin/env bash
. $LAB/scripts/dotbashrc
MYHIST=$PWD/fasta_inputs.txt
history -r $MYHIST
history -a $MYHIST
history -n $MYHIST

QUERY="${PWD}/"
LIBRARY="${PWD}/"

echo "$PBS"
echo "If necessary, change the query library now:"
read -e -i "${QUERY}"
export QUERY=${REPLY%%*( )}
history -s "$QUERY"


echo "If necessary, change the search library now:"
read -e -i "${LIBRARY}"
export LIBRARY=${REPLY%%*( )}
history -s "$LIBRARY"

export PBS=" -V -S /cbcb/lab/nelsayed/local/bin/bash -q long -l walltime=168:00:00 -N fasta -jeo -e ${WD}/fasta.stdout "
echo "If necessary, change the PBS_ARGS now:"
read -e -i "${PBS}"
export MY_PBS_ARGS=${REPLY%%*( )}
history -s "$MY_PBS_ARGS"

export MYCMD="ggsearch36 -b 1 -d 1 $QUERY $LIBRARY 2>>${WD}/ggsearch.err 1>>${WD}/ggsearch.txt"

cat <<"EOF" | qsub $MY_PBS_ARGS -
#!/bin/env bash

. ${LAB}/scripts/dotbashrc

cd $MYWD
start_log
echo "About to run:"
echo "$MYCMD"
sleep 5
eval $MYCMD
end_log
EOF

