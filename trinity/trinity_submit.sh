#!/bin/env bash

. $LAB/scripts/dotbashrc
MYHIST=$HOME/.readline_inputs
history -r $MYHIST
history -a $MYHIST
history -n $MYHIST

echo "This script will attempt to submit a trinity denovo assembly."
echo "Ideally, it is possible to change everything with readline,"
echo "thus you should be able to hit up arrow to look at old arguments..."
## Set initial arguments for pbs, trinity
export PBS=" -V -S /cbcb/lab/nelsayed/local/bin/bash -q large -l mem=120G,walltime=144:00:00 -N trinity -jeo -e ${WD}/trinity.stdout "
#export PBS=" -V -S /cbcb/lab/nelsayed/local/bin/bash -q high_throughput -N trinity -jeo -e ${WD}/trinity.stdout "
export TRINITY=" --seqType fq --min_contig_length 600 --JM 100G  --CPU $(cpus) --no_cleanup "
export MYWD=$WD
export R1=" --left ${MYWD}/"
export R2=" --right ${MYWD}/"
export LOG="./trinity.out"

echo "$PBS"
echo "If necessary, change the PBS_ARGS now:"
read -e -i "${PBS}"
export MY_PBS_ARGS=${REPLY%%*( )}
history -s "$PBS_ARGS"

echo "If necessary, change the Trinity arguments now:"
read -e -i "${TRINITY}"
export TRINITY_ARGS=${REPLY%%*( )}
history -s "$TRINITY_ARGS"

echo "Where are the read1 reads?:"
read -e -i "${R1}"
export R1=${REPLY%%*( )}
history -s "$R1"

echo "Where are the read2 reads?:"
read -e -i "${R2}"
export R2=${REPLY%%*( )}
history -s "$R2"

echo "Redirect STDOUT and STDERR to:"
read -e -i "${LOG}"
export LOG=${REPLY%%*( )}
history -s "$LOG"

history -w $MYHIST
TRINITY="$PREFIX/trinity/Trinity.pl "
export MYCMD="$TRINITY $TRINITY_ARGS $R1 $R2 2>$LOG 1>&2 "

echo "Last chance to stop the madness:"
echo "PBS_ARGS: $MY_PBS_ARGS"
echo "CMD: $MYCMD"

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

echo "JOB completed in $SECONDS"
