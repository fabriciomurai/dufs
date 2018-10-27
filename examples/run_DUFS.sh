#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. >/dev/null && pwd )"
BIN_DIR=$DIR/bin             # BIN_PATH
EXAMPLES_DIR=$DIR/examples   # EXAMPLES_PATH

JUMP_COST=1             # cost of random jump relative to walker transition
NONUNIFORM_SAMPLING=0   # non-uniform sampling? (see options below)
RUNS=1000               # number of runs per simulation
NTHREADS=2              # number of threads

if [[ $# -eq 9 ]]; then
    path=$1
    USE_HYBRID=$2
    AVG_GIVEN=$3
    REDUCE_VAR=$4
    SEE_INCOMING=$5
    NSTEPS=$6
    EXTENSION=$7
    JUMP_WEIGHT=$8
    NOT_RECURSIVE=$9
else

echo "
ERROR! Incorrect number of parameters.

SYNOPSIS
	$0 <DATASET PATH> <USE_HYBRID> <AVG_GIVEN> <REDUCE_VAR> <SEE_INCOMING> <NSTEPS> <EXTENSION> <JUMP_WEIGHT> <NOT_RECURSIVE>

DESCRIPTION
	DATASET PATH is the relative path to a gzipped graph file.
	USE_HYBRID is a binary value indicating whether vertex samples should be used in addition to edge samples.
	When USE_HYBRID is 1:
	 * AVG_GIVEN is a binary value indicating whether the averages are known and can be used to compute normalization constants.
	 * REDUCE_VAR is a binary value indicating whether the variance reduction heuristic should be used.
	 * NOT_RECURSIVE is a binary value indicating whether the simplified version of the estimator should be used.
	When USE_HYBRID is 0:
	 * AVG_GIVEN, REDUCE_VAR and NOT_RECURSIVE are not applicable (but some arbitrary values must be provided anyway).
	 * If NSTEPS is Inf, this results in a single walker; this is the Single Walk with Jumps.
	SEE_INCOMING is a binary value indicating whether the edges incoming to a node are visible.
	NSTEPS is the average budget per walker (same as number of transitions when jump cost is 1).
	EXTENSION is (out | in | jnt) and indicates the metric of interest: out-degree, in-degree or joint degree distribution.
	JUMP_WEIGHT is a real number that controls the weight of the random jump relative to the weight of a regular edge.

EXAMPLES
    # non-hybrid, invisible in-edges, 10 steps per walker, out-degree, jump weight is 100
    $0 datasets/wiki-Vote-lcc.txt.gz 0 X X 0 10 out 100 X

    # hybrid, invisible in-edges, 10 steps per walker, out-degree, jump weight is 100, recursive
    $0 datasets/wiki-Vote-lcc.txt.gz 1 0 0 0 10 out 100 0

    # hybrid, reduce variance, invisible in-edges, 100 steps per walker, out-degree, jump weight is 10, recursive
    $0 datasets/wiki-Vote-lcc.txt.gz 1 0 1 0 100 out 10 0

    # hybrid, average given, visible in-edges, 1000 steps per walker, jnt-degree, jump weight is 100, non-recursive
    $0 datasets/wiki-Vote-lcc.txt.gz 1 1 0 1 1000 jnt 100 1
"

    exit 1
fi


RESULTS_DIR=results_${EXTENSION}


RW_VARIANT=DUFS_main
if [[ $USE_HYBRID -eq 1 ]]; then
    USE_HYBRID="--use_hybrid"
    if [[ $REDUCE_VAR -eq 1 ]]; then
        REDUCE_VAR="--reduce_var"
        SUFFIX="dufs-mle-var_red"
    else
        REDUCE_VAR=""
        SUFFIX="dufs-mle"
    fi

    if [[ $AVG_GIVEN -eq 1 ]]; then
        AVG_GIVEN="--avg_given"
        SUFFIX=${SUFFIX}-avg
    else
        AVG_GIVEN=""
    fi

    if [[ $NOT_RECURSIVE -eq 1 ]]; then
        NOT_RECURSIVE="--not_recursive"
        SUFFIX=${SUFFIX}-not_rec
    else
        NOT_RECURSIVE=""
    fi
else
    SUFFIX="dufs-e"
    if [[ "$NSTEPS" == "Inf" ]]; then
        SUFFIX=rwsw
    fi
    USE_HYBRID=""
    AVG_GIVEN=""
    REDUCE_VAR=""

fi

OUTDIR=$RESULTS_DIR/${SUFFIX}

case $NONUNIFORM_SAMPLING in
    1)
        OUTDIR=${OUTDIR}.propto
        VS_FLAG="--propto_degree"
        ;;
    2)
        OUTDIR=${OUTDIR}.inv_propto
        VS_FLAG="--inv_propto_degree"
        ;;
    *)
        VS_FLAG=""
        ;;
esac


if [[ ! -d $OUTDIR ]]; then
    mkdir -p $OUTDIR
fi


file=`basename $path`
N=`$EXAMPLES_DIR/get_nodes.sh $file`
for BUDGET in 0.1; do
    if [[ "$NSTEPS" == "Inf" ]]; then
        WALKERS=1
    else
        WALKERS=`echo "($BUDGET*$N)/($NSTEPS+$JUMP_COST)" | bc`
    fi
    JUMP_STR=""
    if [[ "$JUMP_WEIGHT" != "0" ]]; then
        JUMP_STR="--jump_weight $JUMP_WEIGHT"
    fi 

    $BIN_DIR/$RW_VARIANT \
    --filename $path \
    --budget $BUDGET \
    --walkers $WALKERS \
    --indep_cost $JUMP_COST \
    --runs $RUNS \
    --threads $NTHREADS \
    --distribution $EXTENSION \
    --visible_in $SEE_INCOMING \
    $USE_HYBRID $AVG_GIVEN $REDUCE_VAR $JUMP_STR $NOT_RECURSIVE --revisit_nocost $VS_FLAG \
    > output.${file%.gz} 2> $OUTDIR/${file%.txt.gz}.${SUFFIX}.c${JUMP_COST}.w${JUMP_WEIGHT}.ns${NSTEPS}.b${BUDGET}.s${SEE_INCOMING}.txt


done
