#!/bin/sh

echo "wphi_run_0.sh"

# Check if configuration file exists.
if [ -e 00-set_env.r ]; then
  echo "00-set_env.r found."
else
  echo "00-set_env.r does not exist."
  exit
fi

# Root of all outputs.
ALL_OUT=`Rscript -e 'source("00-set_env.r");cat(prefix$all.out)'`

# Codes.
CODE_PATH=`Rscript -e 'source("00-set_env.r");cat(prefix$code)'`

# Build the directory and subdirectories.
# rm -rf ${ALL_OUT}
mkdir ${ALL_OUT}
mkdir ${ALL_OUT}/data
mkdir ${ALL_OUT}/subset
mkdir ${ALL_OUT}/log
mkdir ${ALL_OUT}/table

mkdir ${ALL_OUT}/plot
mkdir ${ALL_OUT}/plot/diag
mkdir ${ALL_OUT}/plot/match
mkdir ${ALL_OUT}/plot/single
mkdir ${ALL_OUT}/plot/trace
mkdir ${ALL_OUT}/plot/multi

mkdir ${ALL_OUT}/output
mkdir ${ALL_OUT}/output/tmp

# Make output directories depending on 00-set_env.
Rscript ${CODE_PATH}/01-mkdir_output.r > \
          ${ALL_OUT}/log/01-mkdir_output 2>&1 &
