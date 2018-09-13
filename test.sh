#!/bin/bash
# Test script to start the consistency pipeline on primates [9443] with the limit of 10 tree computations

if [ ! -d "data" ]; then
    echo "data directory not found, required for test!"
    exit 1
fi

if [ ! -d "bin" ]; then
    echo "bin directory not found, required for test!"
    exit 1
fi

# clean up for previous run of test
if [ -d "new_definition" ]; then
    echo "cleaning up previous test results"
    rm -r 9443.*
    rm timings.log
    rm -r new_definition
fi

time python og_consistency_pipeline.py 9443 -l 10

# for a complete run of the pipeline parallel execution is reccomended
# e.g. with 10 cores:
# time python og_consistency_pipeline.py 9443 -c 10