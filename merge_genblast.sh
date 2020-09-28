#!/bin/bash

#myArray=( "$@" )

for arg; do
    cat $arg
    echo ""
done
