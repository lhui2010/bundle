#!/bin/bash

function cdf() {
    THEFILE=$1
    # echo "cd into directory of ${THEFILE}"
    # For Mac, replace find with mdfind to get it a lot faster. And it does not need args ". -name" part.
    THEDIR=$(dirname "$THEFILE")
    echo "cd into directory of ${THEDIR}"
    cd ${THEDIR}
}

