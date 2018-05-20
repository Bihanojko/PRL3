#!/usr/bin/env bash
# Project: Prirazeni poradi preorder vrcholum
# Author:  Nikola Valesova, xvales02
# Date:    25. 4. 2018


# check number of arguments
if [ $# -ne 1 ]; then
    echo "ERROR: Invalid number of arguments!"
    exit 1
fi

# compile source files
mpic++ --prefix /usr/local/share/OpenMPI -o pro pro.cpp

# execute
if [ ${#1} -eq 1 ]; then
    mpirun --prefix /usr/local/share/OpenMPI -np 1 pro $1
else
    mpirun --prefix /usr/local/share/OpenMPI -np $(((${#1} - 1) * 2)) pro $1 2>/dev/null
fi

# clean
rm -f pro
