#!/bin/bash

PWD=$(pwd)

if [[ $PWD == /tthome/* ]]
then
    echo "NWU"
    cd /tthome/andrey/cmssw/CMSSW_6_1_1/src
else
    cd /uscms/home/andreypz/work/CMSSW_7_1_5
fi

cmsenv 

cd -



