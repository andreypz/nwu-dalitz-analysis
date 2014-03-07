#!/bin/bash

PWD=$(pwd)

if [[ $PWD == /tthome/* ]]
then
    echo "NWU"
    cd /tthome/andrey/cmssw/CMSSW_6_1_1/src
else
    cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc462/cms/cmssw/CMSSW_6_1_1
fi

cmsenv 

cd -



