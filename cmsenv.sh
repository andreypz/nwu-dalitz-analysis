#!/bin/bash
# sourse this script in order to set-up root.

PWD=$(pwd)

if [[ $PWD == /tthome/* ]]
then
  echo "NWU"
  #cd /software/tier3/osg/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_11/
  cd /software/tier3/osg/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_2/

else
  cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_11
fi

cmsenv

cd -



