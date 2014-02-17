#!/bin/bash
# sourse this script in order to set-up root.

if [ "${1}" == "nwu" ]
then
  echo "NWU"
  cd /software/tier3/osg/slc6_amd64_gcc462/cms/cmssw/CMSSW_6_0_1/
else
  cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_2_0
fi

cmsenv 

cd -



