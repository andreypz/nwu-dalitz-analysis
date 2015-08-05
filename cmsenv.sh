#!/bin/bash
# sourse this script in order to set-up root.

PWD=$(pwd)

if [[ $PWD == /tthome/* ]]
then
  echo "NWU"
  # cd /software/tier3/osg/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_11/
  cd /software/tier3/osg/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_2/

else
  #cd /uscms/home/andreypz/work/CMSSW_5_3_20
  cd /cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_20
fi

cmsenv

cd -



