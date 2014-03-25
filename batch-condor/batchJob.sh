#!/bin/bash

outDir=$1
count=$2
dataName=$3  #e.g. Photon_Aug05, Run2011B , May10

### Specify addtional arguments here ####
suffix=$4  #e.g. DATA. WZ, ZZ etc
selection=$5
trigger=$6
period=$7
gen=$8
nwu=$9

echo ${7} and ${8} and ${8}

if [ $nwu == "nwu" ]
then
  echo "NWU"
  export OSG_APP=/software/tier3/osg
  export SCRAM_ARCH=slc6_amd64_gcc472
  source /software/tier3/osg/cmsset_default.sh
  cd /software/tier3/osg/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_11
else
  source /uscmst1/prod/sw/cms/bashrc prod
  # this is needed just to set-up root:
  cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_2_0
fi
eval `scramv1 runtime -sh`
cd -


echo "Copy the tarball with the source code from Condor scratch directory"

mkdir nwu-my-analysis
cd nwu-my-analysis
cp ../source.tar.gz .
tar -xzf source.tar.gz

cd zgamma

cp ../../input_${dataName}_${count}.txt input.txt
ls -l

chmod 755 run.py

./run.py --clean

echo ${selection} ${trigger}

if [ ${gen} == "gen" ]
then
    ./run.py ${suffix} ${dataName} -s ${selection} -t ${trigger} -p ${period}  -b
else
    ./run.py ${suffix} ${dataName} -s ${selection} -t ${trigger} -p ${period} --gen -b
fi

echo 'ls in analysis directory'
ls -l

echo "Done. Will copy the files to " $outDir
# for EOS:
cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root
# For condor internal output:
cp hhhh_${dataName}.root hhhh_${dataName}_${count}.root
