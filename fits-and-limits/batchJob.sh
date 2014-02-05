#!/bin/bash

source /uscmst1/prod/sw/cms/bashrc prod
# this is needed just to set-up root:
cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_2_0
eval `scramv1 runtime -sh`
cd -

outDir=$1
jobNumber=$2
nTrials=$3  # trials per job
mass=$4

mkdir nwu-my-analysis
cd nwu-my-analysis
cp ../source.tar.gz .
tar -xzf source.tar.gz
echo "lsing in nwu-my-analysis"
ls
cd fits-and-limits
chmod 755 biasStudy_toyMaker.py

echo "lsing in fits-and-limits"
ls

echo "$outDir $jobNumber $nTrials"
./biasStudy_toyMaker.py -v ${outDir} -j $jobNumber -t $nTrials -m $mass


echo "Done. Will copy the files to " $outDir



cp $outDir/* ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

ls

echo "Now completely done"


