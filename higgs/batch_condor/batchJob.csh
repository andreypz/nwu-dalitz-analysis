#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_2_8
cd CMSSW_4_2_8/src
cmsenv 

set dir = /uscms_data/d2/andreypz/cmssw/higgs7/CMSSW_4_2_8/src/NWU/Higgs/higgs

cp $dir/higgsAnalyzer.h .
cp $dir/*.C .
cp $dir/../src/TC*.cc .
cp $dir/../src/TC*.h .
mkdir sourceFiles
cp $dir/sourceFiles/*.txt ./sourceFiles/

mkdir ../data
cp -r $dir/../data/*root ../data/.

cp $dir/template.C .
cp $dir/run.csh .
chmod 755 run.csh
./run.csh $1 $2 $3 $4

set newDir = $4
mkdir $newDir
cp hhhh_*.root $newDir
#cp events_* $newDir
#cp counts* $newDir

cp -r $newDir $dir/batch_condor/

cd ../../
rm -rf ./*
