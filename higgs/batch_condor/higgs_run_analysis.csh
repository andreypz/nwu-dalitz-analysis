#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_2_3
cd CMSSW_4_2_3/src
cmsenv 

set dir =  /uscms_data/d2/andreypz/cmssw/higgs2/CMSSW_4_2_3/src/UserCode/AndreyPozdnyakov/higgs/
set dir =  /uscms_data/d2/andreypz/cmssw/higgs3/CMSSW_4_2_3/src/NWU/Higgs/higgs/
cp $dir/run_analysis.C .
cp $dir/analyzer_higgs.C .
cp $dir/analyzer_higgs.h .
cp $dir/../src/TC* .

# this wierd string is for 'run_analysis.C("2010A","2010B")' kind of thing:
set run =   "run_analysis.C(${1},"\""${2}"\"","\""${3}"\"")"

echo $run

root -b -q ${run}

cp hhhh.root $dir/batch_condor/hhhh_${1}_${2}_${3}.root

set newDir = dir_${1}_${2}_${3}
mkdir $newDir
cp hhhh.root $newDir
cp events_* $newDir
cp counts* $newDir

cp -r $newDir $dir/batch_condor
