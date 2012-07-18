#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_4_4
cd CMSSW_4_4_4/src
cmsenv 

#### Leave this blank #######

#############################

set srcDir    = $1
set outDir    = $2
set count     = $3
set dataName  = $4  #e.g. Photon_Aug05, Run2011B , May10

### Specify addtional arguments here ####
set suffix    = $5  #e.g. DATA. WZ, ZZ etc
set trigger   = $6
set selection = $7
set period    = $8


set dir = ${srcDir}/..

echo "Copy files"

cp $dir/higgsAnalyzer.h .
cp $dir/*.C .
cp $dir/config.py .
cp $dir/../src/*.cc .
cp $dir/../src/*.h .

mkdir ../plugins
cp $dir/../plugins/*.h  ../plugins
cp $dir/../plugins/*.cc ../plugins
cp $dir/../plugins/libShapeLine.so ../plugins


mkdir ../data
cp -r $dir/../data/*.root ../data/

cp $dir/run.py .
chmod 755 run.py
./run.py ${suffix} ${selection} ${dataName} ${trigger} ${period} b

cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root
cp  events_printout_* $outDir/printouts/
