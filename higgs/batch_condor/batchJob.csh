#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_4_2_8
cd CMSSW_4_2_8/src
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
cp $dir/../src/*.cc .
cp $dir/../src/*.h .

mkdir ../plugins
cp $dir/../plugins/*.h  ../plugins
cp $dir/../plugins/*.cc ../plugins


mkdir ../data
cp -r $dir/../data/*.root ../data/

cp $dir/run.csh .
chmod 755 run.csh
./run.csh ${suffix} ${trigger} ${dataName} ${selection} ${period} b


#set printDir =  printout_${dataName}
#mkdir $printDir
#cp events_* $printDir
#cp -r $printDir $outDir
cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root

#cp counts* $newDir
