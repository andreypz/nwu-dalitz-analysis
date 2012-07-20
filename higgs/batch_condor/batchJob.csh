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


set dir = ${srcDir}/code_dir

echo "Copy all files needed from a working directory"
cd ../
cp -r $dir/* .


cd higgs
mv ../src/input.txt .

#echo 'ls in higgs'
#ls

chmod 755 run.py
./run.py ${suffix} ${selection} ${dataName} ${trigger} ${period} b

cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root
cp  events_printout_* $outDir/printouts/
