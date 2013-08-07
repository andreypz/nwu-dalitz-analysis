#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
scram pro CMSSW CMSSW_5_3_8
cd CMSSW_5_3_8/src
cmsenv 

#### Leave this blank #######

#############################

set outDir    = $1
set count     = $2
set dataName  = $3  #e.g. Photon_Aug05, Run2011B , May10

### Specify addtional arguments here ####
set suffix    = $4  #e.g. DATA. WZ, ZZ etc
set selection = $5
set period    = $6


echo "Copy the tarball with the source code from Condor scratch directory"


mkdir nwu-my-analysis
cd nwu-my-analysis
cp ../../../source.tar.gz .
tar -xzf source.tar.gz

cd zgamma

#ls -l
cp ../../../../input_${dataName}_${count}.txt input.txt
#ls -l ../../../../

chmod 755 run.py

echo $selection
if ( $selection == "electron" ) then
    echo "Electrons"
    ./run.py ${suffix} ${dataName} -e -p ${period} -b
else if ( $selection == "mugamma" ) then
    echo "Mu+photon trigger will be used"
    ./run.py ${suffix} ${dataName} --mugamma -p ${period} -b
else
    echo "Muons"
    ./run.py ${suffix} ${dataName} -p ${period} -b
endif

echo 'ls in analysis directory'
ls -l

echo "Done. Will copy the files to " $outDir
cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root
#cp  events_printout_* $outDir/$selection/printouts/
