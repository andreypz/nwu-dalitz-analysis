#!/bin/csh
source /uscmst1/prod/sw/cms/cshrc prod
# this is needed just to set-up root:
cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_2_0
cmsenv 
cd -

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
cp ../source.tar.gz .
tar -xzf source.tar.gz

cd zgamma

cp ../../input_${dataName}_${count}.txt input.txt

chmod 755 run.py

./run --clean

echo $selection
if ( $selection == "electron" ) then
    echo "Electrons,  double-photon trigger is used"
    ./run.py ${suffix} ${dataName} -e -t pho -p ${period}  -b
else if ( $selection == "mugamma" ) then
    echo "Mu+photon trigger will be used"
    ./run.py ${suffix} ${dataName} -t mugamma -p ${period} -b
else if ( $selection == "single-mu" ) then
    echo "Single - IsoMu trigger will be used"
    ./run.py ${suffix} ${dataName} -t single-mu -p ${period} -b
else 
    echo "Muons"
    ./run.py ${suffix} ${dataName} -t double-mu -p ${period} -b
endif

echo 'ls in analysis directory'
ls -l

echo "Done. Will copy the files to " $outDir
# for EOS:
cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root
# For condor internal output:
cp hhhh_${dataName}.root hhhh_${dataName}_${count}.root
