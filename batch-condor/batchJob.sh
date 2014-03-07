#!/bin/bash

source /uscmst1/prod/sw/cms/bashrc prod
# this is needed just to set-up root:
cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_11
#cd /uscmst1/prod/sw/cmssw/slc5_amd64_gcc472/cms/cmssw/CMSSW_6_2_0
eval `scramv1 runtime -sh`
cd -

outDir=$1
count=$2
dataName=$3  #e.g. Photon_Aug05, Run2011B , May10

### Specify addtional arguments here ####
suffix=$4  #e.g. DATA. WZ, ZZ etc
selection=$5
period=$6
gen=$7

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

echo $selection
if [ $selection == "egamma" ]
    then
    echo "Electrons,  double-photon trigger is used"
    ./run.py ${suffix} ${dataName} -e -t pho26 -p ${period}  -b
elif [ $selection == "mugamma" ] 
    then
    echo "Mu+photon trigger will be used"
    if [ $gen == "gen" ] 
        then
        ./run.py ${suffix} ${dataName} -t mugamma -p ${period} --gen -b
    else
        ./run.py ${suffix} ${dataName} -t mugamma -p ${period}  -b
    fi
elif [ $selection == "single-mu" ] 
    then
    echo "Single - IsoMu trigger will be used"
    ./run.py ${suffix} ${dataName} -t single-mu -p ${period} -b
elif [ $selection == "mumu" }
    then
    echo "Mu-mu selections"
    ./run.py ${suffix} ${dataName} -t double-mu -p ${period} -b
else
    echo "No selection provided!"
    exit
fi

echo 'ls in analysis directory'
ls -l

echo "Done. Will copy the files to " $outDir
# for EOS:
cp hhhh_${dataName}.root $outDir/hhhh_${dataName}_${count}.root
# For condor internal output:
cp hhhh_${dataName}.root hhhh_${dataName}_${count}.root
