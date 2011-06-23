#!/bin/csh

set lepton = 2 # 1- muons, 2 -electrons

set dataset = (\
 2011A_May10_DoubleMu \
 2011A_HZZSkim_DoubleMu \
 2011A_May10_DoubleElectron\
 2011A_HZZSkim_DoubleElectron\
 MC_ZZtoAny\
 MC_Wjets\
 MC_ZllG\
 MC_ZeeGamma\
 MC_tt2l2nu2b\
 MC_tSchannel\
 MC_tTchannel\
 MC_tWchannel\
 MC_WW\
 MC_WZ\
 MC_DYmumu\
 MC_DYee\
 MC_DYtautau\
 MC_Zbb0\
 MC_Zbb1\
 MC_Zbb2\
 MC_Zbb3\
 MC_Zcc0\
 MC_Zcc1\
 MC_Zcc2\
 MC_Zcc3\
 MC_ggH200\
 MC_ggH400\
 MC_ggH600\
)
# MC_ggH200_HWW2l2nu\
# MC_ggH400_HWW2l2nu\
# MC_ggH600_HWW2t2nu\
# MC_ggH200_HWW2t2nu\
# MC_ggH400_HWW2t2nu\
# MC_ggH600_HWW2t2nu\
# )

foreach name ($dataset)
echo $name

cp condor_job temp_job
sed -i "s:<lepton>:${lepton}:g" temp_job
sed -i "s:<dataset>:${name}:g" temp_job 
condor_submit temp_job
rm temp_job
end

