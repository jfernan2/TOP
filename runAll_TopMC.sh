#!/bin/bash

echo "Setting up environment..."
#source /nfs/fanae/root_releases/root.5.34.13.slc5/bin/thisroot.sh
source /nfs/fanae/PAF_releases/PAF_devel/PAF_setup.sh
source /opt/root/bin/thisroot.sh
#cd /nfs/fanae/user/sscruz/TOP

#source /nfs/fanae/user/sscruz/ChargeMiss/MissCode/setup.sh

echo "Now start PoD..." 
source /opt/PoD/PoD_env.sh


resetpaf -a

root -l -b -q 'RunTree_ReReco.C("TestSync", 1, true)'

#root -l -b -q 'RunTree_ReReco.C("TTJets_MadSpin",10, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Pythia"  , 1, true)'
#root -l -b -q 'RunTree_ReReco.C("TW"		 , 1, true)'
#root -l -b -q 'RunTree_ReReco.C("TbarW" 	 , 1, true)'
#root -l -b -q 'RunTree_ReReco.C("ZJets_Madgraph", 1, true)'
#root -l -b -q 'RunTree_ReReco.C("WJets_Madgraph", 1, true)'
#root -l -b -q 'RunTree_ReReco.C("WZTo3LNu"	 , 1, true)'
#root -l -b -q 'RunTree_ReReco.C("ZZ4L"  	 , 1, true)'
#root -l -b -q 'RunTree_ReReco.C("TTWJets"	 , 1, true)'
#root -l -b -q 'RunTree_ReReco.C("TTZJets"	 , 1, true)'

#root -l -b -q 'RunTree_ReReco.C("TTbar_Pythia_scaledown", 1,true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Pythia_scaleup"  , 1,true)'
