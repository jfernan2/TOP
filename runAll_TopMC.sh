#!/bin/bash

if [[ `hostname -d` == *"uniovi"* ]]; then
    echo "Setting up Oviedo environment..."
    source /nfs/fanae/PAF_releases/PAF_devel/PAF_setup.sh
    source /opt/root/bin/thisroot.sh
    source /opt/PoD/PoD_env.sh
else
    echo "Setting up IFCA environment..."
    export PAFPATH=/gpfs/csic_projects/cms/PROOF/paf/
    export PATH=$PAFPATH/bin:$PATH
    source $PAFPATH/PAF_setup.sh
fi

resetpaf -a

root -l -b -q 'RunTree_ReReco.C("TestSync", 4, true)'

#root -l -b -q 'RunTree_ReReco.C("TTJets_MadSpin", 10, true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Pythia"  ,  4, true)'
#root -l -b -q 'RunTree_ReReco.C("TW"		 ,  4, true)'
#root -l -b -q 'RunTree_ReReco.C("TbarW" 	 ,  4, true)'
#root -l -b -q 'RunTree_ReReco.C("ZJets_Madgraph",  4, true)'
#root -l -b -q 'RunTree_ReReco.C("WJets_Madgraph",  4, true)'
#root -l -b -q 'RunTree_ReReco.C("WZTo3LNu"	 ,  4, true)'
#root -l -b -q 'RunTree_ReReco.C("ZZ4L"  	 ,  4, true)'
#root -l -b -q 'RunTree_ReReco.C("TTWJets"	 ,  4, true)'
#root -l -b -q 'RunTree_ReReco.C("TTZJets"	 ,  4, true)'

#root -l -b -q 'RunTree_ReReco.C("TTbar_Pythia_scaledown", 4,true)'
#root -l -b -q 'RunTree_ReReco.C("TTbar_Pythia_scaleup"  , 4,true)'
