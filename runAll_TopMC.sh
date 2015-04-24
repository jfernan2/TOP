#!/bin/bash

echo "Setting up Oviedo environment..."
source /nfs/fanae/PAF_releases/PAF_devel/PAF_setup.sh
source /opt/root/bin/thisroot.sh
source /opt/PoD/PoD_env.sh

#echo "Setting up IFCA environment..."
#export PAFPATH=/gpfs/csic_projects/cms/PROOF/paf/
#export PATH=$PAFPATH/bin:$PATH
#source $PAFPATH/PAF_setup.sh


## Data 
resetpaf -a

root -l -b -q 'RunTree_ReReco.C("DoubleMu",       30)'
root -l -b -q 'RunTree_ReReco.C("DoubleElectron", 30)'
root -l -b -q 'RunTree_ReReco.C("MuEG",           30)'
 
## TTbar (signal+background)
resetpaf -a

root -l -b -q 'RunTree_ReReco.C("TTJetsFullLeptMGtauola", 30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJets_MadSpin",         30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJetsSemiLeptMGtauola", 30, true)'

## TTbar (systematic samples) 
resetpaf -a

root -l -b -q 'RunTree_ReReco.C("TTJets_MadSpinPDF",           30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJets_matchingdown",         30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJets_matchingup",           30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJets_scaledown",            30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJets_scaleup",              30, true)'
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg",                30, true)'
root -l -b -q 'RunTree_ReReco.C("TTbar_Powheg_Herwig",         30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJetsFullLeptMGTuneP11",     30, true)'
root -l -b -q 'RunTree_ReReco.C("TTJetsFullLeptMGTuneP11noCR", 30, true)'

root -l -b -q 'RunTree_ReReco.C("TTbar_PowhegV2_hdampWeights", 30)'
root -l -b -q 'RunTree_ReReco.C("TTJetsHadrMGTuneP11mpiHi",    30)'
root -l -b -q 'RunTree_ReReco.C("TTJetsFullLeptMGTuneP11TeV",  30)'

##root -l -b -q 'RunTree_ReReco.C("TTbar_MCatNLO",                   30, true)
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_mass169",          30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_mass175",          30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_scaledown",        30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_scaleup",          30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTbar_MCatNLO_noCorr",            30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_noCorr_mass169",   30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_noCorr_mass175",   30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_noCorr_scaledown", 30, true)'
##root -l -b -q 'RunTree_ReReco.C("TTLept_mcatnlo_noCorr_scaleup",   30, true)'

## Backgrounds
resetpaf -a

root -l -b -q 'RunTree_ReReco.C("ZJets_Madgraph",  30, true)'
root -l -b -q 'RunTree_ReReco.C("DYJets_Madgraph", 30, true)'
root -l -b -q 'RunTree_ReReco.C("Wbb_Madgraph",    30, true)'
root -l -b -q 'RunTree_ReReco.C("TbarWDilep",      30, true)'
root -l -b -q 'RunTree_ReReco.C("TWDilep",         30, true)'

root -l -b -q 'RunTree_ReReco.C("WgammaToLNuG",       30, true)'
root -l -b -q 'RunTree_ReReco.C("ZgammaToLLG",        30, true)'
root -l -b -q 'RunTree_ReReco.C("WWTo2L2Nu_Madgraph", 30, true)'
root -l -b -q 'RunTree_ReReco.C("WW",                 30, true)'
root -l -b -q 'RunTree_ReReco.C("WZ",                 30, true)'
root -l -b -q 'RunTree_ReReco.C("ZZ",                 30, true)'

root -l -b -q 'RunTree_ReReco.C("TTGJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("TTWJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("TTZJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("TTWWJets", 30, true)'
root -l -b -q 'RunTree_ReReco.C("WWGJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("WWWJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("WWZJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("WZZJets",  30, true)' 
root -l -b -q 'RunTree_ReReco.C("ZZZJets",  30, true)'

## Stop
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 150.0,  1.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 162.5,  1.0)'
 root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 175.0,  1.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 187.5,  1.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 200.0,  1.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 162.5, 12.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 175.0, 12.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 187.5, 12.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 200.0, 12.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 212.5, 12.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 175.0, 25.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 187.5, 25.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 200.0, 25.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 212.5, 25.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 225.0, 25.0)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 187.5, 37.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 200.0, 37.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 212.5, 37.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 225.0, 37.5)'
#root -l -b -q 'RunTree_ReReco.C("T2tt_150to250LSP1to100_LeptonFilter", 10, true, -1, 237.5, 37.5)'
