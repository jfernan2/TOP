#!/bin/bash

echo "Setting up environment..."
source /nfs/fanae/root_releases/root.5.34.13.slc5/bin/thisroot.sh
source /nfs/fanae/PAF_releases/PAF_devel/PAF_setup.sh
cd /nfs/fanae/user/folgueras/TOP

##SEQ PWD=`pwd`
##SEQ mkdir -p jobs
##SEQ 
##SEQ echo "Sent multiple sequential jobs..."
##SEQ SAMPLES="TTGJets TTWJets TTZJets TTWWJets WWGJets WWWJets WWZJets WZZJets ZZZJets DYJets_Madgraph Wbb_Madgraph TWDilep TbarWDilep WgammaToLNuG ZgammaToLLG WWTo2L2Nu_Madgraph WZ ZZ" 
##SEQ 
##SEQ echo "Sending Dibosons, Single Top, Single Boson, TTV and VVV ..."
##SEQ for SAMPLE in $SAMPLES
##SEQ do
##SEQ     echo "Preparing to launch $sample... "
##SEQ     JOBFILE="jobs/$SAMPLE.sh"
##SEQ     cat > $JOBFILE <<EOF
##SEQ #!/bin/bash
##SEQ #PBS -e $PWD/jobs/$SAMPLE.err
##SEQ #PBS -o $PWD/jobs/$SAMPLE.log
##SEQ echo "Setting the environment for ROOT..."
##SEQ source /nfs/fanae/root_releases/root.5.34.13.slc5/bin/thisroot.sh
##SEQ echo "+ ROOTSYS = $ROOTSYS"
##SEQ 
##SEQ echo "Setting the environment for PAF..."
##SEQ source /nfs/fanae/PAF_releases/PAF_devel/PAF_setup.sh 
##SEQ echo "+ PAFPATH = $PAFPATH"
##SEQ 
##SEQ echo "Moving to $PWD..."
##SEQ cd $PWD
##SEQ echo -n "+ "
##SEQ pwd
##SEQ 
##SEQ root -l -b -q "RunTree_ReReco.C(\"${SAMPLE}\", 1, true, false)";
##SEQ 
##SEQ EOF
##SEQ     chmod u+x $JOBFILE
##SEQ     echo "Sending $SAMPLE"
##SEQ     qsub $JOBFILE
##SEQ     sleep 10
##SEQ done

echo "Now start PoD..." 
source /opt/PoD/PoD_env.sh

resetpaf -a

## Data 
root -l -b -q RunTree_ReReco.C\(\"DoubleMu\",30\)
root -l -b -q RunTree_ReReco.C\(\"DoubleElectron\",30\)
root -l -b -q RunTree_ReReco.C\(\"MuEG\",30\)
 
## TTbar (Signal+Bkg)
resetpaf -a
root -l -b -q RunTree_ReReco.C\(\"TTJetsFullLeptMGtauola\",30,true\)
root -l -b -q RunTree_ReReco.C\(\"TTJets_MadSpin\",30,true\)
root -l -b -q RunTree_ReReco.C\(\"TTJets_matchingdown\",30,true\)
root -l -b -q RunTree_ReReco.C\(\"TTJets_matchingup\",30,true\)
root -l -b -q RunTree_ReReco.C\(\"TTJets_scaledown\",30,true\)
root -l -b -q RunTree_ReReco.C\(\"TTJets_scaleup\",30,true\)
root -l -b -q RunTree_ReReco.C\(\"TTJetsSemiLeptMGtauola\",30,true\)

## Single Boson
root -l -b -q RunTree_ReReco.C\(\"ZJets_Madgraph\",30,true\)
resetpaf -a

## Single Top
root -l -b -q RunTree_ReReco.C\(\"DYJets_Madgraph\",15,true\)
root -l -b -q RunTree_ReReco.C\(\"Wbb_Madgraph\",15,true\)
root -l -b -q RunTree_ReReco.C\(\"TbarWDilep\",15,true\)
root -l -b -q RunTree_ReReco.C\(\"TWDilep\",15,true\)

root -l -b -q RunTree_ReReco.C\(\"WgammaToLNuG\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"ZgammaToLLG\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"WWTo2L2Nu_Madgraph\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"WZ\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"ZZ\",15,true\)

root -l -b -q RunTree_ReReco.C\(\"TTGJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"TTWJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"TTZJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"TTWWJets\",15,true\)
root -l -b -q RunTree_ReReco.C\(\"WWGJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"WWWJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"WWZJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"WZZJets\",15,true\) 
root -l -b -q RunTree_ReReco.C\(\"ZZZJets\",15,true\) 


