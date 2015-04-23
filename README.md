git clone https://github.com/folguera/TOP

cd TOP

export PAFPATH=/gpfs/csic_projects/cms/PROOF/paf/
export PATH=$PAFPATH/bin:$PATH
source $PAFPATH/PAF_setup.sh

git version 1.7.1
cd /gpfs/csic_users/piedra/CMSSW_7_3_0/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
git version 1.8.3.1
git clone https://git.cern.ch/reps/IFCA-UO-CMS/Utils

mv ../Utils/PUWeight packages/.

resetpaf -a

# Test
root -l -b -q RunTree_ReReco.C\(\"WZ\",4,true\)


# Make the top trees
./runAll_TopMC.sh

# Compute systematics and top cross section
correr RunTopPlots.C.
