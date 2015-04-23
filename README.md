Everything starts here
====

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /gpfs/csic_projects/cms/piedra/work    


It is time to get the material
====

Go to the master repository (https://github.com/folguera/TOP) and click **Fork** in the top-right corner of the page. Now get the code in your working area.

    git clone https://github.com/piedraj/TOP TOP

    pushd /gpfs/csic_users/piedra/CMSSW_7_3_0/src
    cmsenv
    popd
    git clone https://git.cern.ch/reps/IFCA-UO-CMS/Utils

    cd TOP
    mv ../Utils/PUWeight packages/.


Setup PAF
====

    export PAFPATH=/gpfs/csic_projects/cms/PROOF/paf/
    export PATH=$PAFPATH/bin:$PATH
    source $PAFPATH/PAF_setup.sh


Test the code
====

    resetpaf -a
    root -l -b -q RunTree_ReReco.C\(\"WZ\",4,true\)


Make the top trees
====

    resetpaf -a
    ./runAll_TopMC.sh


Compute the top cross section
====

It includes the systematic uncertainties.

    root -l RunTopPlots.C

