Everything starts here
====

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    cd /gpfs/csic_users/piedra/CMSSW_7_3_0/src
    cmsenv


It is time to get the material
====

Go to the master repository (https://github.com/folguera/TOP) and click **Fork** in the top-right corner of the page. Now get the code in your working area.

    cd /gpfs/csic_projects/cms/piedra/work    

    git clone https://github.com/piedraj/TOP
    git clone https://piedra@git.cern.ch/reps/IFCA-UO-CMS/Utils

    mv Utils/PUWeight TOP/packages/.


Setup PAF
====

    export PAFPATH=/gpfs/csic_projects/cms/PROOF/paf/
    export PATH=$PAFPATH/bin:$PATH
    source $PAFPATH/PAF_setup.sh


Test the code
====

    cd TOP
    resetpaf -a
    root -l -b -q 'RunTree_ReReco.C("WZ",4,true)'


Make the top trees
====

    ./runAll_TopMC.sh


Compute the top cross section
====

    cd TopCode/
    root -l 'RunTopPlots.C("/gpfs/csic_projects/cms/piedra/work/TOP/TopTrees",0)'


It is commit time
====

First get the latest changes in the repository, if any.

    git pull https://github.com/folguera/TOP

And then commit your changes.

    git status
    git add <filepattern>
    git commit -m 'Modified'
    git push origin master

If the changes have been made in a fork of the master, go to https://github.com/piedraj/TOP and click **Pull Request**.

