Everything starts here
====

    ssh -Y gridui.ifca.es -o ServerAliveInterval=240
    cd /gpfs/csic_projects/cms/piedra/work    


It is time to get the material
====

Go to the master repository (https://github.com/folguera/TOP) and click **Fork** in the top-right corner of the page. Now get the code in your working area.

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
    root -l -b -q RunTree_ReReco.C\(\"WZ\",4,true\)


Make the top trees
====

    resetpaf -a
    ./runAll_TopMC.sh


Compute the top cross section
====

    root -l RunTopPlots.C


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

