
### to download current HiPR source code
``
git clone https://github.com/wanglab-upenn/HiPR.git
``

### to recompile HiPR_MCMC
Linux:
``
g++ -I/usr/local/include/ HiPR_MCMC.cpp
``
Mac:
``
g++ -stdlib=libstdc++ HiPR_mapseq.cpp -o HiPR_mapseq.MAC
``

### to install required PERL modules
``
perl -MCPAN -e "install 'File::HomeDir'"
perl -MCPAN -e "install 'Forks::Super'"
``
If *Forks::Super* module fails to install due to failed tests, try skipping testing altogether:
1. using shell
``
perl -MCPAN -e shell
notest install Forks::Super
``
2. or using command line
``
perl -MCPAN -e "notest('install','Forks::Super')"
``
