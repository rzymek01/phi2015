### Message passing simulation on Intel Xeon Phi
See also https://github.com/rzymek01/message-passing-simulation-in-cuda

# Dev
`source /opt/intel/composer_xe_2013_sp1.3.174/bin/compilervars.sh intel64`

`icc main.cpp -fopenmp -vec-report3 -o msg`

`export MIC_ENV_PREFIX=MIC`

`export MIC_KMP_AFFINITY=balanced`

