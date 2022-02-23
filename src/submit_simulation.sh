#!/bin/bash
#-------------input parameters--------------------------------------------------

cd /nfs/dust/ilc/user/fstehr/leap_sims/submit_files
condor_submit simulation${1}.sub
cd -
echo "Simulation ${1} submited"
