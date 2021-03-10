#!/bin/sh
pwdd=`pwd`
./g4_run.sh
cd /w/halla-scifs17exp/sbs/thir/git_master/build_NN/
rm -r *
cmake -DCMAKE_INSTALL_PREFIX=/w/halla-scifs17exp/sbs/thir/git_master/libsbsdig-install/ /w/halla-scifs17exp/sbs/thir/git_master/libsbsdig/ && make install
echo $pwdd
cd $pwdd
sbsdig db/db_genrp_conf_dev.dat gmn13.5_elastic_ex.txt 10000
