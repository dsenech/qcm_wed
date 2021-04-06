#!/bin/sh
make 
mv ~/lib/qcm.so ~/lib/tmp_arm.so
cd ../mac/
arch -x86_64 make
cd ~/lib
mv qcm.so tmp_x86.so
lipo -create -output qcm.so tmp_x86.so tmp_arm.so
rm tmp_x86.so tmp_arm.so
