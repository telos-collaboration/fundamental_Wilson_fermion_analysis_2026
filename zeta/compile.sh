#!/bin/bash

cd zeta/get_wlm/main

make clean
rm get_wlm.so
rm get_wlm.exe
cp Makefile_shared Makefile
make
cp Makefile_orig Makefile
make
cd ../..

cp get_wlm/main/get_wlm.so out/get_wlm.so
cp get_wlm/main/get_wlm.so ../src/get_wlm.so
cp get_wlm/main/get_wlm.exe out/get_wlm.exe