#!/bin/bash

cd zeta/get_wlm/main

make clean
make -f Makefile_shared
make -f Makefile_orig
cd ../..

cp get_wlm/main/get_wlm.so out/get_wlm.so
cp get_wlm/main/get_wlm.exe out/get_wlm.exe