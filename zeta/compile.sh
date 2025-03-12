#!/bin/bash
cd $(dirname "$0")/get_wlm/main

make clean
make -f Makefile_shared
make -f Makefile_orig

cp get_wlm.so  ../../out/get_wlm.so
cp get_wlm.exe ../../out/get_wlm.exe
cd -