#!/bin/bash

bash compile.sh

rm out/C_wlm.out
rm out/py_wlm.out
rm out/diff_wlm.out

touch out/C_wlm.out
touch out/py_wlm.out
touch out/diff_wlm.out

m1=0.5
m2_arr=(0.5 0.7)
q2_arr=(-0.3 0.3 0.9 1.3)
l_arr=(0 1 2)
NL_arr=(10 20)
d1_arr=(0 1)
d2_arr=(0 1)
d3_arr=(0 1)


for m2 in "${m2_arr[@]}"
do
for q2 in "${q2_arr[@]}"
do
for NL in "${NL_arr[@]}"
do
for d1 in "${d1_arr[@]}"
do
for d2 in "${d2_arr[@]}"
do
for d3 in "${d3_arr[@]}"
do
for l in "${l_arr[@]}"
do
for m in $(seq -$l $l);
do
echo m1${m1}m2${m2}q2${q2}l${l}m${m}NL${NL}sz6d${d1}${d2}${d3} 
./get_wlm/main/get_wlm.exe $m1 $m2 $q2 $l $m $NL 6 $d1 $d2 $d3 >> out/C_wlm.out
python3 src/pylink_wlm.py $m1 $m2 $q2 $l $m $NL 6 $d1 $d2 $d3 >> out/py_wlm.out
done
done
done
done
done
done
done
done

diff out/C_wlm.out out/py_wlm.out >> out/diff_wlm.out