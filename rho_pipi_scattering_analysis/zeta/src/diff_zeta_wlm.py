import numpy as np


C_arr = np.genfromtxt("out/C_wlm.out", delimiter = ",")
py_arr = np.genfromtxt("out/py_wlm.out", delimiter = ",")

print(C_arr.shape)
print(py_arr.shape)

m1=0.5
m2_arr=(0.5, 0.7)
q2_arr=(-0.3, 0.3, 0.9, 1.3)
l_arr=(0, 1, 2)
NL_arr=(10, 20)
d1_arr=(0, 1)
d2_arr=(0, 1)
d3_arr=(0, 1)

ind = 0

count1 = 0
count2 = 0

for m2 in m2_arr:
    for q2 in q2_arr:
        for NL in NL_arr:
            for d1 in d1_arr:
                for d2 in d2_arr:
                    for d3 in d3_arr:
                        for l in l_arr:
                            for m in range(-l,l+1):
                                count1+=1
                                diff_re = C_arr[ind][0]-py_arr[ind][0]
                                diff_imag = C_arr[ind][1]-py_arr[ind][1]
                                print(ind, m1, m2, q2, l, m, NL, 6, d1,d2, d3)
                                # print(ind, C_arr[ind][0], py_arr[ind][0], diff_re)
                                # print(ind, C_arr[ind][1], py_arr[ind][1], diff_imag)
                                if abs(diff_re) > 1e-3 or abs(diff_imag) > 1e-3:
                                    count2+=1
                                    # print(m1, m2_arr, q2, l, m, NL, d1,0, d3)
                                    print(ind, C_arr[ind][0], py_arr[ind][0], diff_re)
                                    print(ind, C_arr[ind][1], py_arr[ind][1], diff_imag)
                                    print(ind, "Difference too large!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                                    print()
                                
                                ind += 1
print(count1, count2, count2/count1)


