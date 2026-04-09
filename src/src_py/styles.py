import matplotlib
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [10, 6] 
fontsize = 18
font = {'size' : fontsize}
matplotlib.rc('font', **font)
plt.rcParams.update({"mathtext.fontset": "cm" }) # Computer Modern

c_14 = "olivedrab"
c_10_res = "royalblue"
c_10_non_res = "orangered"