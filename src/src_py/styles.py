import matplotlib
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [10, 6] 
fontsize = 16
font = {'size'   : fontsize}
matplotlib.rc('font', **font)
plt.rcParams.update({
    # "font.family": "serif",
    "mathtext.fontset": "cm",   # Computer Modern
})