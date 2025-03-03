import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
import sys
import os.path

def read_from_hdf(filename):
    res, res_tmp = [{},{}]
    with h5py.File("output/hdf5/"+filename+".hdf5","r") as hfile:
        for key in hfile.keys():
            if key[:4] == "orig":
                res[key[5:]] = hfile[key][()]
            if key[:4] == "samp":
                res_tmp[key[7:]] = hfile[key][()]
    return res, res_tmp

def error_of_array_gauss(array):
    tempo = []
    for i in range(len(array[0])):
        tmp = array[:,i]
        tmp.sort()
        tempo.append(tmp)
    num_perc = math.erf(1/np.sqrt(2))
    length = len(array)
    result = np.transpose(tempo)
    return result[length//2], abs(result[length//2]-result[math.floor(length*(1-num_perc)/2)]), abs(result[length//2]-result[math.ceil(length*(1+num_perc)/2)]),result[math.floor(length*(1-num_perc)/2)],result[math.ceil(length*(1+num_perc)/2)]

def error_of_array_lin(arr):
    length = len(arr)
    obs = len(arr[0])
    return arr[length//2], [abs(arr[length//2][i]-arr[0][i]) for i in range(obs)], [abs(arr[length//2][i]-arr[length-1][i]) for i in range(obs)], arr[0], arr[length-1]

def error_of_array(resampling = "gauss"):
    resampling = bytes.decode(resampling)
    if resampling == "gauss":
        return error_of_array_gauss
    elif resampling == "lin":
        return error_of_array_lin

def en_L(N_L,p2,mpi):                                           # wrong formula for q^2 \elem Z
    return 1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2)
def en_rho_L(N_L,p2,mrho,mpi):
    return np.sqrt(mrho**2+p2*(2*np.pi/N_L)**2)/mpi

def marker(pind):
    if pind == 1:
        return "o"
    else:
        return "x"
def colorf(en_lv_ind):
    if en_lv_ind == 1:
        return "blue"
    else:
        return "red"

# def plot_E(file, show=False, save = True, which = "aE", pref = ""):                         # not up to date
#     plt.rcParams['figure.figsize'] = [10, 6]
#     fontsize = 14
#     font = {'size'   : fontsize}
#     matplotlib.rc('font', **font)
#     fig, ax = plt.subplots()
#     res,  res_sample =read_from_hdf(file)
#     num_gaussian = len(res_sample["En"])

#     mpi = res["mpi"]
#     mrho = res["mrho"]
#     N_Ls = res["N_L"]
#     lvls = res["en_lv"]
#     dvecs = res["dvec"]
#     ds = res["d"]
#     d2s = res["d2"]

#     ax.set_xlabel("N_L")
#     ax.grid()

#     xarr = np.linspace(min(N_Ls)-1, max(N_Ls)+1)
#     if which == "aE":
#         ax.set_ylabel("aE")
#         ax.set_ylim([0.5,1.1])
#         ax.set_ylabel("aE")
#         pipi1arr = [mpi*en_L(x,1,mpi) for x in xarr]
#         pipi2arr = [mpi*en_L(x,2,mpi) for x in xarr]
#         rho1arr = [mpi*en_rho_L(x,1,mrho,mpi) for x in xarr]
#         rho2arr = [mpi*en_rho_L(x,2,mrho,mpi) for x in xarr]
#         ax.axhline(mpi*2, color = "green", label = "pipi")
#         ax.axhline(mpi*4, color = "green", lw = 2)
#         ax.axhline(mrho, color = "orange", label = "rho")
#         ax.plot(xarr, pipi1arr, ls = "dotted", color = "green")
#         ax.plot(xarr, pipi2arr, ls = "dashdot", color = "green")
#         ax.plot(xarr, rho1arr, ls = "dotted", color = "orange")
#         ax.plot(xarr, rho2arr, ls = "dashdot", color = "orange")
        
#         y_plot = error_of_array(res["resampling"])(res_sample["En"])
#     elif which == "sqrt_s":
#         ax.set_ylim([1.9,2.4])
#         ax.set_ylabel("$\sqrt{s}$/$m_\pi$")
#         pipi11arr = [sqrt_s_pi_L(x,1,1) for x in xarr]
#         # pipi12arr = [sqrt_s_pi_L(x,2,1) for x in xarr]
#         # pipi21arr = [sqrt_s_pi_L(x,1,2) for x in xarr]
#         pipi22arr = [sqrt_s_pi_L(x,2,2) for x in xarr]
#         ax.axhline(2, color = "green", label = "pipi")
#         ax.axhline(mrho/mpi, color = "orange", label = "rho")
#         ax.plot(xarr, pipi11arr, ls = "dotted", color = "green")
#         # ax.plot(xarr, pipi12arr, ls = "dashdot", color = "green")
#         # ax.plot(xarr, pipi21arr, ls = "dotted", color = "green")
#         ax.plot(xarr, pipi22arr, ls = "dashdot", color = "green")
        
#         y_plot = error_of_array(res["resampling"])(res_sample["E_cm_prime"])
#     for i in range(len(y_plot[0])):
#         ax.errorbar(N_Ls[i], y_plot[0][i], yerr = [[y_plot[1][i],],[y_plot[2][i],]], capsize=3, ls="", marker = ms_P(dvecs[i]), color = color_NL(N_Ls[i]), markersize = ms_size(lvls[i]), label = "Lv: %i, |P|^2=%i, NL=%i"%(lvls[i],d2s[i],N_Ls[i]))
   
#     ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     if save:
#         fig.savefig("output/plots/%s_L.pdf"%which, bbox_inches='tight')
#     if show:
#         plt.show()
#     fig.clf()


def color_NL(NL):
    if NL == 14:
        return "red"
    elif NL == 16:
        return "green"
    elif NL == 24:
        return "blue"

def ls_P(dvec):
    if list(dvec) == [0,0,1]:
        return "solid"
    elif list(dvec) == [1,1,0]:
        return "dashed"

def ms_P(dvec):
    if list(dvec) == [0,0,1]:
        return "*"
    elif list(dvec) == [1,1,0]:
        return "o"

def ms_size(level):
    if level == 1:
        return 8
    elif level == 2:
        return 15

def delete_steps(arr, sign = 1, delete=True):
    if delete:
        for i in range(len(arr)-1):
            if arr[i+1] < sign*arr[i]: 
                arr[i] = np.nan
        return arr
    else:
        return arr

def p3_cot_PS(file, show=False, save = True, pref = "", x_ax = "sqrt_s", y_ax = "p3cotPS_Ecm", ld="", prime = "", delete = True):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    res,  res_sample = read_from_hdf(file)
    num_gaussian = len(res_sample["E_cm_prime"])

    lvls = res["en_lv"] 
    N_Ls = res["N_L"]    
    plt.grid()

    mpi_pr_s = ""
    if not prime == "":
        mpi_pr_s = r"/m_{\pi}"

    if x_ax == "sqrt_s":
        plt.xlabel("$E_{CM}"+mpi_pr_s+"$")
        x_plot = res["E_cm"+ld+prime]
        x_plot_sam = np.transpose(res_sample["E_cm"+ld+prime])
        ax.set_xlim([2,2.5])
    elif x_ax == "s":
        plt.xlabel("$E_{CM}^2"+mpi_pr_s)                                    # needs fix
        x_plot = res["s"+ld+prime]
        x_plot_sam = np.transpose(res_sample["s"+ld+prime])
        ax.set_xlim([0.1,0.4])
        # ax.set_xlim([4,7])
    elif x_ax == "pstar2":
        plt.xlabel(r"$p^{\star^2}"+mpi_pr_s+"^2$")
        x_plot = res["p2star"+ld+prime]
        x_plot_sam = np.transpose(res_sample["p2star"+ld+prime])
        ax.set_xlim([0,0.5])
    elif x_ax == "aE":
        plt.xlabel("aE")
        x_plot = res["En"]
        x_plot_sam = np.transpose(res_sample["En"])
        ax.set_xlim([0.8, 1.15])
    elif x_ax == "En":
        plt.xlabel("$E"+mpi_pr_s+"$")
        x_plot = res["En"+prime]
        x_plot_sam = np.transpose(res_sample["En"+prime])
        ax.set_xlim([2.1, 3])
        
    if y_ax == "p3cotPS":
        y_plot = np.real(res["p3cotPS"+ld+prime])
        y_plot_sam = np.transpose(np.real(res_sample["p3cotPS"+ld+prime]))
        plt.ylabel(r"$p^3\, \cot(\delta)"+mpi_pr_s+"^3$")    
        ax.set_ylim([-5,5])
    elif y_ax == "p3cotPS_Ecm":
        y_plot = np.real(res["p3cotPS_Ecm"+ld+prime])
        y_plot_sam = np.transpose(np.real(res_sample["p3cotPS_Ecm"+ld+prime]))
        plt.ylabel(r"$p^3\, \cot(\delta)/E_CM"+mpi_pr_s+"^3$")    
        # ax.set_ylim([-0.15,0.05])
        ax.set_ylim([-1,0])
    elif y_ax == "sigma":
        y_plot = np.real(res["sigma"+ld+prime])
        y_plot_sam = np.transpose(np.real(res_sample["sigma"+ld+prime]))
        plt.ylabel(r"$\sigma_1*"+mpi_pr_s+"^2$")    
        # ax.set_ylim([-0.15,0.05])
    elif y_ax == "cot_PS":
        y_plot = np.real(res["cot_PS"+ld])
        y_plot_sam = np.transpose(np.real(res_sample["cot_PS"+ld]))
        plt.ylabel(r"cot($\delta$)")      
        ax.set_ylim([-100,100])
    elif y_ax == "PS":
        y_plot = np.real(res["PS+ld"])
        y_plot_sam = np.transpose(np.real(res_sample["PS"+ld]))
        plt.ylabel("$q^2$")      
        ax.set_ylim([0,180])
    elif y_ax == "q2":
        y_plot = np.real(res["q2"+ld])
        y_plot_sam = np.transpose(np.real(res_sample["q2"+ld]))
        plt.ylabel("$q^2$")    

    N_Ls = res["N_L"]
    dvecs = res["dvec"]
    d2s = res["d2"]

    length = len(x_plot_sam[0])


    num_perc = math.erf(1/np.sqrt(2))
    for i in range(len(N_Ls)):
        if lvls[i] != 1:
            ax.scatter(x_plot[i],y_plot[i], color = color_NL(N_Ls[i]), label = "Lv: %i, |P|^2=%i, NL=%i"%(lvls[i],d2s[i],N_Ls[i]), marker = ms_P(dvecs[i]), s=10*ms_size(lvls[i]))
            if bytes.decode(res["resampling"]) == "gauss":
                sorted_indices = np.argsort(x_plot_sam[i])  
                ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_plot_sam[i][sorted_indices],delete=delete)[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))
                # ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],y_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))
            elif bytes.decode(res["resampling"]) == "lin":
                ax.plot(x_plot_sam[i],delete_steps(y_plot_sam[i],delete=delete), color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))
                # ax.plot(x_plot_sam[i],y_plot_sam[i], color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))

    ax.legend(loc="best")
    if save:
        fig.savefig("output/plots/%s_%s"%(y_ax,x_ax)+ld+prime+"_"+pref+".pdf", bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def print_plymouth_table(file, ld = "_ld"):
    res,  res_sample = read_from_hdf(file)
    mpi = res["mpi"]
    E_n = res["En"]
    E_cm = res["E_cm"+ld]
    q2 = res["q2"+ld]
    p3cotPS_Ecm = res["p3cotPS_Ecm"+ld]
    PS = res["PS"+ld]
    for i in range(len(E_n)):
        print("%f\t%f\t%f\t%f\t%f\t%f"%(E_n[i],E_cm[i],q2[i],p3cotPS_Ecm[i],PS[i].real,PS[i].imag))
    print()
        
def print_Lang_Prelovsek_table(file, ld = ""):
    res,  res_sample = read_from_hdf(file)
    mpi = res["mpi"]
    E_n = res["En"]
    p_star = res["pstar"+ld]
    s = res["s"+ld]
    PS = res["PS"+ld]
    p3cotPS_Ecm_ld = res["p3cotPS_Ecm"+ld]
    for i in range(len(res["En"])):
        print("%f\t%f\t%f\t%f\t%f"%(E_n[i],p_star[i],s[i],PS[i].real,PS[i].imag))
    print()


if __name__ == "__main__":
    args = sys.argv
    name = args[1]    
    
    p3_cot_PS(name, x_ax="pstar2", y_ax="p3cotPS", save=True, show=True, ld = "_ld",prime="_prime", pref=name)
    p3_cot_PS(name, x_ax="pstar2", y_ax="p3cotPS", save=True, show=True, prime="_prime", pref=name)
    p3_cot_PS(name, x_ax="sqrt_s", y_ax="sigma", save=True, show=True, ld = "_ld",prime="_prime", pref=name,delete=False)
    p3_cot_PS(name, x_ax="sqrt_s", y_ax="sigma", save=True, show=True, prime="_prime", pref=name,delete=False)
