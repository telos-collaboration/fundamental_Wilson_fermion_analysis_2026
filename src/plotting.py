import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
import sys
import os.path

# def read_from_hdf(filename):
#     res, res_tmp = [{},{}]
#     with h5py.File("output/hdf5/"+filename+".hdf5","r") as hfile:
#         for key in hfile.keys():
#             if key[:4] == "orig":
#                 res[key[5:]] = hfile[key][()]
#             if key[:4] == "samp":
#                 res_tmp[key[7:]] = hfile[key][()]
#     return res, res_tmp

# def error_of_array_gauss(array):
#     tempo = []
#     for i in range(len(array[0])):
#         tmp = array[:,i]
#         tmp.sort()
#         tempo.append(tmp)
#     num_perc = math.erf(1/np.sqrt(2))
#     length = len(array)
#     result = np.transpose(tempo)
#     return result[length//2], abs(result[length//2]-result[math.floor(length*(1-num_perc)/2)]), abs(result[length//2]-result[math.ceil(length*(1+num_perc)/2)]),result[math.floor(length*(1-num_perc)/2)],result[math.ceil(length*(1+num_perc)/2)]

# def error_of_array_lin(arr):
#     length = len(arr)
#     obs = len(arr[0])
#     return arr[length//2], [abs(arr[length//2][i]-arr[0][i]) for i in range(obs)], [abs(arr[length//2][i]-arr[length-1][i]) for i in range(obs)], arr[0], arr[length-1]

# def error_of_array(resampling = "gauss"):
#     resampling = bytes.decode(resampling)
#     if resampling == "gauss":
#         return error_of_array_gauss
#     elif resampling == "lin":
#         return error_of_array_lin

# def en_L(N_L,p2,mpi):                                           # wrong formula for q^2 \elem Z
#     return 1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2)
# def en_rho_L(N_L,p2,mrho,mpi):
#     return np.sqrt(mrho**2+p2*(2*np.pi/N_L)**2)/mpi

# def marker(pind):
#     if pind == 1:
#         return "o"
#     else:
#         return "x"
# def colorf(en_lv_ind):
#     if en_lv_ind == 1:
#         return "blue"
#     else:
#         return "red"

# def color_NL(NL):
#     if NL == 14:
#         return "red"
#     elif NL == 16:
#         return "green"
#     elif NL == 24:
#         return "blue"

# def ls_P(dvec):
#     if list(dvec) == [0,0,1]:
#         return "solid"
#     elif list(dvec) == [1,1,0]:
#         return "dashed"

# def ms_P(dvec):
#     if list(dvec) == [0,0,1]:
#         return "*"
#     elif list(dvec) == [1,1,0]:
#         return "o"

# def ms_size(level):
#     if level == 1:
#         return 8
#     elif level == 2:
#         return 15

# def delete_steps(arr, sign = 1, delete=True):
#     if delete:
#         for i in range(len(arr)-1):
#             if arr[i+1] < sign*arr[i]: 
#                 arr[i] = np.nan
#         return arr
#     else:
#         return arr

# def p3_cot_PS(file, show=False, save = True, pref = "", x_ax = "sqrt_s", y_ax = "p3cotPS_Ecm", ld="", prime = "", delete = True):
#     plt.rcParams['figure.figsize'] = [10, 6]
#     fontsize = 14
#     font = {'size'   : fontsize}
#     matplotlib.rc('font', **font)
#     fig, ax = plt.subplots()
#     res,  res_sample = read_from_hdf(file)
#     num_gaussian = len(res_sample["E_cm_prime"])

#     lvls = res["en_lv"] 
#     N_Ls = res["N_L"]    
#     plt.grid()

#     mpi_pr_s = ""
#     if not prime == "":
#         mpi_pr_s = r"/m_{\pi}"

#     if x_ax == "sqrt_s":
#         plt.xlabel("$E_{CM}"+mpi_pr_s+"$")
#         x_plot = res["E_cm"+ld+prime]
#         x_plot_sam = np.transpose(res_sample["E_cm"+ld+prime])
#         ax.set_xlim([2,2.5])
#     elif x_ax == "s":
#         plt.xlabel("$E_{CM}^2"+mpi_pr_s)                                    # needs fix
#         x_plot = res["s"+ld+prime]
#         x_plot_sam = np.transpose(res_sample["s"+ld+prime])
#         ax.set_xlim([0.1,0.4])
#     elif x_ax == "pstar2":
#         plt.xlabel(r"$p^{\star^2}"+mpi_pr_s+"^2$")
#         x_plot = res["p2star"+ld+prime]
#         x_plot_sam = np.transpose(res_sample["p2star"+ld+prime])
#         ax.set_xlim([0,0.5])
#     elif x_ax == "aE":
#         plt.xlabel("aE")
#         x_plot = res["En"]
#         x_plot_sam = np.transpose(res_sample["En"])
#         ax.set_xlim([0.8, 1.15])
#     elif x_ax == "En":
#         plt.xlabel("$E"+mpi_pr_s+"$")
#         x_plot = res["En"+prime]
#         x_plot_sam = np.transpose(res_sample["En"+prime])
#         ax.set_xlim([2.1, 3])
        
#     if y_ax == "p3cotPS":
#         y_plot = np.real(res["p3cotPS"+ld+prime])
#         y_plot_sam = np.transpose(np.real(res_sample["p3cotPS"+ld+prime]))
#         plt.ylabel(r"$p^3\, \cot(\delta)"+mpi_pr_s+"^3$")    
#         ax.set_ylim([-5,5])
#     elif y_ax == "p3cotPS_Ecm":
#         y_plot = np.real(res["p3cotPS_Ecm"+ld+prime])
#         y_plot_sam = np.transpose(np.real(res_sample["p3cotPS_Ecm"+ld+prime]))
#         plt.ylabel(r"$p^3\, \cot(\delta)/E_CM"+mpi_pr_s+"^3$")    
#         ax.set_ylim([-1,0])
#     elif y_ax == "sigma":
#         y_plot = np.real(res["sigma"+ld+prime])
#         y_plot_sam = np.transpose(np.real(res_sample["sigma"+ld+prime]))
#         plt.ylabel(r"$\sigma_1*"+mpi_pr_s+"^2$")    
#     elif y_ax == "cot_PS":
#         y_plot = np.real(res["cot_PS"+ld])
#         y_plot_sam = np.transpose(np.real(res_sample["cot_PS"+ld]))
#         plt.ylabel(r"cot($\delta$)")      
#         ax.set_ylim([-100,100])
#     elif y_ax == "PS":
#         y_plot = np.real(res["PS+ld"])
#         y_plot_sam = np.transpose(np.real(res_sample["PS"+ld]))
#         plt.ylabel("$q^2$")      
#         ax.set_ylim([0,180])
#     elif y_ax == "q2":
#         y_plot = np.real(res["q2"+ld])
#         y_plot_sam = np.transpose(np.real(res_sample["q2"+ld]))
#         plt.ylabel("$q^2$")    

#     N_Ls = res["N_L"]
#     dvecs = res["dvec"]
#     d2s = res["d2"]

#     length = len(x_plot_sam[0])

#     num_perc = math.erf(1/np.sqrt(2))
#     for i in range(len(N_Ls)):
#         if lvls[i] != 1:
#             ax.scatter(x_plot[i],y_plot[i], color = color_NL(N_Ls[i]), label = "Lv: %i, |P|^2=%i, NL=%i"%(lvls[i],d2s[i],N_Ls[i]), marker = ms_P(dvecs[i]), s=10*ms_size(lvls[i]))
#             if bytes.decode(res["resampling"]) == "gauss":
#                 sorted_indices = np.argsort(x_plot_sam[i])  
#                 ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_plot_sam[i][sorted_indices],delete=delete)[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))
#             elif bytes.decode(res["resampling"]) == "lin":
#                 ax.plot(x_plot_sam[i],delete_steps(y_plot_sam[i],delete=delete), color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))

#     ax.legend(loc="best")
#     if save:
#         fig.savefig("output/plots/%s_%s"%(y_ax,x_ax)+ld+prime+"_"+pref+".pdf", bbox_inches='tight')
#     if show:
#         plt.show()
#     plt.close(fig)

# def print_plymouth_table(file, ld = "_ld"):
#     res,  res_sample = read_from_hdf(file)
#     mpi = res["mpi"]
#     E_n = res["En"]
#     E_cm = res["E_cm"+ld]
#     q2 = res["q2"+ld]
#     p3cotPS_Ecm = res["p3cotPS_Ecm"+ld]
#     PS = res["PS"+ld]
#     for i in range(len(E_n)):
#         print("%f\t%f\t%f\t%f\t%f\t%f"%(E_n[i],E_cm[i],q2[i],p3cotPS_Ecm[i],PS[i].real,PS[i].imag))
#     print()
        
# def print_Lang_Prelovsek_table(file, ld = ""):
#     res,  res_sample = read_from_hdf(file)
#     mpi = res["mpi"]
#     E_n = res["En"]
#     p_star = res["pstar"+ld]
#     s = res["s"+ld]
#     PS = res["PS"+ld]
#     p3cotPS_Ecm_ld = res["p3cotPS_Ecm"+ld]
#     for i in range(len(res["En"])):
#         print("%f\t%f\t%f\t%f\t%f"%(E_n[i],p_star[i],s[i],PS[i].real,PS[i].imag))
#     print()

def color(d2):
    colors = ["orange", "green", "blueviolet"]
    return colors[d2-1]

def E_pipi(mpi,p12,p22,L):
    return np.sqrt(mpi**2+(2*np.pi/L)**2*p12)+np.sqrt(mpi**2+(2*np.pi/L)**2*p22)  

def E_rho(mrho,p2,L):
    return np.sqrt(mrho**2+(2*np.pi/L)**2*p2)

def get_data_E_L(name, beta, m0, num_lv = 2):
    
    NLs = []
    NL_invs = []
    aEs = []
    aE_ms = []
    aE_ps = []
    d2s = []
    lvs = []

    with h5py.File("../output/scattering/isospin1_scattering"+name+".hdf5","r") as hfile:
        hfile.visit(print)
        for key in hfile:
            print(str(beta))
            if str(beta) in key and str(m0) in key:
                print("cool")
                # print(key)
                for P in hfile[key]:
                    if P[0] == "p":
                        dvec = [int(P[2]),int(P[4]),int(P[6])]
                        d2 = np.dot(dvec,dvec)
                        for irrep in hfile[key][P]:
                            mpi = hfile[key][P][irrep]["lv0"]["info"]["mpi"][()]
                            mrho = hfile[key][P][irrep]["lv0"]["info"]["mrho"][()]
                            for lv in range(num_lv):
                                NL = int(hfile[key][P][irrep]["lv0"]["info"]["NL"][()])
                                NLs.append(NL)
                                NL_invs.append(1/NL)
                                aEs.append(hfile[key][P][irrep]["E%i"%lv][()][0])
                                aE_ms.append(hfile[key][P][irrep]["Delta_E%i"%lv][()][0])
                                aE_ps.append(hfile[key][P][irrep]["Delta_E%i"%lv][()][0])
                                d2s.append(d2)
                                lvs.append(lv)
    return mpi, mrho, d2s, NLs, NL_invs, aEs, aE_ms, aE_ps, lvs

def LminLmax(m0):                   # maybe to be replaced by input file
    if m0 == -0.92:
        return [12,26]
    elif m0 == -0.863:
        return [14,38]
    elif m0 == -0.867:
        return [14,38]
    
def marker(lv):                     # maybe to be replaced by input file
    markers = ["*", "o", "^"]
    return markers[lv-1]

def plot_E_L(name,beta,m0,levels=False,outname=None):
    mpi, mrho, d2s, NLs, NL_invs, En, En_m_err, En_p_err, lvs = get_data_E_L(name, beta, m0)
    
    Lmin, Lmax = LminLmax(m0)
    # Lmin, Lmax = [0.3,2]
    
    # d2s = [x[0]**2+x[1]**2+x[2]**2 for x in dvecs] 

    for i in range(len(En)):
        plt.errorbar([NL_invs[i],],y=[En[i],],yerr=[[En_m_err[i],],[En_p_err[i],]], solid_capstyle="projecting", capsize=5, ls="", color = color(d2s[i]), marker = marker(lvs[i]))   
    # plt.axhline(mpi,c="black", ls="dotted", label = "$m_\pi$")
    plt.axhline(mrho,c="red", ls="dotted", label = "$m_\\rho$")
    plt.axhline(2*mpi,c="black",label = "2$m_\pi$")
    plt.axhline(4*mpi,c="black",label = "4$m_\pi$")
    plt.grid()
    plt.title("$\\beta$ = %f, $m_0$ = %f"%(beta,m0))
    # xarr = np.linspace(14,38)
    xarrinv = np.linspace(1/40,1/13)
    xarr = [1/x for x in xarrinv]
    if levels:
        yarr1_1 = [E_rho(mrho,1,x) for x in xarr]
        yarr1_2 = [E_pipi(mpi,1,0,x) for x in xarr]
        yarr1_3 = [E_pipi(mpi,2,1,x) for x in xarr]
        yarr1_4 = [E_pipi(mpi,3,2,x) for x in xarr]
        yarr1_5 = [E_pipi(mpi,4,1,x) for x in xarr]
        yarr2_1 = [E_rho(mrho,2,x) for x in xarr]
        yarr2_2 = [E_pipi(mpi,2,0,x) for x in xarr]
        yarr2_3 = [E_pipi(mpi,3,1,x) for x in xarr]
        yarr2_4 = [E_pipi(mpi,4,2,x) for x in xarr]
        yarr3_1 = [E_rho(mrho,3,x) for x in xarr]
        yarr3_2 = [E_pipi(mpi,3,0,x) for x in xarr]
        yarr3_3 = [E_pipi(mpi,2,1,x) for x in xarr]
        yarr3_4 = [E_pipi(mpi,4,3,x) for x in xarr]
        plt.plot(xarrinv,yarr1_1, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_2, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_3, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_4, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_5, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr2_1, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr2_2, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr2_3, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr2_4, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr3_1, ls="solid", c=color(3))
        plt.plot(xarrinv,yarr3_2, ls="solid", c=color(3))
        plt.plot(xarrinv,yarr3_3, ls="solid", c=color(3))
        plt.plot(xarrinv,yarr3_4, ls="solid", c=color(3))

    # yarr_low_1 = [np.sqrt(4*mpi**2+(2*np.pi/x)**2*1) for x in xarr]
    # yarr_low_2 = [np.sqrt(4*mpi**2+(2*np.pi/x)**2*2) for x in xarr]
    # yarr_low_3 = [np.sqrt(4*mpi**2+(2*np.pi/x)**2*3) for x in xarr]
    # yarr_hig_1 = [np.sqrt(16*mpi**2+(2*np.pi/x)**2*1) for x in xarr]
    # yarr_hig_2 = [np.sqrt(16*mpi**2+(2*np.pi/x)**2*2) for x in xarr]
    # yarr_hig_3 = [np.sqrt(16*mpi**2+(2*np.pi/x)**2*3) for x in xarr]

    # if ranges:
    #     plt.fill_between(xarrinv,yarr_low_1,yarr_hig_1, color=color(1), alpha = 0.3)
    #     plt.fill_between(xarrinv,yarr_low_2,yarr_hig_2, color=color(2), alpha = 0.3)
    #     plt.fill_between(xarrinv,yarr_low_3,yarr_hig_3, color=color(3), alpha = 0.3)

    # plt.plot([0,0],[0,0],c="grey", label = "non-int for \nP={1,2,3}")
    plt.plot([0,0],[0,0],c="grey", label = "non-int")
    plt.xlim([1/40,1/13])
    plt.ylim([0.3,2])

    for i in range(1,4):
        plt.errorbar([0,],y=[0,],yerr=[[0,],[0,]], solid_capstyle="projecting", capsize=5, ls="", color = color(i), marker = "o", label = "|P|=%i"%(i))
    for i in range(1,4):
        plt.scatter(x=[0,],y=[0,], color = "grey", marker = marker(i), label = "lv=%i"%(i))

    plt.legend(loc='center right', bbox_to_anchor=(1.24, 0.5))

    plt.xticks([1/14,1/16,1/20,1/24,1/36],["1/14","1/16","1/20","1/24","1/36"])
    plt.xlabel("1/$N_L$")
    plt.ylabel("a$E$")
    # plt.show()
    if outname == None:    
        plt.savefig("../output/plots/scattering/E_L_b%f_m0%f_levels_%r.pdf"%(beta,m0,levels), bbox_inches='tight')
    else:    
        plt.savefig("../output/plots/scattering/"+outname+"_levels_%r.pdf"%levels, bbox_inches='tight')
    plt.clf()


if __name__ == "__main__":
    plot_E_L("_evp_deriv_false",6.9,-0.92,True,outname="non_res")
    plot_E_L("_evp_deriv_false",6.9,-0.92,False,outname="non_res")
    # args = sys.argv
    # name = args[1]
    
    # p3_cot_PS(name, x_ax="pstar2", y_ax="p3cotPS", save=True, show=plt.isinteractive(), ld = "_ld",prime="_prime", pref=name)
    # p3_cot_PS(name, x_ax="pstar2", y_ax="p3cotPS", save=True, show=plt.isinteractive(), prime="_prime", pref=name)
    # p3_cot_PS(name, x_ax="sqrt_s", y_ax="sigma",   save=True, show=plt.isinteractive(), ld = "_ld",prime="_prime", pref=name,delete=False)
    # p3_cot_PS(name, x_ax="sqrt_s", y_ax="sigma",   save=True, show=plt.isinteractive(), prime="_prime", pref=name,delete=False)