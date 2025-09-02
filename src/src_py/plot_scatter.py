import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
import fit_scatter
import os.path as op
import os
import sys

num_perc = math.erf(1/np.sqrt(2))

########################################## Plot energy levels in lattice units ##########################################

def color(level):
    colors = ["green","orange", "cyan", "blueviolet"]
    return colors[level]
    
def marker_NL_lv(NL, lv):                     # maybe to be replaced by input file
    if NL == 14:
        if lv == 0:
            return "o"
        elif lv == 1:
            return "s"
    elif NL == 16:
        if lv == 0:
            return "*"
        elif lv == 1:
            return "P"
    elif NL == 20:
        if lv == 0:
            return "h"
        elif lv == 1:
            return "p"
    elif NL == 24:
        if lv == 0:
            return "<"
        elif lv == 1:
            return ">"
    elif NL == 36:
        if lv == 0:
            return "^"
        elif lv == 1:
            return "v"
        raise RuntimeError("Wrong NL or lv in marker: %i, %i"%(NL, lv))

def color_d_irrep(d, irrep):
    if d == 0:
        if irrep == "T1":
            return "sienna"
    elif d == 1:
        if irrep == "A1":
            return "green"
        elif irrep == "E":
            return "darkgreen"
    elif d == 2:
        if irrep == "A1":
            return "blue"
        elif irrep == "B1":
            return "darkblue"
    elif d == 3:
        if irrep == "A1":
            return "red"
        elif irrep == "E":
            return "darkred"
    raise ValueError("wrong d or irrep in color_d_irrep(): %i, %s"%(d,irrep))

def ls_NL(NL):
    if NL == 14:
        return "solid"
    elif NL == 16:
        return (0,(1,1))
    elif NL == 20:
        return "dashed"
    elif NL == 24:
        return "dashdot"
    elif NL == 36:
        return "dotted"
    else:
        raise ValueError("Wrong NL given to ls_NL()")

def E_pipi(mpi,p12,p22,L):
    return np.sqrt(mpi**2+(2*np.pi/L)**2*p12)+np.sqrt(mpi**2+(2*np.pi/L)**2*p22)  

def E_rho(mrho,p2,L):
    return np.sqrt(mrho**2+(2*np.pi/L)**2*p2)

########################################## Plot com energies unitless ##########################################


def get_data_E_CM_L(h5file_fit, beta, m0):
    NLs,NL_invs,aEs,aE_ms,aE_ps,d2s,lvs,irreps = [[],[],[],[],[],[],[],[]]

    with h5py.File(h5file_fit,"r") as hfile:
        for key in hfile:
            if str(beta) in key and str(m0) in key:
                for P in hfile[key]:
                    if P[0] == "p":
                        dvec = [int(P[2]),int(P[4]),int(P[6])]
                        d2 = np.dot(dvec,dvec)
                        for irrep in hfile[key][P]:
                            if irrep != "pi":
                                mpi = hfile[key][P][irrep]["lv0"]["info"]["mpi"][()]
                                mrho = hfile[key][P][irrep]["lv0"]["info"]["mrho"][()]
                                for lv in hfile[key][P][irrep]:
                                    # print(hfile[key][P][irrep].keys())
                                    if lv[:2] == "lv":
                                        # print(key+P+irrep,lv)
                                        i = int(lv[2:])
                                        NL = int(hfile[key][P][irrep][lv]["info"]["NL"][()])
                                        NLs.append(NL)
                                        aE = hfile[key][P][irrep][lv]["mean"]["E_cm_prime"][()]
                                        aEs.append(aE)
                                        # print(aE)
                                        # print()
                                        E_spl = sorted(hfile[key][P][irrep][lv]["sample"]["E_cm_prime"][()])
                                        length = len(E_spl)
                                        aE_ms.append(abs(aE-E_spl[math.floor(length*(1-num_perc)/2)]))
                                        aE_ps.append(abs(aE-E_spl[math.ceil(length*(1+num_perc)/2)]))
                                        NL_invs.append(1/NL)
                                        d2s.append(d2)
                                        lvs.append(i)
                                        irreps.append(irrep)
    return mpi, mrho, d2s, NLs, NL_invs, aEs, aE_ms, aE_ps, lvs, irreps

def plot_E_CM_L(h5file_scatter,beta,m0,levels=False,outname=None,show=False):
    
    mpi, mrho, d2s, NLs, NL_invs, ECMs, ECM_errms, ECM_errps, lvs, irreps = get_data_E_CM_L(h5file_scatter, beta, m0)
    

    for i in range(len(ECMs)):
        ECM_errms[i] = 0 if ECM_errms[i] > 1 else ECM_errms[i]
        ECM_errps[i] = 0 if ECM_errps[i] > 1 else ECM_errps[i]
        # plt.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, ls="", color = color(lvs[i]), marker = marker(d2s[i],irreps[i]))   
        plt.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, ls="", color = color_d_irrep(d2s[i],irreps[i]), marker = marker_NL_lv(NLs[i],lvs[i]))   
    plt.axhline(1,c="black", ls="dotted", label = r"$m_\pi$")
    plt.axhline(mrho/mpi,c="red", ls="dotted", label = r"$m_\rho$")
    plt.axhline(2,c="black",label = r"2$m_\pi$")
    plt.axhline(4,c="black",label = r"4$m_\pi$")
    plt.grid()
    plt.title("$\\beta$ = %f, $m_0$ = %f"%(beta,m0))
    xarrinv = np.linspace(1/40,1/13)
    xarr = [1/x for x in xarrinv]
    if levels:
        yarr1_2 = [np.sqrt(E_pipi(mpi,1,0,x)**2-(2*np.pi/x)**2*1)/mpi for x in xarr]
        yarr1_3 = [np.sqrt(E_pipi(mpi,2,1,x)**2-(2*np.pi/x)**2*1)/mpi for x in xarr]
        yarr1_4 = [np.sqrt(E_pipi(mpi,3,2,x)**2-(2*np.pi/x)**2*1)/mpi for x in xarr]
        yarr1_5 = [np.sqrt(E_pipi(mpi,4,1,x)**2-(2*np.pi/x)**2*1)/mpi for x in xarr]
        yarr2_2 = [np.sqrt(E_pipi(mpi,2,0,x)**2-(2*np.pi/x)**2*2)/mpi for x in xarr]
        yarr2_3 = [np.sqrt(E_pipi(mpi,3,1,x)**2-(2*np.pi/x)**2*2)/mpi for x in xarr]
        yarr2_4 = [np.sqrt(E_pipi(mpi,4,2,x)**2-(2*np.pi/x)**2*2)/mpi for x in xarr]
        yarr3_2 = [np.sqrt(E_pipi(mpi,3,0,x)**2-(2*np.pi/x)**2*3)/mpi for x in xarr]
        yarr3_3 = [np.sqrt(E_pipi(mpi,2,1,x)**2-(2*np.pi/x)**2*3)/mpi for x in xarr]
        yarr3_4 = [np.sqrt(E_pipi(mpi,4,3,x)**2-(2*np.pi/x)**2*3)/mpi for x in xarr]
        plt.plot(xarrinv,yarr1_2, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_3, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_4, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr1_5, ls="dashed", c=color(1))
        plt.plot(xarrinv,yarr2_2, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr2_3, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr2_4, ls="dashdot", c=color(2))
        plt.plot(xarrinv,yarr3_2, ls="solid", c=color(3))
        plt.plot(xarrinv,yarr3_3, ls="solid", c=color(3))
        plt.plot(xarrinv,yarr3_4, ls="solid", c=color(3))
    plt.plot([0,0],[0,0],c="grey", label = "non-int")
    plt.xlim([1/40,1/13])
    plt.ylim([1,5])

    for tmp in [[0,"T1"],[1,"A1"],[1,"E"],[2,"A1"],[2,"B1"],[3,"A1"],[3,"E"]]:
        plt.scatter(x=[-100,],y=[-100,], color = color_d_irrep(tmp[0],tmp[1]), marker = "o", label = "p=%i, %s"%(tmp[0],tmp[1]))
    for tmp in [[14,0],[14,1],[16,0],[16,1],[24,0],[24,1],[36,0],[36,1]]:
        plt.scatter(x=[-100,],y=[-100,], color = "grey", marker = marker_NL_lv(tmp[0],tmp[1]), label = "$N_L$=%i, lv%s"%(tmp[0],tmp[1]))

    # for i in range(1,4):
    #     plt.errorbar([0,],y=[0,],yerr=[[0,],[0,]], solid_capstyle="projecting", capsize=5, ls="", color = color(i), marker = "", label = "|P|=%i"%(i))

    plt.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))

    plt.xlabel("1/$N_L$")
    plt.ylabel("$E_{CM}$/$m_\\pi$")
    plt.xticks([1/14,1/16,1/20,1/24,1/36],["1/14","1/16","1/20","1/24","1/36"])
    if outname == None:    
        plt.savefig(op.join(PLTDIR, "E_CM_L_b%f_m0%f_levels_%r.pdf"%(beta,m0,levels)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR, "E_CM_L_"+outname+"_levels_%r.pdf"%levels), bbox_inches='tight')
    if show:
        plt.show()
    plt.clf()

########################################## Plot p3cotPS ##########################################

def delete_steps(arr, sign = 1, delete=False):
    # for i in range(1,len(arr)-1):
    #     if abs(arr[i]-arr[i+1]) > 10*abs(arr[i-1]-arr[i]):
    #         arr[i] = np.nan
    for i in range(1,len(arr)-1):
        if abs(arr[i] > 1):
            if np.sign(arr[i]) != np.sign(arr[i+1]):
                arr[i-1] = np.nan
                arr[i] = np.nan
                arr[i+1] = np.nan
    for i in range(len(arr)):
        if arr[i] == 0:
            arr[i] = np.nan
    return arr
    # if delete:
    #     for i in range(len(arr)-1):
    #         if arr[i+1] < sign*arr[i]: 
    #             arr[i] = np.nan
    #     return arr
    # else:
    #     return arr

def get_data_p3cotPS(h5file_scatter_fit, beta, m0):
    fit_param_mean = {}
    fit_param_spl = {}
    scat_fit_mean = {}
    scat_fit_spl = {}
    info = {}
    # scat_not_fit_mean = {}
    # scat_not_fit_spl = {}
    # fit_beta_m = "fit_b%f_m%f"%(beta,m0)
    with h5py.File(h5file_scatter_fit,"r") as hfile:
        for ens in hfile:
            if str(beta) in ens and str(m0) in ens:
                # for key in hfile[fit_beta_m]["mean"]:
                #     fit_param_mean[key] = hfile[fit_beta_m]["mean"][key][()]
                #     fit_param_spl[key] = hfile[fit_beta_m]["sample"][key][()]
                for P in hfile[ens]:
                    if P[0] == "p":
                        for irrep in hfile[ens][P]:
                            for lv in hfile[ens][P][irrep]:
                                if lv[:2] == "lv":
                                    info.setdefault("lv",[]).append(int(lv[2:]))
                                    info.setdefault("irrep",[]).append(irrep)
                                    # print(ens+P+irrep+lv)
                                    # if hfile[ens][P][irrep][lv]["fit"][()]:
                                    for key in hfile[ens][P][irrep][lv]["mean"]:
                                        if key == "dvec" or key == "N_L":
                                            scat_fit_mean.setdefault(key,[]).append(hfile[ens][P][irrep][lv]["mean"][key][()])
                                        else:
                                            scat_fit_mean.setdefault(key,[]).append(float(np.real(hfile[ens][P][irrep][lv]["mean"][key][()])))
                                            scat_fit_spl.setdefault(key,[]).append([float(x) for x in np.real(hfile[ens][P][irrep][lv]["sample"][key][()])])
                                    # else:
                                    #     scat_not_fit_mean.setdefault("irrep",[]).append(irrep)
                                    #     for key in hfile[ens][P][irrep][lv]["mean"]:
                                    #         if key == "dvec" or key == "N_L":
                                    #             scat_not_fit_mean.setdefault(key,[]).append(hfile[ens][P][irrep][lv]["mean"][key][()])
                                    #         else:
                                    #             scat_not_fit_mean.setdefault(key,[]).append(float(np.real(hfile[ens][P][irrep][lv]["mean"][key][()])))
                                    #             scat_not_fit_spl.setdefault(key,[]).append([float(x) for x in np.real(hfile[ens][P][irrep][lv]["sample"][key][()])])
    return info, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl# , scat_not_fit_mean, scat_not_fit_spl

def plot_p3cotPS(h5file_scatter_fit,beta,m0,fit=False,outname=None,show=False):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    plt.grid()
    info, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl = get_data_p3cotPS(h5file_scatter_fit, beta, m0)
    xlim = [0,0.3] if beta==6.9 else [0,3]
    ax.set_xlim(xlim)
    # ylim = [-1,1] if beta==6.9 else [-4,4]
    ylim = [-4,4]
    ax.set_ylim(ylim)

    plt.xlabel(r"$p^{\star^2}/m_\pi^2$")
    x_plot       = np.asarray(scat_fit_mean["p2star_prime"])
    x_plot_sam   = np.asarray(scat_fit_spl["p2star_prime"])
    y_plot       = np.asarray(scat_fit_mean["p3cotPS_prime"])
    y_plot_sam   = np.asarray(scat_fit_spl["p3cotPS_prime"])
    # x_n_plot     = np.asarray(scat_not_fit_mean["p2star_prime"])                                        # n marks that it was not fitted
    # x_n_plot_sam = np.asarray(scat_not_fit_spl["p2star_prime"])
    # y_n_plot     = np.asarray(scat_not_fit_mean["p3cotPS_prime"])
    # y_n_plot_sam = np.asarray(scat_not_fit_spl["p3cotPS_prime"])
    plt.ylabel(r"$p^3\, \cot(\delta)/m_\pi^3$")
    
    length = len(x_plot_sam[0])
    
    N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    
    # N_L_ns = scat_not_fit_mean["N_L"]
    # dvec_ns = scat_not_fit_mean["dvec"]
    # dvec_ns = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvec_ns]
    # d2_ns = [np.dot(d,d) for d in dvec_ns]

    lvs = info["lv"]
    irreps = info["irrep"]

    for i in  range(len(x_plot)):
        # if 0<x_plot[i]<3: 
        if y_plot[i] != 0:
            ax.scatter(x_plot[i],y_plot[i], color = color_d_irrep(d2s[i],irreps[i]), ls = ls_NL(N_Ls[i]), marker = marker_NL_lv(N_Ls[i],lvs[i]),s=60)   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
            sorted_indices = np.argsort(x_plot_sam[i])
            ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_plot_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color_d_irrep(d2s[i],irreps[i]), ls = ls_NL(N_Ls[i]))
        # else:
        #     raise ValueError("Fitted momentum is not in elastic threshhold in plotting.py!!!")

    # for i in  range(len(x_n_plot)):
    #     ax.scatter(x_n_plot[i],y_n_plot[i], color = "grey", ls = ls_P(dvec_ns[i]), marker = marker(scat_not_fit_mean["irrep"][i]),s=60)
    #     sorted_indices = np.argsort(x_n_plot_sam[i])
    #     ax.plot(x_n_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_n_plot_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = "grey", ls = ls_P(dvec_ns[i]))

    # xarr = np.linspace(xlim[0], xlim[1])
    
    # if fit:
    #     a1_1 = fit_param_mean["a1_1"]
    #     r1_1 = fit_param_mean["r1_1"]

    #     a1_1_smp = fit_param_spl["a1_1"]
    #     r1_1_smp = fit_param_spl["r1_1"]

    #     yarr = [fit_scatter.ERE_1(x,a1_1,r1_1) for x in xarr]
    #     yarr_smp = [sorted([fit_scatter.ERE_1(x,a1_1_smp[i],r1_1_smp[i]) for i in range(len(a1_1_smp))]) for x in xarr]

    #     yarr_m = [yarr_smp[i][math.floor(length*(1-num_perc)/2)] for i in range(len(yarr_smp))]
    #     yarr_p = [yarr_smp[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(yarr_smp))]

    #     plt.plot(xarr,yarr, color = "blue")
    #     plt.fill_between(xarr, yarr_m, yarr_p, alpha = 0.3, color = "blue")

    for tmp in [[0,"T1"],[1,"A1"],[1,"E"],[2,"A1"],[2,"B1"],[3,"A1"],[3,"E"]]:
        plt.scatter(x=[-100,],y=[-100,], color = color_d_irrep(tmp[0],tmp[1]), marker = "o", label = "p=%i, %s"%(tmp[0],tmp[1]))
    for tmp in [[14,0],[14,1],[16,0],[16,1],[24,0],[24,1],[36,0],[36,1]]:
        plt.scatter(x=[-100,],y=[-100,], color = "grey", marker = marker_NL_lv(tmp[0],tmp[1]), label = "$N_L$=%i, lv%s"%(tmp[0],tmp[1]))

    ax.legend(loc='center right', bbox_to_anchor=(1.35, 0.5))
    if outname == None:    
        plt.savefig(op.join(PLTDIR, "p3cotPS_b%f_m0%f_%r.pdf"%(beta,m0,fit)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR, "p3cotPS_"+outname+"_%r.pdf"%fit), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def p2_s_prime(s):
    return s/4-1

if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file_scatter  = args[2]

    os.makedirs(PLTDIR, exist_ok=True)

    # plot_E_CM_L(h5file_scatter,6.9,-0.92,False,outname="non_res",show=False)
    plot_E_CM_L(h5file_scatter,6.9,-0.92,True,outname="non_res")
    # plot_E_CM_L(h5file_scatter,7.05,-0.863,False,outname="close_res",show=False)
    plot_E_CM_L(h5file_scatter,7.05,-0.863,True,outname="close_res")
    # plot_E_CM_L(h5file_scatter,7.05,-0.867,False,outname="res",show=False)
    plot_E_CM_L(h5file_scatter,7.05,-0.867,True,outname="res",show=False)
    
    plot_p3cotPS(h5file_scatter,6.9,-0.92,False,outname="non_res",show=True)
    plot_p3cotPS(h5file_scatter,7.05,-0.863,False,outname="close_res",show=True)
    plot_p3cotPS(h5file_scatter,7.05,-0.867,False,outname="res",show=True)