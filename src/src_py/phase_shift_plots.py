import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import math
# import fit_scatter
import os.path as op
import os
import sys
import plotting_functions as pf
import fit_models as fm

import styles

mpl.rcParams['lines.markersize'] = 9

figs1, figs2 = 12,7

plt.rcParams['figure.figsize'] = [10, 6] 
fontsize = 14
font = {'size'   : fontsize}
matplotlib.rc('font', **font)
plt.rcParams.update({
    # "font.family": "serif",
    "mathtext.fontset": "cm",   # Computer Modern
})

num_perc = math.erf(1/np.sqrt(2))

def delete_steps(arr, sign = 1, delete=False):
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

def get_data(h5file_scatter_fit, beta, m0, fit):                # wont work with current get_data()
    fit_param_mean = {}
    fit_param_spl = {}
    scat_fit_mean = {}
    scat_fit_spl = {}
    info = {}
    scat_nf_mean = {}
    scat_nf_spl = {}
    info_nf = {}
    fit_beta_m = "fit_b%f_m%f"%(beta,m0)
    with h5py.File(h5file_scatter_fit,"r") as hfile:
        for ens in hfile:
            if str(beta) in ens and str(m0) in ens:
                if fit:
                    for key in hfile[fit_beta_m]["mean"]:
                        fit_param_mean[key] = hfile[fit_beta_m]["mean"][key][()]
                        fit_param_spl[key] = hfile[fit_beta_m]["sample"][key][()]
                for P in hfile[ens]:
                    if P[0] == "p":
                        for irrep in hfile[ens][P]:
                            for lv in hfile[ens][P][irrep]:
                                if lv[:2] == "lv":
                                    fit_this = False
                                    if fit:
                                        fit_this = hfile[ens][P][irrep][lv]["fit"][()]
                                    if fit_this or not fit:
                                        info.setdefault("lv",[]).append(int(lv[2:]))
                                        info.setdefault("irrep",[]).append(irrep)
                                        for key in hfile[ens][P][irrep][lv]["info"]:
                                            info.setdefault(key, []).append(hfile[ens][P][irrep][lv]["info"][key][()])
                                        scat_fit_mean.setdefault("irrep",[]).append(irrep)
                                        for key in hfile[ens][P][irrep][lv]["mean"]:
                                            if key == "dvec" or key == "N_L":
                                                scat_fit_mean.setdefault(key,[]).append(hfile[ens][P][irrep][lv]["mean"][key][()])
                                            else:
                                                scat_fit_mean.setdefault(key,[]).append(float(np.real(hfile[ens][P][irrep][lv]["mean"][key][()])))
                                                scat_fit_spl.setdefault(key,[]).append([float(x) for x in np.real(hfile[ens][P][irrep][lv]["sample"][key][()])])
                                    else:
                                        info_nf.setdefault("lv",[]).append(int(lv[2:]))
                                        info_nf.setdefault("irrep",[]).append(irrep)
                                        for key in hfile[ens][P][irrep][lv]["info"]:
                                            info_nf.setdefault(key, []).append(hfile[ens][P][irrep][lv]["info"][key][()])
                                        scat_nf_mean.setdefault("irrep",[]).append(irrep)
                                        for key in hfile[ens][P][irrep][lv]["mean"]:
                                            if key == "dvec" or key == "N_L":
                                                scat_nf_mean.setdefault(key,[]).append(hfile[ens][P][irrep][lv]["mean"][key][()])
                                            else:
                                                scat_nf_mean.setdefault(key,[]).append(float(np.real(hfile[ens][P][irrep][lv]["mean"][key][()])))
                                                scat_nf_spl.setdefault(key,[]).append([float(x) for x in np.real(hfile[ens][P][irrep][lv]["sample"][key][()])])
    return info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl

def s_p2(p2):
    return 4*(1+p2)

def p2_s(s):
    return s/4-1

def delta_x(x,p2):
    return np.arctan(p2**(3/2)/x)*360/(2*np.pi)

def delta_res_x(x,p2):
    return np.arctan(p2**(3/2)/(2*np.sqrt(1+p2)*x))*360/(2*np.pi)

def plot_PS_ERE_non_res(h5file,show=False):
    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(figs1,figs2))
    plt.subplots_adjust(wspace=0, hspace=0.05)   

    ax1.xaxis.grid()
    ax2.xaxis.grid()
    ax1.yaxis.grid()
    ax2.yaxis.grid()
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, 6.9, -0.92, True)

    slow, shigh = 4,5.07

    x_s   = np.asarray(scat_fit_spl["s_prime"])
    PS_m   = np.asarray(scat_fit_mean["PS"])
    PS_s   = np.asarray(scat_fit_spl["PS"])
    p3cotPS_m   = np.asarray(scat_fit_mean["p3cotPS_prime"])
    p3cotPS_s   = np.asarray(scat_fit_spl["p3cotPS_prime"])

    ax2.set_xlabel(r"$s/m_\pi^2$", fontsize=styles.fontsize)
    ax2.set_ylabel(r"$p^3\, \cot(\delta_1)/m_\pi^3$", fontsize=styles.fontsize)
    ax1.set_ylabel(r"$\delta_1$", fontsize=styles.fontsize)
    
    length = len(x_s[0])
    
    N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(N_Ls,d2s,irreps,lvs))
    
    sarr = np.linspace(slow,shigh,100)
    p2arr = [p2_s(s) for s in sarr]
    fit_model = fm.ERE_0_model
    
    fit_param_m = np.asarray([fit_param_mean[fp] for fp in fit_model.param_names])
    yarr_m = np.asarray([fit_model.model(x,*fit_param_m) for x in p2arr])

    fit_param_s = np.transpose(np.asarray([fit_param_spl[fp] for fp in fit_model.param_names]))
    yarr_tmp = np.asarray([sorted([fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))]) for x in p2arr])

    yarr_e_m_plot = np.asarray([yarr_tmp[i][math.floor(length*(1-num_perc)/2)] for i in range(len(p2arr))])
    yarr_e_p_plot = np.asarray([yarr_tmp[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(p2arr))])

    ax2.plot(sarr,yarr_m, color = styles.c_10_non_res)
    ax2.fill_between(sarr, yarr_e_m_plot, yarr_e_p_plot, alpha = 0.5, color = styles.c_10_non_res)

    PS_m_plot = [delta_x(yarr_m[i], p2arr[i]) for i in range(len(p2arr))]
    PS_e_m_plot = [delta_x(yarr_e_m_plot[i], p2arr[i]) for i in range(len(p2arr))]
    PS_e_p_plot = [delta_x(yarr_e_p_plot[i], p2arr[i]) for i in range(len(p2arr))]

    ax1.plot(sarr,PS_m_plot, color = styles.c_10_non_res)
    ax1.fill_between(sarr, PS_e_m_plot, PS_e_p_plot, alpha = 0.5, color = styles.c_10_non_res)

    x_m = [sorted(x_s[i])[length//2-1] for i in range(len(x_s))]

    for i in  range(len(x_m)):
        ax1.scatter(x_m[i],PS_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax1.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(PS_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))

    for i in  range(len(x_m)):
        ax2.scatter(x_m[i],p3cotPS_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax2.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(p3cotPS_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))
   

    for tmp in [[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
        ax1.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "^", label = r"$|p|=%i$"%(tmp[1]))
    # ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
    ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")
    # ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")

    ax1.set_xlim([slow,shigh])
    ax1.set_ylim([0,50])
    ax2.set_ylim([0,0.55])

    xticks = np.linspace(4,5,6)
    ax1.set_xticks(xticks, [r"$%1.1f$"%x for x in xticks])
    yticks = np.linspace(0,50,6)
    ax1.set_yticks(yticks, [r"$%i$"%x for x in yticks])
    yticks = np.linspace(0,0.5,6)
    ax2.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    
    ax1.legend(loc='upper left', fontsize=styles.fontsize)

    plt.savefig(op.join(PLTDIR, "phase_shift_plot_non_res.pdf"), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def plot_PS_ERE_close_res(h5file,show=False):
    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(figs1,figs2))
    plt.subplots_adjust(wspace=0, hspace=0.05)   

    ax1.xaxis.grid()
    ax2.xaxis.grid()
    ax1.yaxis.grid()
    ax2.yaxis.grid()
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, 7.05, -0.863, False)

    x_s   = np.asarray(scat_fit_spl["s_prime"])
    PS_m   = np.asarray(scat_fit_mean["PS"])
    PS_s   = np.asarray(scat_fit_spl["PS"])
    p3cotPS_m   = np.asarray(scat_fit_mean["p3cotPS_prime"])
    p3cotPS_s   = np.asarray(scat_fit_spl["p3cotPS_prime"])

    ax2.set_xlabel(r"$s/m_\pi^2$", fontsize=styles.fontsize)
    ax2.set_ylabel(r"$p^3\, \cot(\delta_1)/m_\pi^3$", fontsize=styles.fontsize)
    ax1.set_ylabel(r"$\delta_1$", fontsize=styles.fontsize)
    
    length = len(x_s[0])
    
    N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(N_Ls,d2s,irreps,lvs))

    ax1.set_xlim([4,15])
    ax1.set_ylim([0,180])
    ax2.set_ylim([-8,8])

    xticks = np.linspace(4,14,6)
    ax1.set_xticks(xticks, [r"$%i$"%x for x in xticks])
    yticks = np.linspace(0,180,7)
    ax1.set_yticks(yticks, [r"$%i$"%x for x in yticks])
    yticks = np.linspace(-6,6,7)
    ax2.set_yticks(yticks, [r"$%1.2f$"%x for x in yticks])

    x_m = [sorted(x_s[i])[length//2-1] for i in range(len(x_s))]

    for i in  range(len(x_m)):
        # print(x_m[i])
        ax1.scatter(x_m[i],PS_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax1.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(PS_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))

    for i in  range(len(x_m)):
        ax2.scatter(x_m[i],p3cotPS_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax2.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(p3cotPS_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))
   

    for tmp in [[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
        ax1.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
    # ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
    ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")
    ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")

    ax1.legend(loc='lower right', fontsize=styles.fontsize)

    plt.savefig(op.join(PLTDIR, "phase_shift_plot_close_res.pdf"), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def plot_PS_ERE_res(h5file,show=False):
    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(figs1,figs2))
    plt.subplots_adjust(wspace=0, hspace=0.05)   

    ax1.xaxis.grid()
    ax2.xaxis.grid()
    ax1.yaxis.grid()
    ax2.yaxis.grid()
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, 7.05, -0.867, True)

    slow, shigh = 4,11.1
    ax1.set_xlim([slow,shigh])
    ax1.set_ylim([0,180])
    ax2.set_ylim([-0.75,0.3])

    # x_m   = np.asarray(scat_fit_mean["s_prime"])
    x_s   = np.asarray(scat_fit_spl["s_prime"])
    PS_m   = np.asarray(scat_fit_mean["PS"])
    PS_s   = np.asarray(scat_fit_spl["PS"])
    p3cotPS_m   = np.asarray(scat_fit_mean["p3cotPS_Ecm_prime"])
    p3cotPS_s   = np.asarray(scat_fit_spl["p3cotPS_Ecm_prime"])

    ax2.set_xlabel(r"$s/m_\pi^2$", fontsize=styles.fontsize)
    ax2.set_ylabel(r"$p^3\, \cot(\delta_1)/E_{cm}/m_\pi^2$", fontsize=styles.fontsize)
    ax1.set_ylabel(r"$\delta_1$", fontsize=styles.fontsize)
    
    length = len(x_s[0])
    
    N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(N_Ls,d2s,irreps,lvs))
    
    sarr = np.linspace(slow,shigh,100)
    p2arr = [p2_s(s) for s in sarr]
    fit_model = fm.BW_I_model
    
    fit_param_m = np.asarray([fit_param_mean[fp] for fp in fit_model.param_names])
    yarr_m = np.asarray([fit_model.model(x,*fit_param_m) for x in sarr])

    fit_param_s = np.transpose(np.asarray([fit_param_spl[fp] for fp in fit_model.param_names]))
    yarr_tmp = np.asarray([sorted([fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))]) for x in sarr])

    yarr_e_m_plot = np.asarray([yarr_tmp[i][math.floor(length*(1-num_perc)/2)] for i in range(len(sarr))])
    yarr_e_p_plot = np.asarray([yarr_tmp[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(sarr))])

    ax2.plot(sarr,yarr_m, color = styles.c_10_res)
    ax2.fill_between(sarr, yarr_e_m_plot, yarr_e_p_plot, alpha = 0.5, color = styles.c_10_res)

    PS_m_plot = [delta_res_x(yarr_m[i], p2arr[i])%180 for i in range(len(p2arr))]
    PS_e_m_plot = [delta_res_x(yarr_e_m_plot[i], p2arr[i])%180 for i in range(len(p2arr))]
    PS_e_p_plot = [delta_res_x(yarr_e_p_plot[i], p2arr[i])%180 for i in range(len(p2arr))]

    ax1.plot(sarr,PS_m_plot, color = styles.c_10_res)
    ax1.fill_between(sarr, PS_e_m_plot, PS_e_p_plot, alpha = 0.5, color = styles.c_10_res)

    x_m = [sorted(x_s[i])[length//2-1] for i in range(len(x_s))]

    for i in  range(len(x_m)):
        ax1.scatter(x_m[i],PS_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax1.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(PS_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))

    for i in  range(len(x_m)):
        ax2.scatter(x_m[i],p3cotPS_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax2.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(p3cotPS_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))
    
    for tmp in [[None,0,None,None],[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
        ax1.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
        ax2.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
    # ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
    # ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")
    ax1.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")

    xticks = np.linspace(4,11,8)
    ax1.set_xticks(xticks, [r"$%i$"%x for x in xticks])
    yticks = np.linspace(0,180,7)
    ax1.set_yticks(yticks, [r"$%i$"%x for x in yticks])
    yticks = np.linspace(-0.75,0.25,5)
    ax2.set_yticks(yticks, [r"$%1.2f$"%x for x in yticks])
    
    ax2.legend(loc='upper right', fontsize=styles.fontsize)
    plt.savefig(op.join(PLTDIR, "phase_shift_plot_res.pdf"), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file  = args[2]

    plot_PS_ERE_non_res(h5file, False)
    plot_PS_ERE_close_res(h5file, False)
    plot_PS_ERE_res(h5file, False)