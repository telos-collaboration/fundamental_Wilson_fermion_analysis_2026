import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
# import fit_scatter
import os.path as op
import os
import sys
import plotting_functions as pf
import fit_models as fm

import styles

def nth(num):
    if num <= 10:
        return 1
    elif num <= 500:
        return num//10
    else:
        return num//100
    
color_fit = "olivedrab"

num_perc = math.erf(1/np.sqrt(2))

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

def E_pipi(mpi,p12,p22,L):
    return np.sqrt(mpi**2+(2*np.pi/L)**2*p12)+np.sqrt(mpi**2+(2*np.pi/L)**2*p22)  

def color(d2):
    colors = ["cyan","orange", "green", "blueviolet"]
    return colors[d2]

def plot_E_CM_L_multi_non_res(h5file,beta,m0,levels=False,outname=None,show=False):
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, beta, m0, False)

    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=False,figsize=(10,7))
    plt.subplots_adjust(wspace=0, hspace=0.05)   

    NLs = [int(x) for x in info["NL"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(NLs,d2s,irreps,lvs))
    mpi = info["mpi"][0]
    mrho = info["mrho"][0]

    ECMs = np.asarray(scat_fit_mean["E_cm_prime"])
    ECM_spl = np.asarray(scat_fit_spl["E_cm_prime"])
    length = len(ECM_spl[0])
    ECM_errms = [abs(ECMs[i]-sorted(ECM_spl[i])[math.floor(length*(1-num_perc)/2)]) for i in range(len(ECMs))]
    ECM_errps = [abs(ECMs[i]-sorted(ECM_spl[i])[math.ceil(length*(1+num_perc)/2)]) for i in range(len(ECMs))]
    NL_invs = [1/x for x in NLs]
    
    for i in range(len(ECMs)):
        ECM_errms[i] = 0 if ECM_errms[i] > 1 else ECM_errms[i]
        ECM_errps[i] = 0 if ECM_errps[i] > 1 else ECM_errps[i]
        if irreps[i] == "T1":
            print(beta, m0, NLs[i], ECMs[i],mpi)
        if lvs[i] == 1:
            ax1.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
        else:
            if irreps[i] == "A1":
                ax2.errorbar([NL_invs[i]+0.0006,],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
            else:
                ax2.errorbar([NL_invs[i]-0.0006,],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   

    # plt.axhline(1,c="black", ls="dotted", label = r"$m_\pi$")
    # plt.axhline(mrho/mpi,c="red", ls="dotted", label = r"$m_\rho$")
    ax1.axhline(2,c="black",label = r"2$m_\pi$")
    ax2.axhline(2,c="black",label = r"2$m_\pi$")
    # ax2.axhline(4,c="black",label = r"4$m_\pi$")
    ax1.yaxis.grid()
    ax2.yaxis.grid()
    # ax1.set_title(r"$\\beta$ = %1.1f, $am_0$ = %1.2f$"%(beta,m0))
    # if levels:
    #     plt.plot([0,0],[0,0],c="grey", label = "non-int")
    ax2_y_top = 1.52
    ax1_y_bot = 1.97
    ax1.set_ylim([ax1_y_bot,2.29])
    ax1.set_xlim([1/26,1/13])
    ax2.set_xlim([1/26,1/13])
    ax2.set_ylim([1.3,ax2_y_top])


    for tmp in [[None,0,None,None],[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
        ax2.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")

    ax2.legend(loc='center right', bbox_to_anchor=(1.25, 0.95))

    ax2.set_xlabel(r"$a/L$")
    # plt.ylabel("$E_{CM}$/$m_\\pi$")

    ax2.text(s="$E_{CM}$/$m_\\pi$", rotation = "vertical", x=0.0347, y = 1.5, fontsize = 14)

    props = dict(facecolor = "white")
    ax1.text(0.02, 0.95, r'$\text{non-resonant}$',
     horizontalalignment='left',
     verticalalignment='top',
     transform = ax1.transAxes,
     fontsize=18,
     bbox=props)

    NL_label = []
    for x in NLs:
        if x not in NL_label:
            NL_label.append(x)
    NL_inv_label = [1/x for x in NL_label]

    for x in NL_inv_label:
        ax1.axvline(x+0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
        ax1.axvline(x-0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
        ax2.axvline(x+0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
        ax2.axvline(x-0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")

    ax1.set_xticks([])
    ax2.set_xticks(NL_inv_label, [r"$1/%i$"%x for x in NL_label])
    y1ticks = np.linspace(1.3,1.5,5)
    ax2.set_yticks(y1ticks, [r"$%1.2f$"%x for x in y1ticks])
    y2ticks = np.linspace(2,2.3,5)
    ax1.set_yticks(y2ticks, [r"$%1.1f$"%x for x in y2ticks])

    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
    ax1.axhline(ax1_y_bot, ls = "dashed", color = "black")
    ax2.axhline(ax2_y_top, ls = "dashed", color = "black")

    ax1.spines[['bottom',]].set_visible(False)
    ax2.spines[['top',]].set_visible(False)
    
    plt.savefig(op.join(PLTDIR, "E_CM_L_b_%1.3f_m0_%1.3f_multi_levels_%r.pdf"%(beta,m0,levels)), bbox_inches='tight')
    if show:
        plt.show()
    plt.clf()

def plot_E_CM_L_multi_705(h5file,levels=False,outname=None,show=False):
    fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(10,10))
    plt.subplots_adjust(wspace=0, hspace=0.02)   
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, 7.05, 0.863, False)

    NLs = [int(x) for x in info["NL"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(NLs,d2s,irreps,lvs))
    mpi = info["mpi"][0]
    mrho = info["mrho"][0]

    ECMs = np.asarray(scat_fit_mean["E_cm_prime"])
    ECM_spl = np.asarray(scat_fit_spl["E_cm_prime"])
    length = len(ECM_spl[0])
    ECM_errms = [abs(ECMs[i]-sorted(ECM_spl[i])[math.floor(length*(1-num_perc)/2)]) for i in range(len(ECMs))]
    ECM_errps = [abs(ECMs[i]-sorted(ECM_spl[i])[math.ceil(length*(1+num_perc)/2)]) for i in range(len(ECMs))]
    NL_invs = [1/x for x in NLs]
    
    for i in range(len(ECMs)):
        ECM_errms[i] = 0 if ECM_errms[i] > 1 else ECM_errms[i]
        ECM_errps[i] = 0 if ECM_errps[i] > 1 else ECM_errps[i]
        if lvs[i] == 1:
            ax1.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
        else:
            if irreps[i] == "A1":
                ax1.errorbar([NL_invs[i]+0.0006,],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
            else:
                ax1.errorbar([NL_invs[i]-0.0006,],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   


    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, 7.05, 0.867, False)

    plt.subplots_adjust(wspace=0, hspace=0.1)   

    NLs = [int(x) for x in info["NL"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(NLs,d2s,irreps,lvs))
    mpi = info["mpi"][0]
    mrho = info["mrho"][0]

    ECMs = np.asarray(scat_fit_mean["E_cm_prime"])
    ECM_spl = np.asarray(scat_fit_spl["E_cm_prime"])
    length = len(ECM_spl[0])
    ECM_errms = [abs(ECMs[i]-sorted(ECM_spl[i])[math.floor(length*(1-num_perc)/2)]) for i in range(len(ECMs))]
    ECM_errps = [abs(ECMs[i]-sorted(ECM_spl[i])[math.ceil(length*(1+num_perc)/2)]) for i in range(len(ECMs))]
    NL_invs = [1/x for x in NLs]
    
    for i in range(len(ECMs)):
        ECM_errms[i] = 0 if ECM_errms[i] > 1 else ECM_errms[i]
        ECM_errps[i] = 0 if ECM_errps[i] > 1 else ECM_errps[i]
        if lvs[i] == 1:
            ax2.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
        else:
            if irreps[i] == "A1":
                ax2.errorbar([NL_invs[i]+0.0006,],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
            else:
                ax2.errorbar([NL_invs[i]-0.0006,],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   

    ax2.axhline(2,c="black",label = r"2$m_\pi$")
    ax2.yaxis.grid()
    ax2.axhline(4,c="black",label = r"4$m_\pi$", ls= "dashed")
    ax2.set_ylim([1,6])
    ax1.axhline(2,c="black",label = r"2$m_\pi$")
    ax1.yaxis.grid()
    ax1.axhline(4,c="black",label = r"4$m_\pi$", ls= "dashed")
    ax1.set_ylim([1,6])


    for tmp in [[None,0,None,None],[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
        ax2.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")

    ax2.legend(loc='center right', bbox_to_anchor=(1.25, 1))

    ax2.set_xlabel(r"$a/L$")

    ax1.set_ylabel("$E_{CM}$/$m_\\pi$")
    ax2.set_ylabel("$E_{CM}$/$m_\\pi$")

    props = dict(facecolor = "white")
    ax1.text(0.02, 0.96, r'$\text{close to resonant}$',
     horizontalalignment='left',
     verticalalignment='top',
     transform = ax1.transAxes,
     fontsize=18,
     bbox=props)
    ax2.text(0.02, 0.96, r'$\text{resonant}$',
     horizontalalignment='left',
     verticalalignment='top',
     transform = ax2.transAxes,
     fontsize=18,
     bbox=props)

    NL_label = [16,20,24,36]
    NL_inv_label = [1/x for x in NL_label]

    for x in NL_inv_label:
        ax1.axvline(x+0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
        ax1.axvline(x-0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
        ax2.axvline(x+0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
        ax2.axvline(x-0.0012, alpha = 0.5, color = "grey", lw=1, ls = "dotted")
    
    ax1.set_ylim([1.3,3.9])
    ax2.set_xlim([1/40,1/15])
    ax2.set_ylim([0.8,6])

    print(NL_label)

    ax2.set_xticks(NL_inv_label, [r"$1/%i$"%x for x in NL_label])
    yticks = np.linspace(1.5,3.5,5)
    ax1.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    yticks = np.linspace(1,5.5,10)
    ax2.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    
    plt.savefig(op.join(PLTDIR, "E_CM_L_705_multi_levels_%r.pdf"%(levels)), bbox_inches='tight')
    if show:
        plt.show()
    plt.clf()

if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file  = args[2]

    os.makedirs(PLTDIR, exist_ok=True)

    plot_E_CM_L_multi_non_res(h5file, 6.9, -0.92, show=False)
    plot_E_CM_L_multi_705(h5file, show=False)