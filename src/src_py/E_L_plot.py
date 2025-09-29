import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
# import fit_scatter
import os.path as op
import os
import sys
import plotting_functions_thesis as pf
import fit_models as fm


def nth(num):
    if num <= 10:
        return 1
    elif num <= 500:
        return num//10
    else:
        return num//100
    
# color_points = "green"
color_fit = "olivedrab"

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
    # else:
    #     plt.xlim([1/40,1/13])
    #     plt.ylim([1,6])


    for tmp in [[None,0,None,None],[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
        ax2.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")
    ax2.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")

    ax2.legend(loc='center right', bbox_to_anchor=(1.2, 0.95))

    ax2.set_xlabel(r"$a/N_L$")
    # plt.ylabel("$E_{CM}$/$m_\\pi$")

    # ax2.text(s="$E/m_\\pi^\\infty$", rotation = "vertical", x=0.0347, y = 1.55, fontsize = 14)
    ax2.text(s="$E_{CM}$/$m_\\pi$", rotation = "vertical", x=0.0347, y = 1.5, fontsize = 14)

    props = dict(facecolor = "white")
    ax1.text(0.02, 0.95, r'$\text{non-resonant}$',
     horizontalalignment='left',
     verticalalignment='top',
     transform = ax1.transAxes,
     fontsize=18,
     bbox=props)

    # print(NL_invs)

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
    y2ticks = np.linspace(2,2.3,4)
    ax1.set_yticks(y2ticks, [r"$%1.1f$"%x for x in y2ticks])

    # plt.xticks([1/14,1/16,1/20,1/24,1/36],["1/14","1/16","1/20","1/24","1/36"])
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

    # plt.subplots_adjust(wspace=0, hspace=0.1)   

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

    ax2.legend(loc='center right', bbox_to_anchor=(1.2, 1))

    ax2.set_xlabel(r"$a/N_L$")

    ax1.set_ylabel("$E_{CM}$/$m_\\pi$")
    ax2.set_ylabel("$E_{CM}$/$m_\\pi$")
    # ax2.text(s="$E/m_\\pi^\\infty$", rotation = "vertical", x=0.0347, y = 1.55, fontsize = 14)

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
    # for x in NLs:
    #     if x not in NL_label:
    #         NL_label.append(x)
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
    # ax1.set_xticks(NL_inv_label, ["" for x in NL_label])
    yticks = np.linspace(1.5,3.5,5)
    ax1.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    yticks = np.linspace(1,5.5,10)
    ax2.set_yticks(yticks, [r"$%1.1f$"%x for x in yticks])
    
    plt.savefig(op.join(PLTDIR, "E_CM_L_705_multi_levels_%r.pdf"%(levels)), bbox_inches='tight')
    if show:
        plt.show()
    plt.clf()

# def plot_any(h5file,beta,m0,xaxis="p2star_prime",yaxis="p3cotPS_prime",fit_model=None,outname=None,show=False):
#     fit = fit_model != None
#     plt.rcParams['figure.figsize'] = [10, 6]
#     fontsize = 14
#     font = {'size'   : fontsize}
#     matplotlib.rc('font', **font)
#     fig, ax = plt.subplots()
#     plt.grid()
#     info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, beta, m0, fit)
#     xlim = xlim_f(m0,xaxis)
#     if xlim == None:
#         ax.set_xlim(auto=True)
#     else: 
#         ax.set_xlim(xlim)
#     ylim = ylim_f(m0,yaxis)
#     if ylim == None:
#         ax.set_ylim(auto=True)
#     else: 
#         ax.set_ylim(ylim)

#     x_m   = np.asarray(scat_fit_mean[xaxis])
#     x_s   = np.asarray(scat_fit_spl[xaxis])
#     y_m   = np.asarray(scat_fit_mean[yaxis])
#     y_s   = np.asarray(scat_fit_spl[yaxis])

#     xlabel = xlabel_f(xaxis)
#     plt.xlabel(xlabel)
#     ylabel = ylabel_f(yaxis)
#     plt.ylabel(ylabel)
    
#     length = len(x_s[0])
    
#     N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
#     dvecs = scat_fit_mean["dvec"]
#     dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
#     d2s = [np.dot(d,d) for d in dvecs]
#     lvs = info["lv"]
#     irreps = info["irrep"]
#     plot_args = list(zip(N_Ls,d2s,irreps,lvs))
    
#     for i in  range(len(x_m)):
#         ax.scatter(x_m[i],y_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))#, s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
#         sorted_indices = np.argsort(x_s[i])
#         ax.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))

#     if fit:
#         x_nf_m       = np.asarray(scat_nf_mean[xaxis])
#         x_nf_s   = np.asarray(scat_nf_spl[xaxis])
#         y_nf_m       = np.asarray(scat_nf_mean[yaxis])
#         y_nf_s   = np.asarray(scat_nf_spl[yaxis])
        
#         N_Ls_nf = [int(x) for x in scat_nf_mean["N_L"]]
#         dvecs_nf = scat_nf_mean["dvec"]
#         dvecs_nf = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs_nf]
#         d2s_nf = [np.dot(d,d) for d in dvecs_nf]
#         lvs_nf = info_nf["lv"]
#         irreps_nf = info_nf["irrep"]
#         plot_args_nf = list(zip(N_Ls_nf,d2s_nf,irreps_nf,lvs_nf))
        
#         xarr = np.linspace(xlim[0], xlim[1], 600)
#         for i in  range(len(x_nf_m)):
#             ax.scatter(x_nf_m[i],y_nf_m[i], color = "grey", ls = pf.ls(*plot_args_nf[i]), marker = pf.marker(*plot_args_nf[i]))#, s = 10*pf.ms(*plot_args_nf[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
#             sorted_indices = np.argsort(x_nf_s[i])
#             ax.plot(x_nf_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_nf_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = "grey", ls = pf.ls(*plot_args_nf[i]))
#         plt.plot([-1,-1],[-1,-1], color = "grey", label = "not fitted")
        
#         fit_param_m = np.asarray([fit_param_mean[fp] for fp in fit_model.param_names])
#         yarr_m = np.asarray([fit_model.model(x,*fit_param_m) for x in xarr])

#         fit_param_s = np.transpose(np.asarray([fit_param_spl[fp] for fp in fit_model.param_names]))
#         yarr_tmp = np.asarray([sorted([fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))]) for x in xarr])

#         if fit_model.yaxis == yaxis:
#             yarr_s = yarr_tmp
#             yarr_m_plot = yarr_m
#         else:
#             yarr_s = np.asarray([sorted([from_to(fit_model.yaxis,yaxis)(xarr[i],yarr_tmp[i,j]) for j in range(len(yarr_tmp[0]))]) for i in range(len(yarr_tmp))])
#             yarr_m_plot = np.vectorize(from_to(fit_model.yaxis,yaxis))(xarr,yarr_m)
            
#         yarr_med_plot = np.asarray([yarr_s[i][length//2-1] for i in range(len(xarr))])
#         yarr_e_m_plot = np.asarray([yarr_s[i][math.floor(length*(1-num_perc)/2)] for i in range(len(xarr))])
#         yarr_e_p_plot = np.asarray([yarr_s[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(xarr))])

#         plt.plot(xarr,yarr_m_plot, color = color_fit)
#         # plt.plot(xarr,yarr_med_plot, color = "blue")
#         plt.fill_between(xarr, yarr_e_m_plot, yarr_e_p_plot, alpha = 0.3, color = color_fit)

#         if xaxis == "s_prime" and yaxis == "sigma_prime":
#             yarr_14 = np.asarray([sigma_14(x,0.56,4.4) for x in xarr])
#             plt.plot(xarr,yarr_14, color = "green", label = "14-dim")

#         if fit_model.yaxis == yaxis:
#             y_e = np.asarray([sorted(y_s[i])[length//2] for i in range(len(y_s))])
#             y_pred = np.asarray([fit_model.model(x, *fit_param_m) for x in x_m])
#             print("beta = %1.3f, m0 = %1.3f"%(beta, m0))
#             print(fit_model.name)
#             for i in range(fit_model.num_params):
#                 param_med = (sorted(np.transpose(fit_param_s)[i])[length//2])
#                 param_e_m = (sorted(np.transpose(fit_param_s)[i])[math.floor(length*(1-num_perc)/2)])
#                 param_e_p = (sorted(np.transpose(fit_param_s)[i])[math.ceil(length*(1+num_perc)/2)])
#                 print("%s  =  %.3f^{+%.3f}{-%.3f}"%(fit_model.param_names[i], fit_param_m[i], param_med-param_e_m, param_e_p-param_med))
#             chi2 = np.sum(((y_m - y_pred) / y_e) ** 2)
#             ndof = len(y_m) - fit_model.num_params  # degrees of freedom
#             chi2_ndof = chi2 / ndof
#             print("chi^2 = %1.3f,  dof = %i,  chi2/dof = %f"%(chi2,ndof,chi2_ndof),end="\n\n")
#             # ax.text(0.05, 0.85, "chi^2/dof=%1.3f"%chi2_ndof, transform=ax.transAxes, fontsize=12, verticalalignment="top", bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
#     # ax.text(0.05, 0.95, "$\\beta$=%1.3f, $m_0$=%1.3f"%(beta,m0), transform=ax.transAxes, fontsize=12, verticalalignment="top", bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))

#     for tmp in [[None,0,None,None],[None,1,None,None],[None,2,None,None],[None,3,None,None]]:
#         ax.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = r"$|p|=%i$"%(tmp[1]))
#     ax.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",0), label = r"$E^{A_1}_0$")
#     ax.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"A1",1), label = r"$E^{A_1}_1$")
#     ax.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(None,None,"B1",0), label = r"$E^{\rho}$")
#     # for tmp in [[None,0,"T1",0],[None,1,"E",0],[None,2,"B1",0],[None,3,"E",0],[None,1,"A1",0],[None,1,"A1",1],[None,2,"A1",0],[None,2,"A1",1],[None,3,"A1",0],[None,3,"A1",1]]:
#     #     plt.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = "p=%i, %s, lv=%i"%(tmp[1],tmp[2],tmp[3]))
#     # for tmp in [[14,None,None,None],[16,None,None,None],[20,None,None,None],[24,None,None,None],[36,None,None,None]]:
#     #     plt.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(*tmp), label = "$N_L$=%i"%(tmp[0]))
#     ax.legend(loc='center right', bbox_to_anchor=(1.35, 0.5))
#     fit_str = "" if fit_model == None else "_fit_%s"%fit_model.name
#     out_str = "b%1.3f_m0%1.3f"%(beta,m0) if outname == None else outname
#     plt.savefig(op.join(PLTDIR, "%s_%s%s__%s.pdf"%(yaxis,xaxis,fit_str,out_str)), bbox_inches='tight')
#     if show:
#         plt.show()
#     plt.close(fig)

if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file  = args[2]
    # fit = args[3] == "True"

    os.makedirs(PLTDIR, exist_ok=True)

    plot_E_CM_L_multi_non_res(h5file, 6.9, -0.92, show=False)
    plot_E_CM_L_multi_705(h5file, show=False)