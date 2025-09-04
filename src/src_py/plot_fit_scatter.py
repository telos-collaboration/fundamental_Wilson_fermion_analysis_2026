import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
import fit_scatter
import os.path as op
import os
import sys
import plotting_functions as pf

num_perc = math.erf(1/np.sqrt(2))

########################################## Plot p3cotPS ##########################################

    
# def marker_NL(NL):                     # maybe to be replaced by input file
#     if NL == 14:
#         return "D"
#     elif NL == 16:
#         return "p"
#     elif NL == 20:
#         return "X"
#     elif NL == 24:
#         return "o"
#     elif NL == 36:
#         return "*"
#     else:
#         raise RuntimeError("Wrong NL in marker: NL=%i"%(NL))
    
# def ms_p(p):                     # maybe to be replaced by input file
#     if p == 0:
#         return 3
#     elif p == 1:
#         return 4
#     elif p == 2:
#         return 5
#     elif p == 3:
#         return 6
#     else:
#         raise RuntimeError("Wrong p in ms_p: %i, %i"%(p))

# def color_irrep_lv(irrep,lv,p):
#     if irrep == "A1":
#         if p == 1:
#             if lv == 0:
#                 return "red"
#             elif lv == 1:
#                 return "darkred"
#         elif p == 2:
#             if lv == 0:
#                 return "yellow"
#             elif lv == 1:
#                 return "gold"
#         elif p == 3:
#             if lv == 0:
#                 return "fuchsia"
#             elif lv == 1:
#                 return "purple"
#     elif irrep == "E":
#         if p == 1:
#             if lv == 0:
#                 return "blue"
#             elif lv == 1:
#                 return "darkblue"
#         elif p == 3:
#             if lv == 0:
#                 return "lightseagreen"
#             elif lv == 1:
#                 return "mediumturquise"
#     elif irrep == "B1":
#         if lv == 0:
#             return "green"
#         elif lv == 1:
#             return "darkgreen"
#     elif irrep == "T1":
#         if lv == 0:
#             return "peru"
#     raise ValueError("wrong irrep or lv in color_irrep_lv(): %i, %i"%(irrep,lv))

# def ls_NL(NL):
#     if NL == 14:
#         return "solid"
#     elif NL == 16:
#         return (0,(1,1))
#     elif NL == 20:
#         return "dashed"
#     elif NL == 24:
#         return "dashdot"
#     elif NL == 36:
#         return "dotted"
#     else:
#         raise ValueError("Wrong NL given to ls_NL()")

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
    scat_nf_mean = {}
    scat_nf_spl = {}
    info_nf = {}
    fit_beta_m = "fit_b%f_m%f"%(beta,m0)
    with h5py.File(h5file_scatter_fit,"r") as hfile:
        for ens in hfile:
            if str(beta) in ens and str(m0) in ens:
                for key in hfile[fit_beta_m]["mean"]:
                    fit_param_mean[key] = hfile[fit_beta_m]["mean"][key][()]
                    fit_param_spl[key] = hfile[fit_beta_m]["sample"][key][()]
                for P in hfile[ens]:
                    if P[0] == "p":
                        for irrep in hfile[ens][P]:
                            for lv in hfile[ens][P][irrep]:
                                if lv[:2] == "lv":
                                    # print(ens+P+irrep+lv)
                                    if hfile[ens][P][irrep][lv]["fit"][()]:
                                        info.setdefault("lv",[]).append(int(lv[2:]))
                                        info.setdefault("irrep",[]).append(irrep)
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
                                        scat_nf_mean.setdefault("irrep",[]).append(irrep)
                                        for key in hfile[ens][P][irrep][lv]["mean"]:
                                            if key == "dvec" or key == "N_L":
                                                scat_nf_mean.setdefault(key,[]).append(hfile[ens][P][irrep][lv]["mean"][key][()])
                                            else:
                                                scat_nf_mean.setdefault(key,[]).append(float(np.real(hfile[ens][P][irrep][lv]["mean"][key][()])))
                                                scat_nf_spl.setdefault(key,[]).append([float(x) for x in np.real(hfile[ens][P][irrep][lv]["sample"][key][()])])
    return info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl

def plot_p3cotPS(h5file_scatter_fit,beta,m0,fit=False,outname=None,show=False):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    plt.grid()
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data_p3cotPS(h5file_scatter_fit, beta, m0)
    xlim = [0,0.3] if beta==6.9 else [0,3]
    ax.set_xlim(xlim)
    ylim = [-1,1] if beta==6.9 else [-4,4]
    ax.set_ylim(ylim)

    plt.xlabel(r"$p^{\star^2}/m_\pi^2$")
    x_plot       = np.asarray(scat_fit_mean["p2star_prime"])
    x_plot_sam   = np.asarray(scat_fit_spl["p2star_prime"])
    y_plot       = np.asarray(scat_fit_mean["p3cotPS_prime"])
    y_plot_sam   = np.asarray(scat_fit_spl["p3cotPS_prime"])
    x_n_plot     = np.asarray(scat_nf_mean["p2star_prime"])                                        # n marks that it was not fitted
    x_n_plot_sam = np.asarray(scat_nf_spl["p2star_prime"])
    y_n_plot     = np.asarray(scat_nf_mean["p3cotPS_prime"])
    y_n_plot_sam = np.asarray(scat_nf_spl["p3cotPS_prime"])
    plt.ylabel(r"$p^3\, \cot(\delta)/m_\pi^3$")

    print(y_n_plot_sam.shape)
    
    length = len(x_plot_sam[0])
    
    N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(N_Ls,d2s,irreps,lvs))
    
    N_Ls_nf = [int(x) for x in scat_nf_mean["N_L"]]
    dvecs_nf = scat_nf_mean["dvec"]
    dvecs_nf = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs_nf]
    d2s_nf = [np.dot(d,d) for d in dvecs_nf]
    lvs_nf = info_nf["lv"]
    irreps_nf = info_nf["irrep"]
    plot_args_nf = list(zip(N_Ls_nf,d2s_nf,irreps_nf,lvs_nf))


    for i in  range(len(x_plot)):
        if 0<x_plot[i]<3: 
            ax.scatter(x_plot[i],y_plot[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]), s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
            sorted_indices = np.argsort(x_plot_sam[i])
            ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_plot_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))
        else:
            raise ValueError("Fitted momentum is not in elastic threshhold in plotting.py!!!")

    for i in  range(len(x_n_plot)):
        ax.scatter(x_n_plot[i],y_n_plot[i], color = "grey", ls = pf.ls(*plot_args_nf[i]), marker = pf.marker(*plot_args_nf[i]), s = 10*pf.ms(*plot_args_nf[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_n_plot_sam[i])
        ax.plot(x_n_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_n_plot_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = "grey", ls = pf.ls(*plot_args_nf[i]))


    xarr = np.linspace(xlim[0], xlim[1])
    
    if fit:
        a1_1 = fit_param_mean["a1_1"]
        r1_1 = fit_param_mean["r1_1"]

        a1_1_smp = fit_param_spl["a1_1"]
        r1_1_smp = fit_param_spl["r1_1"]

        yarr = [fit_scatter.ERE_1(x,a1_1,r1_1) for x in xarr]
        yarr_smp = [sorted([fit_scatter.ERE_1(x,a1_1_smp[i],r1_1_smp[i]) for i in range(len(a1_1_smp))]) for x in xarr]

        yarr_m = [yarr_smp[i][math.floor(length*(1-num_perc)/2)] for i in range(len(yarr_smp))]
        yarr_p = [yarr_smp[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(yarr_smp))]

        plt.plot(xarr,yarr, color = "blue")
        plt.fill_between(xarr, yarr_m, yarr_p, alpha = 0.3, color = "blue")

    plt.plot([-1,-1],[-1,-1], color = "grey", label = "not fitted")
    for tmp in [[None,0,"T1",0],[None,1,"E",0],[None,2,"B1",0],[None,3,"E",0],[None,1,"A1",0],[None,1,"A1",1],[None,2,"A1",0],[None,2,"A1",1],[None,3,"A1",0],[None,3,"A1",1]]:
        plt.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = "p=%i, %s, lv=%i"%(tmp[1],tmp[2],tmp[3]))
    for tmp in [[14,None,None,None],[16,None,None,None],[20,None,None,None],[24,None,None,None],[36,None,None,None]]:
        plt.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(*tmp), label = "$N_L$=%i"%(tmp[0]))
    ax.legend(loc='center right', bbox_to_anchor=(1.35, 0.5))
    if outname == None:    
        plt.savefig(op.join(PLTDIR, "p3cotPS_b%f_m0%f_fit_%r.pdf"%(beta,m0,fit)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR, "p3cotPS_"+outname+"_fit_%r.pdf"%fit), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def p2_s_prime(s):
    return s/4-1

# def plot_p3cotPS_ECM(h5file_scatter_fit,beta,m0,fit=False,outname=None,show=False):
#     plt.rcParams['figure.figsize'] = [10, 6]
#     fontsize = 14
#     font = {'size'   : fontsize}
#     matplotlib.rc('font', **font)
#     fig, ax = plt.subplots()
#     plt.grid()
#     fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data_p3cotPS(h5file_scatter_fit, beta, m0)
#     xlim = [4,16]
#     ax.set_xlim(xlim)
#     ylim = [-2,2]
#     ax.set_ylim(ylim)

#     plt.xlabel(r"$s/m_\pi^2$")
#     x_plot       = np.asarray(scat_fit_mean["s_prime"])
#     x_plot_sam   = np.asarray(scat_fit_spl["s_prime"])
#     y_plot       = np.asarray(scat_fit_mean["p3cotPS_Ecm_prime"])
#     y_plot_sam   = np.asarray(scat_fit_spl["p3cotPS_Ecm_prime"])
#     x_n_plot     = np.asarray(scat_nf_mean["s_prime"])                                        # n marks that it was not fitted
#     x_n_plot_sam = np.asarray(scat_nf_spl["s_prime"])
#     y_n_plot     = np.asarray(scat_nf_mean["p3cotPS_Ecm_prime"])
#     y_n_plot_sam = np.asarray(scat_nf_spl["p3cotPS_Ecm_prime"])
#     plt.ylabel(r"$p^3\, \cot(\delta)/(E_{cm}m_\pi^2$")
    
#     length = len(x_plot_sam[0])
#     num_perc = math.erf(1/np.sqrt(2))
    
#     N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
#     dvecs = scat_fit_mean["dvec"]
#     dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
#     d2s = [np.dot(d,d) for d in dvecs]
    
#     dvec_ns = scat_nf_mean["dvec"]
#     dvec_ns = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvec_ns]
#     for i in  range(len(x_plot)):
#         if 4<x_plot[i]<16: 
#             ax.scatter(x_plot[i],y_plot[i], label = "%s |P|=%i, NL=%i"%(scat_fit_mean["irrep"][i],d2s[i],N_Ls[i]), color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]), marker = marker(scat_fit_mean["irrep"][i]),s=60)
#             sorted_indices = [int(x) for x in np.argsort(x_plot_sam[i])]
#             ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_plot_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))
#         else:
#             raise ValueError("Fitted momentum is not in elastic threshhold in plotting.py!!!")

#     for i in  range(len(x_n_plot)):
#         ax.scatter(x_n_plot[i],y_n_plot[i], color = "grey", ls = ls_P(dvec_ns[i]), marker = marker(scat_nf_mean["irrep"][i]),s=60)
#         sorted_indices = np.argsort(x_n_plot_sam[i])
#         ax.plot(x_n_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_n_plot_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = "grey", ls = ls_P(dvec_ns[i]))

#     xarr = np.linspace(xlim[0], xlim[1])
    
#     if fit:
#         p2arr = [p2_s_prime(s) for s in xarr]
#         m_R = fit_param_mean["m_R_D"]
#         gVPP2 = fit_param_mean["gVPP2_D"]

#         m_R_smp = fit_param_spl["m_R_D"]
#         gVPP2_smp = fit_param_spl["gVPP2_D"]

#         yarr = [fit_scatter.RES_Drach(x,m_R,gVPP2) for x in p2arr]
#         yarr_smp = [sorted([fit_scatter.RES_Drach(x,m_R_smp[i],gVPP2_smp[i]) for i in range(len(m_R_smp))]) for x in p2arr]

#         yarr_m = [yarr_smp[i][math.floor(length*(1-num_perc)/2)] for i in range(len(yarr_smp))]
#         yarr_p = [yarr_smp[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(yarr_smp))]

#         plt.plot(xarr,yarr, color = "blue")
#         plt.fill_between(xarr, yarr_m, yarr_p, alpha = 0.3, color = "blue")

#     if m0 == -0.867:
#         plt.axvline(5.756, label="naive rho")

#     ax.legend(loc='center right', bbox_to_anchor=(1.35, 0.5))
#     if outname == None:    
#         plt.savefig(op.join(PLTDIR, "p3cotPS_Ecm_b%f_m0%f_fit_%r.pdf"%(beta,m0,fit)), bbox_inches='tight')
#     else:    
#         plt.savefig(op.join(PLTDIR, "p3cotPS_Ecm_"+outname+"_fit_%r.pdf"%fit), bbox_inches='tight')
#     if show:
#         plt.show()
#     plt.close(fig)

########################################## Plot sigma1 ##########################################

# def sigma_ERE_s_wave(s, a, r):
#     p2 = p2_s_prime(s)
#     if p2 == 0:
#         return 0
#     cot_PS = (-1/a+p2*r/2)/np.sqrt(p2)
#     return 4*np.pi/(cot_PS**2+1)/p2
    
# def sigma_of_P3cotPS(P3cotPS, p2):
#     if p2 == 0:
#         return 0
#     else:
#         cot_PS = P3cotPS/(p2**(3/2))
#         return 4*np.pi*3/(cot_PS**2+1)/p2

# def plot_sigma_1(h5file_scatter_fit,beta,m0,fit=False,outname=None,show=False):         # HAS TO BE FIXED WITH NEW data format!!!
    
#     plt.rcParams['figure.figsize'] = [10, 6]
#     fontsize = 14
#     font = {'size'   : fontsize}
#     matplotlib.rc('font', **font)
#     fig, ax = plt.subplots()
#     plt.grid()
#     res, res_smp = get_data_p3cotPS(h5file_scatter_fit, beta, m0)

#     xlim = [4,6.5]

#     plt.xlabel(r"$s/m_\pi^2$")
#     x_plot_sam = np.transpose(res_smp["s_prime"])

#     plt.ylabel(r"$\sigma m_\pi^2$")        
    
#     length = len(x_plot_sam[0])
#     num_perc = math.erf(1/np.sqrt(2))

#     sarr = np.linspace(xlim[0],xlim[1],5000)
#     p2arr = [x/4-1 for x in sarr]

#     yarr_14 = [sigma_ERE_s_wave(s, 0.52, 6.7) for s in sarr]

#     plt.plot(sarr, yarr_14, color = "red", label = "2405.06506")
#     a1_1 = res["a1_1"]
#     r1_1 = res["r1_1"]

#     a1_1_smp = res_smp["a1_1"]
#     r1_1_smp = res_smp["r1_1"]

#     yarr = [sigma_of_P3cotPS(fit_scatter.ERE_1(x,a1_1,r1_1), x) for x in p2arr]
#     plt.plot(sarr,yarr,label="This work")

#     yarr_smp = []
#     for p2 in p2arr:
#         p3cotPS_smp = [fit_scatter.ERE_1(p2,a1_1_smp[i],r1_1_smp[i]) for i in range(len(a1_1_smp))]
#         sorted_indices = np.argsort(p3cotPS_smp)[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)]
#         yarr_smp.append(np.asarray([sigma_of_P3cotPS(fit_scatter.ERE_1(p2,a1_1_smp[i],r1_1_smp[i]), p2) for i in range(len(a1_1_smp))])[sorted_indices])


#     yarr_m = [min(yarr_smp[i]) for i in range(len(yarr_smp))]
#     yarr_p = [max(yarr_smp[i]) for i in range(len(yarr_smp))]

#     plt.fill_between(sarr, yarr_m, yarr_p, alpha = 0.3)

#     if outname == None:    
#         plt.savefig(op.join(PLTDIR,"scattering/sigma1_b%f_m0%f_fit_%r.pdf"%(beta,m0,fit)), bbox_inches='tight')
#     else:    
#         plt.savefig(op.join(PLTDIR,"sigma1_"+outname+"_fit_%r.pdf"%fit), bbox_inches='tight')
#     if show:
#         plt.show()
#     plt.close(fig)

if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file_scatter_fit  = args[2]

    os.makedirs(PLTDIR, exist_ok=True)
    
    plot_p3cotPS(h5file_scatter_fit,6.9,-0.92,True,outname="non_res",show=False)
    plot_p3cotPS(h5file_scatter_fit,7.05,-0.863,True,outname="close_res",show=False)
    plot_p3cotPS(h5file_scatter_fit,7.05,-0.867,True,outname="res",show=False)