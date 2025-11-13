import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import math
import os.path as op
import os
import sys
import plotting_functions as pf
import fit_models as fm

import styles

mpl.rcParams['lines.markersize'] = 9

figs1, figs2 = 10,7

def nth(num):
    if num <= 10:
        return 1
    elif num <= 500:
        return num//10
    else:
        return num//100

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

def read_from_hdf(filename):
    res, res_tmp = [{},{}]
    with h5py.File("data_assets/14_dim/"+filename+".hdf5","r") as hfile:
        for key in hfile.keys():
            if key[:4] == "orig":
                res[key[5:]] = hfile[key][()]
            if key[:4] == "samp":
                res_tmp[key[7:]] = hfile[key][()]
    return res, res_tmp

def s_p2(p2):
    return 4*(1+p2)

def p2_s(s):
    return s/4-1

def delta_x(x,p2):
    return np.arctan(p2**(3/2)/x)*360/(2*np.pi)

def delta_res_x(x,p2):
    return np.arctan(p2**(3/2)/(2*np.sqrt(1+p2)*x))*360/(2*np.pi)

def sigma_of_p3cotPS_Ecm_prime(s, p3cotPS_Ecm_prime):
    cot_PS = p3cotPS_Ecm_prime*np.sqrt(s)/(s/4-1)**(3/2)
    return 12*np.pi/((s/4-1)*(1+cot_PS**2))

def sigma_of_p3cotPS_prime(s, p3cotPS_prime):
    cot_PS = p3cotPS_prime/(s/4-1)**(3/2)
    return 12*np.pi/((s/4-1)*(1+cot_PS**2))

def UTE(P2, a, b):
    """
    Second order univesal threshold expansion
    """
    return a+P2*b

def sigma_of_pcotPS_prime(s, pcotPS_prime):
    # cot_PS = pcotPS_prime/(s/4-1)**(1/2)
    return 4*np.pi/(((s/4-1)+pcotPS_prime**2))

def plot_sigma(h5file,show=False,units=False):
    fig, ax = plt.subplots(figsize=(figs1,figs2))
    plt.subplots_adjust(wspace=0, hspace=0.05)   

    plt.grid()
    slow, shigh = 4,11.1
    ax.set_xlim([slow,shigh])
    ax.set_ylim([0,100])

    info, info_nf, fit_param_mean_NR, fit_param_spl_NR, scat_fit_mean_NR, scat_fit_spl_NR, scat_nf_mean_NR, scat_nf_spl_NR = get_data(h5file, 6.9, -0.92, True)

    x_NR   = np.asarray(scat_fit_mean_NR["s_prime"])
    x_NR_s   = np.asarray(scat_fit_spl_NR["s_prime"])
    sigma_NR_m   = np.asarray(scat_fit_mean_NR["sigma_prime"])
    sigma_NR_s   = np.asarray(scat_fit_spl_NR["sigma_prime"])

    info, info_nf, fit_param_mean_R, fit_param_spl_R, scat_fit_mean_R, scat_fit_spl_R, scat_nf_mean_R, scat_nf_spl_R = get_data(h5file, 7.05, -0.867, True)

    x_R   = np.asarray(scat_fit_mean_R["s_prime"])
    x_R_s   = np.asarray(scat_fit_spl_R["s_prime"])
    sigma_R_m   = np.asarray(scat_fit_mean_R["sigma_prime"])
    sigma_R_s   = np.asarray(scat_fit_spl_R["sigma_prime"])
    
    length = len(x_NR_s[0])
    
    sarr = np.linspace(slow-0.1,shigh,200)
    p2arr = [p2_s(s) for s in sarr]
    fit_model = fm.ERE_0_model
    
    fit_param_NR_m = np.asarray([fit_param_mean_NR[fp] for fp in fit_model.param_names])
    yarr_NR_m = np.asarray([fit_model.model(x,*fit_param_NR_m) for x in p2arr])
    sigma_NR_m = np.asarray([sigma_of_p3cotPS_prime(sarr[i],yarr_NR_m[i]) for i in range(len(sarr))])

    fit_param_s = np.transpose(np.asarray([fit_param_spl_NR[fp] for fp in fit_model.param_names]))
    yarr_NR_s = np.asarray([[fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))] for x in p2arr])
    sigma_NR_s = np.asarray([sorted([sigma_of_p3cotPS_prime(sarr[i],yarr_NR_s[i][j]) for j in range(len(yarr_NR_s[0]))]) for i in range(len(sarr))])

    sigma_NR_mean = [sigma_NR_s[i][length//2-1] for i in range(len(sarr))]
    sigma_NR_em = [sigma_NR_s[i][math.floor(length*(1-num_perc)/2)] for i in range(len(sarr))]
    sigma_NR_ep = [sigma_NR_s[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(sarr))]

    ax.plot(sarr,sigma_NR_mean, color = styles.c_10_non_res, label = r"$\boldsymbol{10}:\,\beta=6.9,\,am_0=-0.92$")
    ax.fill_between(sarr, sigma_NR_em, sigma_NR_ep, alpha = 0.5, color = styles.c_10_non_res)

    #####################################################################################################################    

    fit_model = fm.BW_I_model


    fit_param_R_m = np.asarray([fit_param_mean_R[fp] for fp in fit_model.param_names])
    yarr_R_m = np.asarray([fit_model.model(x,*fit_param_R_m) for x in sarr])
    sigma_R_m = np.asarray([sigma_of_p3cotPS_Ecm_prime(sarr[i],yarr_R_m[i]) for i in range(len(sarr))])

    fit_param_s = np.transpose(np.asarray([fit_param_spl_R[fp] for fp in fit_model.param_names]))
    yarr_R_s = np.asarray([[fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))] for x in sarr])
    sigma_R_s = np.asarray([sorted([sigma_of_p3cotPS_Ecm_prime(sarr[i],yarr_R_s[i][j])for j in range(len(yarr_R_s[0]))]) for i in range(len(sarr))])

    sigma_R_mean = [sigma_R_s[i][length//2-1] for i in range(len(sarr))]
    sigma_R_em = [sigma_R_s[i][math.floor(length*(1-num_perc)/2)] for i in range(len(sarr))]
    sigma_R_ep = [sigma_R_s[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(sarr))]

    ax.plot(sarr,sigma_R_mean, color = styles.c_10_res, label = r"$\boldsymbol{10}:\,\beta=7.05,\,am_0=-0.867$")
    ax.fill_between(sarr, sigma_R_em, sigma_R_ep, alpha = 0.5, color = styles.c_10_res)

    ###############################################################

    res,  res_sample = read_from_hdf("scattering_Fig5.3_b6.900_m-0.920")

    a_14_s = np.asarray(res_sample["a2"][:,0])
    b_14_s = np.asarray(res_sample["b2"][:,0])

    length = len(a_14_s)

    pcot_PS_14_s = np.asarray([[UTE(p2arr[j], a_14_s[i], b_14_s[i]) for i in range(len(a_14_s))] for j in range(len(p2arr))])

    sigma_14_s = np.asarray([sorted([sigma_of_pcotPS_prime(sarr[i], pcot_PS_14_s[i,j]) for j in range(len(a_14_s))]) for i in range(len(p2arr))])

    sigma_14_mean = [sigma_14_s[i][length//2-1] for i in range(len(sarr))]
    sigma_14_em = [sigma_14_s[i][math.floor(length*(1-num_perc)/2)] for i in range(len(sarr))]
    sigma_14_ep = [sigma_14_s[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(sarr))]

    ax.plot(sarr,sigma_14_mean, color = styles.c_14, label = r"$\boldsymbol{14}:\,\beta=6.9,\,am_0=-0.92$")
    ax.fill_between(sarr, sigma_14_em, sigma_14_ep, alpha = 0.5, color = styles.c_14)

    ax.set_xlabel(r"$s/m_\pi^2$")
    ax.set_ylabel(r"$\sigma_1 m_\pi^2$")
    xticks = np.linspace(4,11,8)
    yticks = np.linspace(0,100,6)
    ax.set_xticks(xticks, [r"$%i$"%x for x in xticks])
    ax.set_yticks(yticks, [r"$%i$"%x for x in yticks])


    ax.legend(loc='upper right')
    plt.savefig(op.join(PLTDIR, "sigma_comb.pdf"), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def plot_sigma_units(h5file,show=False):
    fig, ax = plt.subplots(figsize=(figs1,figs2))

    mDM = 100 # MeV
    Ecm_conv = 1/mDM
    sigma_conv = mDM**3/218426#/(10**(10))

    plt.grid()
    slow, shigh = 4+0.0001,11.1

    info, info_nf, fit_param_mean_NR, fit_param_spl_NR, scat_fit_mean_NR, scat_fit_spl_NR, scat_nf_mean_NR, scat_nf_spl_NR = get_data(h5file, 6.9, -0.92, True)

    x_NR   = np.asarray(scat_fit_mean_NR["s_prime"])
    x_NR_s   = np.asarray(scat_fit_spl_NR["s_prime"])
    sigma_NR_m   = np.asarray(scat_fit_mean_NR["sigma_prime"])
    sigma_NR_s   = np.asarray(scat_fit_spl_NR["sigma_prime"])

    info, info_nf, fit_param_mean_R, fit_param_spl_R, scat_fit_mean_R, scat_fit_spl_R, scat_nf_mean_R, scat_nf_spl_R = get_data(h5file, 7.05, -0.867, True)

    x_R   = np.asarray(scat_fit_mean_R["s_prime"])
    x_R_s   = np.asarray(scat_fit_spl_R["s_prime"])
    sigma_R_m   = np.asarray(scat_fit_mean_R["sigma_prime"])
    sigma_R_s   = np.asarray(scat_fit_spl_R["sigma_prime"])
    
    length = len(x_NR_s[0])
    
    sarr = np.linspace(slow,shigh,200)
    ECMarr = [np.sqrt(s)*mDM for s in sarr]
    p2arr = [p2_s(s) for s in sarr]
    fit_model = fm.ERE_0_model
    
    fit_param_NR_m = np.asarray([fit_param_mean_NR[fp] for fp in fit_model.param_names])
    yarr_NR_m = np.asarray([fit_model.model(x,*fit_param_NR_m) for x in p2arr])
    sigma_NR_m = np.asarray([sigma_of_p3cotPS_prime(sarr[i],yarr_NR_m[i]) for i in range(len(sarr))])

    fit_param_s = np.transpose(np.asarray([fit_param_spl_NR[fp] for fp in fit_model.param_names]))
    yarr_NR_s = np.asarray([[fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))] for x in p2arr])
    sigma_NR_s = np.asarray([sorted([sigma_of_p3cotPS_prime(sarr[i],yarr_NR_s[i][j]) for j in range(len(yarr_NR_s[0]))]) for i in range(len(sarr))])

    sigma_NR_mean = [sigma_NR_s[i][length//2-1]/sigma_conv for i in range(len(sarr))]
    sigma_NR_em = [sigma_NR_s[i][math.floor(length*(1-num_perc)/2)]/sigma_conv for i in range(len(sarr))]
    sigma_NR_ep = [sigma_NR_s[i][math.ceil(length*(1+num_perc)/2)]/sigma_conv for i in range(len(sarr))]

    ax.plot(ECMarr,sigma_NR_mean, color = styles.c_10_non_res, label = r"$\boldsymbol{10}:\,\beta=6.9,\,am_0=-0.92$", ls = "solid")
    ax.fill_between(ECMarr, sigma_NR_em, sigma_NR_ep, alpha = 0.5, color = styles.c_10_non_res)

    #####################################################################################################################    

    fit_model = fm.BW_I_model


    fit_param_R_m = np.asarray([fit_param_mean_R[fp] for fp in fit_model.param_names])
    yarr_R_m = np.asarray([fit_model.model(x,*fit_param_R_m) for x in sarr])
    sigma_R_m = np.asarray([sigma_of_p3cotPS_Ecm_prime(sarr[i],yarr_R_m[i]) for i in range(len(sarr))])

    fit_param_s = np.transpose(np.asarray([fit_param_spl_R[fp] for fp in fit_model.param_names]))
    yarr_R_s = np.asarray([[fit_model.model(x,*fit_param_s[i]) for i in range(len(fit_param_s))] for x in sarr])
    sigma_R_s = np.asarray([sorted([sigma_of_p3cotPS_Ecm_prime(sarr[i],yarr_R_s[i][j])for j in range(len(yarr_R_s[0]))]) for i in range(len(sarr))])

    sigma_R_mean = [sigma_R_s[i][length//2-1]/sigma_conv for i in range(len(sarr))]
    sigma_R_em = [sigma_R_s[i][math.floor(length*(1-num_perc)/2)]/sigma_conv for i in range(len(sarr))]
    sigma_R_ep = [sigma_R_s[i][math.ceil(length*(1+num_perc)/2)]/sigma_conv for i in range(len(sarr))]

    # ax.plot(sarr,sigma_R_m, color = "red")
    ax.plot(ECMarr,sigma_R_mean, color = styles.c_10_res, label = r"$\boldsymbol{10}:\,\beta=7.05,\,am_0=-0.867$", ls = "dashed")
    ax.fill_between(ECMarr, sigma_R_em, sigma_R_ep, alpha = 0.5, color = styles.c_10_res)

    ###############################################################

    res,  res_sample = read_from_hdf("scattering_Fig5.3_b6.900_m-0.920")

    a_14_s = np.asarray(res_sample["a2"][:,0])
    b_14_s = np.asarray(res_sample["b2"][:,0])

    length = len(a_14_s)

    pcot_PS_14_s = np.asarray([[UTE(p2arr[j], a_14_s[i], b_14_s[i]) for i in range(len(a_14_s))] for j in range(len(p2arr))])

    sigma_14_s = np.asarray([sorted([sigma_of_pcotPS_prime(sarr[i], pcot_PS_14_s[i,j]) for j in range(len(a_14_s))]) for i in range(len(p2arr))])

    sigma_14_mean = [sigma_14_s[i][length//2-1]/sigma_conv for i in range(len(sarr))]
    sigma_14_em = [sigma_14_s[i][math.floor(length*(1-num_perc)/2)]/sigma_conv for i in range(len(sarr))]
    sigma_14_ep = [sigma_14_s[i][math.ceil(length*(1+num_perc)/2)]/sigma_conv for i in range(len(sarr))]

    ax.plot(ECMarr,sigma_14_mean, color = styles.c_14, label = r"$\boldsymbol{14}:\,\beta=6.9,\,am_0=-0.92$", ls = "dashdot")
    ax.fill_between(ECMarr, sigma_14_em, sigma_14_ep, alpha = 0.5, color = styles.c_14)

    ax.set_xlabel(r"$E_{\text{cm}}\,[\text{MeV}]$", fontsize=styles.fontsize)
    ax.set_ylabel(r"$\sigma / m_{\text{DM}}\,[\text{cm}^2/g]$", fontsize=styles.fontsize)
    xticks = np.linspace(200,320,7)
    yticks = np.linspace(0,35,8)
    ax.set_xticks(xticks, [r"$%i$"%(x) for x in xticks])
    ax.set_yticks(yticks, [r"$%i$"%(x) for x in yticks])
    ax.set_xlim([200,320])
    ax.set_ylim([0,36])

    handles, labels = ax.get_legend_handles_labels()
    ins = [2,0,1]
    handles = [handles[ins[0]],handles[ins[1]],handles[ins[2]]]
    labels = [labels[ins[0]],labels[ins[1]],labels[ins[2]]]
    ax.legend(handles, labels, title=r'$m_{\text{DM}}=%i\,\text{MeV}$'%(mDM), loc="upper right", fontsize=styles.fontsize)

    plt.savefig(op.join(PLTDIR, "sigma_comb_units.pdf"), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file  = args[2]

    plot_sigma_units(h5file, False)
    # plot_sigma(h5file, False)
