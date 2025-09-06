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
import fit_models as fm

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

def get_data(h5file_scatter_fit, beta, m0, fit):
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
                                    # print(ens+P+irrep+lv)
                                    if hfile[ens][P][irrep][lv]["fit"][()] or not fit:
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

def plot_E_CM_L(h5file,beta,m0,levels=False,outname=None,show=False):
    
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, beta, m0, False)
    
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
        plt.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]))   
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
    plt.ylim([1,6])

    for tmp in [[None,0,"T1",0],[None,1,"E",0],[None,2,"B1",0],[None,3,"E",0],[None,1,"A1",0],[None,1,"A1",1],[None,2,"A1",0],[None,2,"A1",1],[None,3,"A1",0],[None,3,"A1",1]]:
        plt.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = "p=%i, %s, lv=%i"%(tmp[1],tmp[2],tmp[3]))
    for tmp in [[14,None,None,None],[16,None,None,None],[20,None,None,None],[24,None,None,None],[36,None,None,None]]:
        plt.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(*tmp), label = "$N_L$=%i"%(tmp[0]))

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

def xlim_f(m0, xaxis="p2star_prime"):
    if xaxis == "p2star_prime":
        if m0 == -0.92:
            # return [0,0.3]
            return [0,0.6]
        elif m0 == -0.863:
            return [0,3]
        elif m0 == -0.867:
            return [0,3]
    elif xaxis == "s_prime":
        if m0 == -0.92:
            return [4,6]
        elif m0 == -0.863:
            return [4,16]
        elif m0 == -0.867:
            return [4,16]
    raise ValueError("x- or y-lim not defined for %s, %s"%(m0, xaxis))

def ylim_f(m0, yaxis="p3cotPS_prime"):
    if yaxis == "p3cotPS_prime":
        if m0 == -0.92:
            # return [0,2]
            return [-3,3]
        elif m0 == -0.863:
            return [-4,4]
        elif m0 == -0.867:
            return [-4,4]
    elif yaxis == "p3cotPS_Ecm_prime":
        if m0 == -0.92:
            return [0,2]
        elif m0 == -0.863:
            return [-2,2]
        elif m0 == -0.867:
            return [-2,2]
    elif yaxis == "sigma_prime":
        if m0 == -0.92:
            return [0,100]
        elif m0 == -0.863:
            return [0,100]
        elif m0 == -0.867:
            return [0,100]
    raise ValueError("x- or y-lim not defined for %s, %s"%(m0, yaxis))

def xlabel_f(xaxis):
    if xaxis == "p2star_prime":
        return r"$p^{\star^2}/m_\pi^2$"
    elif xaxis == "s_prime":
        return r"$s/m_\pi^2$"
    raise ValueError("Label not defined for %s"%(xaxis))

def ylabel_f(yaxis):
    if yaxis == "p3cotPS_prime":
        return r"$p^3\, \cot(\delta)/m_\pi^3$"
    elif yaxis == "p3cotPS_Ecm_prime":
        return r"$p^3\, \cot(\delta)/E_{cm}/m_\pi$"
    elif yaxis == "sigma_prime":
        return r"$\sigma_1 m_\pi^2$"
    raise ValueError("Label not defined for %s"%(yaxis))

# def p2_p2(p2):
#     return p2
def p2_s(s):
    return s/4-1

def x_axis_func(x, xaxis):
    if xaxis == "p2star_prime":
        return x
    elif xaxis == "s_prime":
        return p2_s(x)
    raise ValueError("x axis func not defined for %s"%(xaxis))

def p3_cot_PS_ECM(p2, p3cotPS):
    ECM = np.sqrt(4+p2)
    return p3cotPS/ECM
def sigma_1(p2, p3cotPS):
    return 12*np.pi*p2**2/(p2**3+p3cotPS**2)
def y_axis_func(p2, p3cotPS, yaxis):
    if yaxis == "p3cotPS_prime":
        return p3cotPS
    elif yaxis == "p3cotPS_Ecm_prime":
        return p3_cot_PS_ECM(p2,p3cotPS)
    elif yaxis == "sigma_prime":
        return sigma_1(p2,p3cotPS)
    raise ValueError("y axis func not defined for %s"%(yaxis))

def plot_any(h5file,beta,m0,xaxis="p2star_prime",yaxis="p3cotPS_prime",fit_model=None,outname=None,show=False):
    fit = fit_model != None
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    plt.grid()
    info, info_nf, fit_param_mean, fit_param_spl, scat_fit_mean, scat_fit_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, beta, m0, fit)
    xlim = xlim_f(m0,xaxis)
    if xlim == None:
        ax.set_xlim(auto=True)
    else: 
        ax.set_xlim(xlim)
    ylim = ylim_f(m0,yaxis)
    if ylim == None:
        ax.set_ylim(auto=True)
    else: 
        ax.set_ylim(ylim)

    x_m   = np.asarray(scat_fit_mean[xaxis])
    x_s   = np.asarray(scat_fit_spl[xaxis])
    y_m   = np.asarray(scat_fit_mean[yaxis])
    y_s   = np.asarray(scat_fit_spl[yaxis])
    if fit:
        x_nf_m       = np.asarray(scat_nf_mean[xaxis])
        x_nf_s   = np.asarray(scat_nf_spl[xaxis])
        y_nf_m       = np.asarray(scat_nf_mean[yaxis])
        y_nf_s   = np.asarray(scat_nf_spl[yaxis])

    xlabel = xlabel_f(xaxis)
    plt.xlabel(xlabel)
    ylabel = ylabel_f(yaxis)
    plt.ylabel(ylabel)
    
    length = len(x_s[0])
    
    N_Ls = [int(x) for x in scat_fit_mean["N_L"]]
    dvecs = scat_fit_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]
    lvs = info["lv"]
    irreps = info["irrep"]
    plot_args = list(zip(N_Ls,d2s,irreps,lvs))
    
    if fit:
        N_Ls_nf = [int(x) for x in scat_nf_mean["N_L"]]
        dvecs_nf = scat_nf_mean["dvec"]
        dvecs_nf = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs_nf]
        d2s_nf = [np.dot(d,d) for d in dvecs_nf]
        lvs_nf = info_nf["lv"]
        irreps_nf = info_nf["irrep"]
        plot_args_nf = list(zip(N_Ls_nf,d2s_nf,irreps_nf,lvs_nf))

    for i in  range(len(x_m)):
        ax.scatter(x_m[i],y_m[i], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]), marker = pf.marker(*plot_args[i]), s = 10*pf.ms(*plot_args[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
        sorted_indices = np.argsort(x_s[i])
        ax.plot(x_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = pf.color(*plot_args[i]), ls = pf.ls(*plot_args[i]))

    if fit:
        for i in  range(len(x_nf_m)):
            ax.scatter(x_nf_m[i],y_nf_m[i], color = "grey", ls = pf.ls(*plot_args_nf[i]), marker = pf.marker(*plot_args_nf[i]), s = 10*pf.ms(*plot_args_nf[i]))   #, label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i])
            sorted_indices = np.argsort(x_nf_s[i])
            ax.plot(x_nf_s[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_nf_s[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = "grey", ls = pf.ls(*plot_args_nf[i]))


    xarr = np.linspace(xlim[0]+1e-3, xlim[1], 200)
    
    if fit:
        p2_arr = [x_axis_func(x, xaxis) for x in xarr]
        fit_param_s = np.asarray([fit_param_spl[fp] for fp in fit_model.param_names])

        p3cotPS_s = [[fit_model.model(p2,*fit_param_s[:,i]) for i in range(len(fit_param_s[0]))] for p2 in p2_arr]

        # y_f_m = [y_axis_func(p2_arr[i],p3cotPS_m[i],yaxis) for i in range(len(p2_arr))]
        y_f_s = [sorted([y_axis_func(p2_arr[i],p3cotPS_s[i][j],yaxis) for j in range(len(p3cotPS_s[0]))]) for i in range(len(p2_arr))]

        y_f_m_m = [y_f_s[i][length//2-1] for i in range(len(p2_arr))]
        y_f_e_m = [y_f_s[i][math.floor(length*(1-num_perc)/2)] for i in range(len(p2_arr))]
        y_f_e_p = [y_f_s[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(p2_arr))]

        plt.plot(xarr,y_f_m_m, color = "blue")
        plt.fill_between(xarr, y_f_e_m, y_f_e_p, alpha = 0.3, color = "blue")

    if fit:
        p2_m = np.asarray(scat_fit_mean["p2star_prime"])
        
        # print(len(p2_m), len(p2_e))
        p3cotPS_mean = np.asarray(scat_fit_mean["p3cotPS_prime"])
        p3cotPS_s = np.asarray(scat_fit_spl["p3cotPS_prime"])
        p3cotPS_e = [abs((p3cotPS_s[i][math.ceil(length*(1+num_perc)/2)]-p3cotPS_s[i][math.floor(length*(1-num_perc)/2)]))/2 for i in range(len(p2_m))]
        fit_param_m = np.asarray([fit_param_mean[fp] for fp in fit_model.param_names])
        p3cotPS_pred = [fit_model.model(p2, *fit_param_m) for p2 in p2_m]

        chi2 = np.sum(((p3cotPS_mean - p3cotPS_pred) / p3cotPS_e) ** 2)
        ndof = len(p3cotPS_mean) - len(fit_param_m)  # degrees of freedom
        chi2_ndof = chi2 / ndof
        print(beta, m0, fit_model.name, "chi2 = %f,\tdof = %f,\tchi2/dof = %f"%(chi2,ndof,chi2_ndof))
        plt.plot([-1,-1],[-1,-1], color = "grey", label = "not fitted")

    for tmp in [[None,0,"T1",0],[None,1,"E",0],[None,2,"B1",0],[None,3,"E",0],[None,1,"A1",0],[None,1,"A1",1],[None,2,"A1",0],[None,2,"A1",1],[None,3,"A1",0],[None,3,"A1",1]]:
        plt.scatter(x=[-1,],y=[-1,], color = pf.color(*tmp), marker = "o", label = "p=%i, %s, lv=%i"%(tmp[1],tmp[2],tmp[3]))
    for tmp in [[14,None,None,None],[16,None,None,None],[20,None,None,None],[24,None,None,None],[36,None,None,None]]:
        plt.scatter(x=[-1,],y=[-1,], color = "grey", marker = pf.marker(*tmp), label = "$N_L$=%i"%(tmp[0]))
    ax.legend(loc='center right', bbox_to_anchor=(1.35, 0.5))
    fit_str = "" if fit_model == None else "_fit_%s"%fit_model.name
    if outname == None:    
        plt.savefig(op.join(PLTDIR, "%s_%s__b%f_m0%f%s.pdf"%(yaxis,xaxis,beta,m0,fit_str)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR, "%s_%s__"%(yaxis,xaxis)+outname+"%s.pdf"%(fit_str)), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

if __name__ == "__main__":

    args = sys.argv
    PLTDIR = args[1]
    h5file  = args[2]

    os.makedirs(PLTDIR, exist_ok=True)

    plot_E_CM_L(h5file, 6.9, -0.92)
    plot_E_CM_L(h5file, 7.05, -0.863)
    plot_E_CM_L(h5file, 7.05, -0.867)

    plot_any(h5file, 6.9, -0.92, "p2star_prime", "p3cotPS_prime", None)
    plot_any(h5file, 7.05, -0.863, "p2star_prime", "p3cotPS_prime", None)
    plot_any(h5file, 7.05, -0.863, "p2star_prime", "p3cotPS_Ecm_prime", None)
    plot_any(h5file, 7.05, -0.867, "p2star_prime", "p3cotPS_Ecm_prime", None)
    plot_any(h5file, 7.05, -0.863, "s_prime", "p3cotPS_Ecm_prime", None)
    plot_any(h5file, 7.05, -0.867, "s_prime", "p3cotPS_Ecm_prime", None)

    plot_any(h5file, 6.9, -0.92, "p2star_prime", "p3cotPS_prime", fm.ERE_0_model)
    plot_any(h5file, 6.9, -0.92, "p2star_prime", "p3cotPS_prime", fm.ERE_1_model)
    plot_any(h5file, 6.9, -0.92, "p2star_prime", "p3cotPS_prime", fm.ERE_2_model)
    plot_any(h5file, 6.9, -0.92, "p2star_prime", "p3cotPS_prime", fm.NR_I_model)
    # plot_any(h5file, 6.9, -0.92, "p2star_prime", "p3cotPS_prime", fm.NR_II_model)
    plot_any(h5file, 7.05, -0.863, "s_prime", "p3cotPS_Ecm_prime", fm.BW_I_model)
    plot_any(h5file, 7.05, -0.863, "s_prime", "p3cotPS_Ecm_prime", fm.BW_II_model)
    plot_any(h5file, 7.05, -0.867, "s_prime", "p3cotPS_Ecm_prime", fm.BW_I_model)
    plot_any(h5file, 7.05, -0.867, "s_prime", "p3cotPS_Ecm_prime", fm.BW_II_model)

    plot_any(h5file, 6.9, -0.92, "s_prime", "sigma_prime", fm.NR_I_model)
    plot_any(h5file, 6.9, -0.92, "s_prime", "sigma_prime", fm.ERE_0_model)
    plot_any(h5file, 7.05, -0.863, "s_prime", "sigma_prime", fm.BW_I_model)
    plot_any(h5file, 7.05, -0.863, "s_prime", "sigma_prime", fm.BW_II_model)
    plot_any(h5file, 7.05, -0.867, "s_prime", "sigma_prime", fm.BW_I_model)
    plot_any(h5file, 7.05, -0.867, "s_prime", "sigma_prime", fm.BW_II_model)