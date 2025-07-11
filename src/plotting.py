import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import math
import fit_scatter
import os.path as op
import os

########################################## Plot energy levels in lattice units ##########################################

def color(d2):
    colors = ["orange", "green", "blueviolet"]
    return colors[d2-1]

def E_pipi(mpi,p12,p22,L):
    return np.sqrt(mpi**2+(2*np.pi/L)**2*p12)+np.sqrt(mpi**2+(2*np.pi/L)**2*p22)  

def E_rho(mrho,p2,L):
    return np.sqrt(mrho**2+(2*np.pi/L)**2*p2)
    
def marker(lv):                     # maybe to be replaced by input file
    markers = ["*", "o", "^"]
    return markers[lv-1]

def get_data_E_L(name, beta, m0, num_lv = 2):
    NLs,NL_invs,aEs,aE_ms,aE_ps,d2s,lvs = [[],[],[],[],[],[],[]]

    with h5py.File(op.join(OUTDIR, "isospin1_scattering"+name+".hdf5"),"r") as hfile:
        for key in hfile:
            if str(beta) in key and str(m0) in key:
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

def plot_E_L(name,beta,m0,levels=False,outname=None,show=False):
    mpi, mrho, d2s, NLs, NL_invs, En, En_m_err, En_p_err, lvs = get_data_E_L(name, beta, m0)
    
    for i in range(len(En)):
        plt.errorbar([NL_invs[i],],y=[En[i],],yerr=[[En_m_err[i],],[En_p_err[i],]], solid_capstyle="projecting", capsize=5, ls="", color = color(d2s[i]), marker = marker(lvs[i]))   
    plt.axhline(mrho,c="red", ls="dotted", label = "$m_\\rho$")
    plt.axhline(2*mpi,c="black",label = r"2$m_\pi$")
    plt.axhline(4*mpi,c="black",label = r"4$m_\pi$")
    plt.grid()
    plt.title("$\\beta$ = %f, $m_0$ = %f"%(beta,m0))
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
    plt.plot([0,0],[0,0],c="grey", label = "non-int")
    plt.xlim([1/40,1/13])
    plt.ylim([0.3,2])

    for i in range(1,max(lvs)+2):
        plt.errorbar([0,],y=[0,],yerr=[[0,],[0,]], solid_capstyle="projecting", capsize=5, ls="", color = color(i), marker = "o", label = "|P|=%i"%(i))
    for i in range(1,max(lvs)+2):
        plt.scatter(x=[0,],y=[0,], color = "grey", marker = marker(i), label = "lv=%i"%(i))

    plt.legend(loc='center right', bbox_to_anchor=(1.24, 0.5))

    plt.xticks([1/14,1/16,1/20,1/24,1/36],["1/14","1/16","1/20","1/24","1/36"])
    plt.xlabel("1/$N_L$")
    plt.ylabel("a$E$")
    if outname == None:    
        plt.savefig(op.join(PLTDIR, "E_L_b%f_m0%f_levels_%r.pdf"%(beta,m0,levels)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR, "E_L_"+outname+"_levels_%r.pdf"%levels), bbox_inches='tight')
    if show:
        plt.show()
    plt.clf()

########################################## Plot com energies unitless ##########################################

def get_data_E_CM_L(name, beta, m0, num_lv = 2):
    NLs,NL_invs,aEs,aE_ms,aE_ps,d2s,lvs = [[],[],[],[],[],[],[]]

    with h5py.File(op.join(OUTDIR, "isospin1_scattering"+name+".hdf5"),"r") as hfile:
        for key in hfile:
            if str(beta) in key and str(m0) in key:
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

def plot_E_CM_L(name,beta,m0,levels=False,outname=None,show=False):
    mpi, mrho, d2s, NLs, NL_invs, En, En_m_err, En_p_err, lvs = get_data_E_L(name, beta, m0)

    ECMs = [np.sqrt(En[i]**2-(2*np.pi/NLs[i])**2*d2s[i]) for i in range(len(En))]
    ECM_errms = [abs(np.sqrt((En[i]-En_m_err[i])**2-(2*np.pi/NLs[i])**2*d2s[i])-ECMs[i])/mpi  for i in range(len(En))]
    ECM_errps = [abs(np.sqrt((En[i]+En_p_err[i])**2-(2*np.pi/NLs[i])**2*d2s[i])-ECMs[i])/mpi  for i in range(len(En))]
    ECMs = [ECMs[i]/mpi for i in range(len(ECMs))]
    
    for i in range(len(ECMs)):
        plt.errorbar([NL_invs[i],],y=[ECMs[i],],yerr=[[ECM_errms[i],],[ECM_errps[i],]], solid_capstyle="projecting", capsize=5, ls="", color = color(d2s[i]), marker = marker(lvs[i]))   
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

    for i in range(1,4):
        plt.errorbar([0,],y=[0,],yerr=[[0,],[0,]], solid_capstyle="projecting", capsize=5, ls="", color = color(i), marker = "o", label = "|P|=%i"%(i))
    for i in range(1,4):
        plt.scatter(x=[0,],y=[0,], color = "grey", marker = marker(i), label = "lv=%i"%(i))

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

def color_NL(NL):
    if NL == 14:
        return "red"
    elif NL == 16:
        return "green"
    elif NL == 24:
        return "c"
    elif NL == 36:
        return "m"

def ls_P(dvec):
    if list(dvec) == [0,0,1]:
        return "solid"
    elif list(dvec) == [1,1,0] or list(dvec) == [0,1,1]:
        return "dashed"
    elif list(dvec) == [1,1,1]:
        return "dashdot"

def ms_P(dvec):
    if list(dvec) == [0,0,1]:
        return "*"
    elif list(dvec) == [1,1,0] or list(dvec) == [0,1,1]:
        return "o"
    elif list(dvec) == [1,1,1]:
        return "^"

def delete_steps(arr, sign = 1, delete=True):
    if delete:
        for i in range(len(arr)-1):
            if arr[i+1] < sign*arr[i]: 
                arr[i] = np.nan
        return arr
    else:
        return arr

def get_data_p3cotPS(name, beta, m0):
    res = {}
    res_spl = {}

    with h5py.File(op.join(OUTDIR, "isospin1_fit_scatter"+name+".hdf5"),"r") as hfile:
        for key in hfile["fit_scatter_b%f_m%f"%(beta,m0)]["mean"]:
            res[key] = hfile["fit_scatter_b%f_m%f"%(beta,m0)]["mean"][key][()]
            res_spl[key] = hfile["fit_scatter_b%f_m%f"%(beta,m0)]["sample"][key][()]
    return res, res_spl

def plot_p3cotPS(name,beta,m0,fit=False,outname=None,show=False):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    plt.grid()
    res, res_smp = get_data_p3cotPS(name, beta, m0)
    xlim = [0,3]
    ax.set_xlim(xlim)

    plt.xlabel(r"$p^{\star^2}/m_\pi^2$")
    x_plot = res["p2star_prime"]
    x_plot_sam = np.transpose(res_smp["p2star_prime"])
    y_plot = np.real(res["p3cotPS_prime"])
    y_plot_sam = np.transpose(np.real(res_smp["p3cotPS_prime"]))
    plt.ylabel(r"$p^3\, \cot(\delta)/m_\pi^3$")
    
    length = len(x_plot_sam[0])
    num_perc = math.erf(1/np.sqrt(2))
    N_Ls = res["N_L"]
    dvecs = res["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = [np.dot(d,d) for d in dvecs]

    for i in  range(14):
        ax.scatter(x_plot[i],y_plot[i], label = "|P|=%i, NL=%i"%(d2s[i],N_Ls[i]), color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]), marker = ms_P(dvecs[i]),s=60)
        sorted_indices = np.argsort(x_plot_sam[i])
        ax.plot(x_plot_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(y_plot_sam[i][sorted_indices],delete=True)[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = color_NL(N_Ls[i]), ls = ls_P(dvecs[i]))
        
    xarr = np.linspace(xlim[0], xlim[1])
    
    if fit:
        a1_1 = res["a1_1"]
        r1_1 = res["r1_1"]

        a1_1_smp = res_smp["a1_1"]
        r1_1_smp = res_smp["r1_1"]

        yarr = [fit_scatter.ERE_1(x,a1_1,r1_1) for x in xarr]
        yarr_smp = [sorted([fit_scatter.ERE_1(x,a1_1_smp[i],r1_1_smp[i]) for i in range(len(a1_1_smp))]) for x in xarr]

        yarr_m = [yarr_smp[i][math.floor(length*(1-num_perc)/2)] for i in range(len(yarr_smp))]
        yarr_p = [yarr_smp[i][math.ceil(length*(1+num_perc)/2)] for i in range(len(yarr_smp))]

        plt.plot(xarr,yarr, color = "blue")
        plt.fill_between(xarr, yarr_m, yarr_p, alpha = 0.3, color = "blue")

    ax.legend(loc='center right', bbox_to_anchor=(1.35, 0.5))
    if outname == None:    
        plt.savefig(op.join(PLTDIR, "p3cotPS_b%f_m0%f_fit_%r.pdf"%(beta,m0,fit)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR, "p3cotPS_"+outname+"_fit_%r.pdf"%fit), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

########################################## Plot sigma1 ##########################################

def sigma_ERE_s_wave(s, a, r):
    p2 = s/4-1
    cot_PS = (-1/a+p2*r/2)/np.sqrt(p2)
    return 4*np.pi/(cot_PS**2+1)/p2
    
def sigma_of_P3cotPS(P3cotPS, p2):
    if p2 == 0:
        return 0
    else:
        cot_PS = P3cotPS/(p2**(3/2))
        return 4*np.pi*3/(cot_PS**2+1)/p2

def plot_sigma_1(name,beta,m0,fit=False,outname=None,show=False):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    plt.grid()
    res, res_smp = get_data_p3cotPS(name, beta, m0)

    xlim = [4,6.5]

    plt.xlabel(r"$s/m_\pi^2$")
    x_plot_sam = np.transpose(res_smp["s_prime"])

    plt.ylabel(r"$\sigma m_\pi^2$")        
    
    length = len(x_plot_sam[0])
    num_perc = math.erf(1/np.sqrt(2))

    sarr = np.linspace(xlim[0],xlim[1],5000)
    p2arr = [x/4-1 for x in sarr]

    yarr_14 = [sigma_ERE_s_wave(s, 0.52, 6.7) for s in sarr]

    plt.plot(sarr, yarr_14, color = "red", label = "2405.06506")
    a1_1 = res["a1_1"]
    r1_1 = res["r1_1"]

    print("smax = ", 4*(1+1/r1_1**2))
    a1_1_smp = res_smp["a1_1"]
    r1_1_smp = res_smp["r1_1"]

    yarr = [sigma_of_P3cotPS(fit_scatter.ERE_1(x,a1_1,r1_1), x) for x in p2arr]
    plt.plot(sarr,yarr,label="This work")

    yarr_smp = []
    for p2 in p2arr:
        p3cotPS_smp = [fit_scatter.ERE_1(p2,a1_1_smp[i],r1_1_smp[i]) for i in range(len(a1_1_smp))]
        sorted_indices = np.argsort(p3cotPS_smp)[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)]
        yarr_smp.append(np.asarray([sigma_of_P3cotPS(fit_scatter.ERE_1(p2,a1_1_smp[i],r1_1_smp[i]), p2) for i in range(len(a1_1_smp))])[sorted_indices])


    yarr_m = [min(yarr_smp[i]) for i in range(len(yarr_smp))]
    yarr_p = [max(yarr_smp[i]) for i in range(len(yarr_smp))]

    plt.fill_between(sarr, yarr_m, yarr_p, alpha = 0.3)

    if outname == None:    
        plt.savefig(op.join(PLTDIR,"scattering/sigma1_b%f_m0%f_fit_%r.pdf"%(beta,m0,fit)), bbox_inches='tight')
    else:    
        plt.savefig(op.join(PLTDIR,"sigma1_"+outname+"_fit_%r.pdf"%fit), bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

if __name__ == "__main__":

    # avod hard-coding of names outside of main
    OUTDIR = "../output/scattering/"
    PLTDIR = "../output/plots/scattering/"

    # create directories if they do not exist
    os.makedirs("../output/scattering", exist_ok=True)
    os.makedirs("../output/plots/scattering", exist_ok=True)

    # name = "_evp_deriv_false"
    name = "_evp_deriv_true"
    plot_E_L(name,6.9,-0.92,True,outname="non_res")
    plot_E_L(name,6.9,-0.92,False,outname="non_res")
    plot_E_L(name,7.05,-0.863,True,outname="close_res")
    plot_E_L(name,7.05,-0.863,False,outname="close_res")
    plot_E_L(name,7.05,-0.867,True,outname="res")
    plot_E_L(name,7.05,-0.867,False,outname="res")

    plot_E_CM_L(name,6.9,-0.92,True,outname="non_res")
    plot_E_CM_L(name,6.9,-0.92,False,outname="non_res")
    plot_E_CM_L(name,7.05,-0.863,True,outname="close_res")
    plot_E_CM_L(name,7.05,-0.863,False,outname="close_res")
    plot_E_CM_L(name,7.05,-0.867,True,outname="res")
    plot_E_CM_L(name,7.05,-0.867,False,outname="res")
    
    plot_p3cotPS(name,6.9,-0.92,True,outname="non_res",show=False)
    plot_p3cotPS(name,7.05,-0.863,True,outname="close_res",show=False)
    plot_p3cotPS(name,7.05,-0.867,True,outname="res",show=False)
    
    plot_sigma_1(name,6.9,-0.92,True,outname="non_res",show=False)
    plot_sigma_1(name,7.05,-0.863,True,outname="close_res",show=False)
    plot_sigma_1(name,7.05,-0.867,True,outname="res",show=False)
