import numpy as np
import h5py
import os.path as op
import os
from scipy.optimize import curve_fit

def curve_fit_try(func, x, y, num_res):
    try:
        popt, pcov = curve_fit(func, x, y)
        return popt
    except:
        return np.zeros(num_res)

def ERE_0(p2, a1_0):
    return -1/a1_0**3+0*p2
def ERE_1(p2, a1_1, r1_1):
    return -1/a1_1**3+p2/(2*r1_1)
def ERE_2(p2, a1_2, r1_2, c1_2):
    return -1/a1_2**3+p2/(2*r1_2)+c1_2*p2**2
def ERE_1_lin(p2, a, b):
    return a+b*p2
def ERE_2_lin(p2, a, b, c):
    return a+b*p2+c*p2**2
def ERE_fit(p3cotPS, p2):                                   # all primed
    result = {}
    popt = curve_fit_try(ERE_0, p2, p3cotPS,1)
    result["a1_0"] = popt[0]
    popt = curve_fit_try(ERE_1_lin, p2, p3cotPS,2)
    result["a1_1"]=-np.sign(popt[0])*np.power(1/abs(popt[0]),1/3.)
    result["r1_1"]=1/(2*popt[1])
    popt = curve_fit_try(ERE_2_lin, p2, p3cotPS,3)
    result["a1_2"]=-np.sign(popt[0])*np.power(1/abs(popt[0]),1/3.)
    result["r1_2"]=1/(2*popt[1])
    result["c1_2"]=popt[2]
    return result

def RES_Drach(ECM, m_R, gVPP2):
    return 6*np.pi*(m_R**2-ECM**2)/gVPP2
def ECM_p(p):
    return 2*np.sqrt(1+p**2)
def RES_Alex_BWI(p, m_R, gVPP2):
    ECM = ECM_p(p)
    return ECM*gVPP2*p**3/(6*np.pi*ECM**2*(m_R**2-ECM**2))
def RES_Alex_BWII(p, m_R, gVPP2, r0):
    ECM = ECM_p(p)
    return ECM*gVPP2*p**3*(1+np.sqrt((m_R/2)**2-1))/(6*np.pi*ECM**2*(m_R**2-ECM**2)*(1+(p*r0)**2))
def RES_fit(p3cotPS_ECM, ECM, tan_PS, p):                                   # all primed
    result = {}
    popt = curve_fit_try(RES_Drach, ECM, p3cotPS_ECM,2)
    result["m_R_D"], result["gVPP2_D"] = popt
    popt = curve_fit_try(RES_Alex_BWI, p, tan_PS,2)
    result["m_R_BWI"], result["gVPP2_BWI"] = popt
    popt = curve_fit_try(RES_Alex_BWII, p, tan_PS,2)
    result["m_R_BWII"], result["gVPP2_BWII"], result["r0_BWII"] = popt
    return result

def get_fits(res, res_spl):
    res_tmp={}
    res_spl_tmp={}
    win_ind = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    for key, val in ERE_fit(res["p3cotPS_prime"][win_ind],res["p2star_prime"][win_ind]).items():
        res_tmp[key] = val
    for key, val in RES_fit(res["p3cotPS_Ecm_prime"][win_ind],res["E_cm_prime"][win_ind],res["tan_PS"][win_ind],res["pstar_prime"][win_ind]).items():
        res_tmp[key] = val
        
    for key in res_tmp.keys():
        res_spl_tmp[key] = []

    for i in range(len(res_spl["p3cotPS_prime"])):
        for key, val in ERE_fit(res_spl["p3cotPS_prime"][i][win_ind],res_spl["p2star_prime"][i][win_ind]).items():
            res_spl_tmp[key].append(val)
        for key, val in RES_fit(res_spl["p3cotPS_Ecm_prime"][i][win_ind],res_spl["E_cm_prime"][i][win_ind],res_spl["tan_PS"][i][win_ind],res_spl["pstar_prime"][i][win_ind]).items():
            res_spl_tmp[key].append(val)
    return res_tmp, res_spl_tmp

def fit_one_phaseshift(name, beta, m0):
    res_calc = {}
    res_spl_calc = {}

    with h5py.File(op.join(OUTDIR,"isospin1_scattering"+name+".hdf5"),"r") as hfile:
        for ens in hfile:
            if str(beta) in ens and str(m0) in ens:
                for P in hfile[ens]:
                    if P[0] == "p":
                        dvec = [int(P[2]),int(P[4]),int(P[6])]
                        for irrep in hfile[ens][P]:
                            for lv in hfile[ens][P][irrep]:
                                if lv[:2] == "lv":
                                    if res_calc == {}:
                                        num_resample = hfile[ens][P][irrep][lv]["info"]["num_resample"][()]
                                        for key in hfile[ens][P][irrep][lv]["mean"]:
                                            res_calc[key] = []
                                            res_spl_calc[key] = []
                                    for tmp in hfile[ens][P][irrep][lv]["mean"]:
                                        res_calc[tmp].append(hfile[ens][P][irrep][lv]["mean"][tmp][()])
                                        res_spl_calc[tmp].append(hfile[ens][P][irrep][lv]["sample"][tmp][()])

    for key in res_calc:
        res_calc[key] = np.asarray(res_calc[key])
        res_spl_calc[key] = np.transpose(np.asarray(res_spl_calc[key]))

    res_fit, res_spl_fit = get_fits(res_calc,res_spl_calc)
    with h5py.File(op.join(OUTDIR,"isospin1_fit_scatter"+name+".hdf5"),"a") as hfile:
        for key, val in res_fit.items():
            hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"mean/"+key, data = val)
        for key, val in res_spl_fit.items():
            hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"sample/"+key, data = val)
        for key, val in res_calc.items():
            hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"mean/"+key, data = val)
        for key, val in res_spl_calc.items():
            hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"sample/"+key, data = val)

def fit_all_phase_shifts(name):
    fit_one_phaseshift(name,6.9,-0.92)
    fit_one_phaseshift(name,7.05,-0.863)
    fit_one_phaseshift(name,7.05,-0.867)


if __name__ == "__main__":
    # avod hard-coding of names outside of main
    # create directories if they do not exist
    OUTDIR = "../data_assets/scattering/"
    os.makedirs("../data_assets/scattering", exist_ok=True)

    # name = "_evp_deriv_false"
    name = "_evp_deriv_true"
    fit_all_phase_shifts(name)