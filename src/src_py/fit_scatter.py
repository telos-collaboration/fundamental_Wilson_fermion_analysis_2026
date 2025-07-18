import numpy as np
import h5py
import sys
from scipy.optimize import curve_fit
from tqdm import tqdm

import warnings
warnings.simplefilter("error")


def curve_fit_try(func, x, y, num_res):
    # popt, pcov = curve_fit(func, x, y)
    # return popt
    try:
        # print("hey")
        # print(x)
        # print(y)
        popt, pcov = curve_fit(func, x, y)
        # print("\tyay")
        return popt
    except:
        # print("\t\tnay")
        return np.zeros(num_res)

def ERE_0(p2, a1_0):
    return -1/a1_0**3+0*p2
def ERE_1(p2, a1_1, r1_1):
    return 0 if a1_1 == 0 or r1_1 == 0 else  -1/a1_1**3+p2/(2*r1_1)
    # return result
def ERE_2(p2, a1_2, r1_2, c1_2):
    return -1/a1_2**3+p2/(2*r1_2)+c1_2*p2**2
def ERE_1_lin(p2, a, b):
    return a+b*p2
def ERE_2_lin(p2, a, b, c):
    return a+b*p2+c*p2**2
    # return result
def ERE_fit(p3cotPS, p2):                                   # all primed
    result = {}
    popt = curve_fit_try(ERE_0, p2, p3cotPS,1)
    result["a1_0"] = popt[0]
    popt = curve_fit_try(ERE_1_lin, p2, p3cotPS,2)
    result["a1_1"] = 0 if popt[0] == 0 else -np.sign(popt[0])*np.power(1/abs(popt[0]),1/3.)
    result["r1_1"] = 0 if popt[1] == 0 else 1/(2*popt[1])
    popt = curve_fit_try(ERE_2_lin, p2, p3cotPS,3)
    result["a1_2"] = 0 if popt[0] == 0 else -np.sign(popt[0])*np.power(1/abs(popt[0]),1/3.)
    result["r1_2"] = 0 if popt[1] == 0 else 1/(2*popt[1])
    result["c1_2"]=popt[2]
    return result

####################### Ab hier wurden Änderungen gemacht. Bitte vorsichtig sein ###########################

def ECM_p2(p2):
    return 2*np.sqrt(1+p2)
def RES_Drach(p2, m_R, gVPP2):
    ECM = ECM_p2(p2)
    if gVPP2 == 0:
        return 0
    else:
        return 6*np.pi*(m_R**2-ECM**2)/gVPP2
def RES_Alex_BWI(p2, m_R, gVPP2):                 # ident zu "RES_Drach"
    ECM = ECM_p2(p2)
    return ECM*gVPP2*p2**(3/2)/(6*np.pi*ECM**2*(m_R**2-ECM**2))
def RES_Alex_BWII(p2, m_R, gVPP2, r0):                                      # possibly include first guess for mR2 or swap around
    ECM = ECM_p2(p2)
    m_R2 = m_R**2
    k_R2 = m_R2 - 1
    return ECM/(m_R2-ECM**2)*(gVPP2*ECM**2/6*np.pi)*p2**(3/2)*(1+k_R2*r0**2)/(1+p2*r0**2)



    # return ECM*gVPP2*p2**(3/2)*(1+np.sqrt((k_R/2)**2-1))/(6*np.pi*ECM**2*(k_R**2-ECM**2)*(1+p2*r0**2))
def RES_fit(p3cotPS_ECM, tan_PS, p2):                                   # all primed
    result = {}
    popt = curve_fit_try(RES_Drach, p2, p3cotPS_ECM,2)
    result["m_R_D"], result["gVPP2_D"] = popt
    popt = curve_fit_try(RES_Alex_BWI, p2, tan_PS,2)
    result["m_R_BWI"], result["gVPP2_BWI"] = popt
    # popt = curve_fit_try(RES_Alex_BWII, p2, tan_PS,3)                   # maybe wrong. check later
    # result["m_R_BWII"], result["gVPP2_BWII"], result["r0_BWII"] = popt
    return result

def get_fits(res, res_spl):
    res_tmp={}
    res_spl_tmp={}
    for key, val in ERE_fit(res["p3cotPS_prime"],res["p2star_prime"]).items():
        res_tmp[key] = val
    for key, val in RES_fit(res["p3cotPS_Ecm_prime"],res["tan_PS"],res["p2star_prime"]).items():
        res_tmp[key] = val
        
    for key in res_tmp.keys():
        res_spl_tmp[key] = []

    for i in range(len(res_spl["p3cotPS_prime"])):
        for key, val in ERE_fit(res_spl["p3cotPS_prime"][i],res_spl["p2star_prime"][i]).items():
            res_spl_tmp[key].append(val)
        for key, val in RES_fit(res_spl["p3cotPS_Ecm_prime"][i],res_spl["tan_PS"][i],res_spl["p2star_prime"][i]).items():
            res_spl_tmp[key].append(val)

    # print(res_tmp)
    return res_tmp, res_spl_tmp

####################### Bis hier wurden Änderungen gemacht. Bitte vorsichtig sein ###########################

def only(arr):
    if len(arr) == 1:
        return arr[0]
    raise ValueError(f"Expected exactly one element, got {len(arr)}")

def genfromtxt_skip_empty(filename, **kwargs):
    with open(filename) as f:
        lines = [line for line in f if line.strip()]
        return np.genfromtxt(lines, **kwargs)

def fit_one_phaseshift(h5file_in, h5file_out, input_file, beta, m0):
    res_scat = {}
    res_spl_scat = {}

    infile = np.transpose(genfromtxt_skip_empty(input_file,delimiter=";",skip_header=1,dtype=str))

    with h5py.File(h5file_out,"a") as hfile:
        # with h5py.File(h5file_in,"r") as hfilein:
        for ens in hfile:
            if str(beta) in ens and str(m0) in ens:
                # ens_here = ens
                for P in hfile[ens]:
                    if P[0] == "p":
                        for irrep in hfile[ens][P]:
                            for lv in hfile[ens][P][irrep]:
                                if lv[:2] == "lv":
                                    fit_in = infile[1][infile[0] == ens+P+irrep+lv]
                                    fit = None
                                    if len(fit_in) == 0:
                                        fit = "False"
                                    else:
                                        fit = fit_in[0]
                                    if fit == "True":
                                        hfile[ens][P][irrep][lv]["fit"] = True
                                        p2star_prime = hfile[ens][P][irrep][lv]["mean"]["p2star_prime"][()]
                                        if 0 < p2star_prime < 15:
                                            if res_scat == {}:
                                                for key in hfile[ens][P][irrep][lv]["mean"]:
                                                    res_scat[key] = []
                                                    res_spl_scat[key] = []
                                            for tmp in hfile[ens][P][irrep][lv]["mean"]:
                                                if type(hfile[ens][P][irrep][lv]["mean"][tmp][()]) == np.complex128:                                                              
                                                    res_scat[tmp].append(hfile[ens][P][irrep][lv]["mean"][tmp][()].real)
                                                    res_spl_scat[tmp].append(hfile[ens][P][irrep][lv]["sample"][tmp][()].real)
                                                else:
                                                    res_scat[tmp].append(hfile[ens][P][irrep][lv]["mean"][tmp][()])
                                                    res_spl_scat[tmp].append(hfile[ens][P][irrep][lv]["sample"][tmp][()])
                                        else:
                                            raise RuntimeError("Energy not in elastic window at: %s"%(ens+P+irrep+lv))
                                    elif fit == "False":
                                        hfile[ens][P][irrep][lv]["fit"] = False
                                        pass
                                    else:
                                        raise RuntimeError("Wrong assignment in 'fit_scatter_input.csv' at: %s"%(ens+P+irrep+lv))

        
        if res_scat == {} or len(res_scat["p2star_prime"]) < 3:
            raise RuntimeError("Less than 3 energy levels selected in 'fit_scatter_input.csv' for beta=%f and m0%f"%(beta, m0))

        for key in res_scat:
            res_scat[key] = np.asarray(res_scat[key])
            res_spl_scat[key] = np.transpose(np.asarray(res_spl_scat[key]))

        
        fit_beta_m = "fit_b%f_m%f"%(beta,m0)

        res_fit, res_spl_fit = get_fits(res_scat,res_spl_scat)
        for key, val in res_fit.items():
            mean_group = hfile.require_group(fit_beta_m+"/mean")
            mean_group.create_dataset(key, data=val)
        for key, val in res_spl_fit.items():
            spl_group = hfile.require_group(fit_beta_m+"/sample")
            spl_group.create_dataset(key, data=val)
        # for key, val in res_fit.items():
        #     hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"mean/"+key, data = val)
        # for key, val in res_spl_fit.items():
        #     hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"sample/"+key, data = val)
        # for key, val in res_scat.items():
        #     hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"mean/"+key, data = val)
        # for key, val in res_spl_scat.items():
        #     hfile.create_dataset("fit_scatter_b%f_m%f/"%(beta,m0)+"sample/"+key, data = val)

def fit_all_phase_shifts(h5file_in, h5file_out, input_file):
    print("Fitting phase shifts...")
    fit_one_phaseshift(h5file_in, h5file_out, input_file,6.9,-0.92)
    fit_one_phaseshift(h5file_in, h5file_out, input_file,7.05,-0.863)
    fit_one_phaseshift(h5file_in, h5file_out, input_file,7.05,-0.867)
    print("Done!")


if __name__ == "__main__":

    args = sys.argv
    h5file_in = args[1]
    h5file_out = args[2]
    input_file = args[3]

    fit_all_phase_shifts(h5file_in, h5file_out,input_file)