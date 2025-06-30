import numpy as np
from tqdm import tqdm
import h5py
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def save_to_hdf(res,res_sample, filename):
    with h5py.File("rho_pipi_scattering_analysis/output/hdf5/"+filename+".hdf5","w") as hfile:
    # with h5py.File("output/hdf5/"+filename+".hdf5","w") as hfile:
        for key, val in res.items():
            hfile.create_dataset("orig_"+key, data = val)
        for key, val in res_sample.items():
            hfile.create_dataset("sample_"+key, data = val)

def read_from_hdf(filename):
    res, res_spl = [{},{}]
    with h5py.File("rho_pipi_scattering_analysis/output/hdf5/"+filename+".hdf5","r") as hfile:
        for key in hfile.keys():
            if key[:5] == "orig_":
                res[key[5:]] = hfile[key][()]
            if key[:7] == "sample_":
                res_spl[key[7:]] = hfile[key][()]
    return res, res_spl

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
    popt, pcov = curve_fit(ERE_0, p2, p3cotPS)
    result["a1_0"] = popt[0]
    popt, pcov = curve_fit(ERE_1_lin, p2, p3cotPS)
    result["a1_1"]=-np.sign(popt[0])*np.power(1/abs(popt[0]),1/3.)
    result["r1_1"]=1/(2*popt[1])
    popt, pcov = curve_fit(ERE_2_lin, p2, p3cotPS)
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
    # popt, pcov = curve_fit(RES_Drach, ECM, p3cotPS_ECM)
    # result["m_R_D"], result["gVPP2_D"] = popt
    # popt, pcov = curve_fit(RES_Alex_BWI, p, tan_PS)
    # result["m_R_BWI"], result["gVPP2_BWI"] = popt
    # popt, pcov = curve_fit(RES_Alex_BWII, p, tan_PS)
    # result["m_R_BWII"], result["gVPP2_BWII"], result["r0_BWII"] = popt
    return result

def w_ind_m0(m0):
    return [1,2,3,4,5,6,7,8,9,10,11,12,13]
    # if m0 == 0.92:
    #     return [4, 10, 13, 16]
    # elif m0 == -0.863:
    #     return [4,7,10,13,16]


def get_fits(res, res_spl, ld = "_ld",m0 = -0.92):
    res_tmp={}
    res_spl_tmp={}
    # win_ind=[]
    # for i in range(len(res["E_cm"+ld+"_prime"])):
    #     if res["E_cm"+ld+"_prime"][i] > 2 and res["E_cm"+ld+"_prime"][i] < 4:
    #         win_ind.append(i)

    # win_ind = [4, 10, 13, 16]
    win_ind = w_ind_m0(m0)

    for key, val in ERE_fit(res["p3cotPS"+ld+"_prime"][win_ind],res["p2star"+ld+"_prime"][win_ind]).items():
        res_tmp[key] = val
    for key, val in RES_fit(res["p3cotPS_Ecm"+ld+"_prime"][win_ind],res["E_cm"+ld+"_prime"][win_ind],res["tan_PS"+ld][win_ind],res["pstar"+ld+"_prime"][win_ind]).items():
        res_tmp[key] = val
        
    for key in res_tmp.keys():
        res_spl_tmp[key] = []

    for i in range(len(res_spl["p3cotPS"+ld+"_prime"])):
        for key, val in ERE_fit(res_spl["p3cotPS"+ld+"_prime"][i][win_ind],res_spl["p2star"+ld+"_prime"][i][win_ind]).items():
            res_spl_tmp[key].append(val)
        for key, val in RES_fit(res_spl["p3cotPS_Ecm"+ld+"_prime"][i][win_ind],res_spl["E_cm"+ld+"_prime"][i][win_ind],res_spl["tan_PS"+ld][i][win_ind],res_spl["pstar"+ld+"_prime"][i][win_ind]).items():
            res_spl_tmp[key].append(val)
    
    for key, val in res_tmp.items():
        res[key] = val
    for i in range(len(res_spl["p3cotPS"+ld+"_prime"])):
        for key, val in res_spl_tmp.items():
            res_spl[key] = val

    return res, res_spl

def fit_one_phaseshift(name, beta, m0):
    res_calc = {}
    res_spl_calc = {}

    with h5py.File("../output/scattering/isospin1_fit_scatter"+name+".hdf5","r") as hfile:
        for ens in hfile:
            print(ens)
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
                                        
                                        
                                        
                                        
                                            # print(res_calc)
    # for key, val in res_calc.items():
    #     print(key)
    #     print(val)
    # for key, val in res_spl_calc.items():
    #     print(key)
    #     # print(val)
    #     print(len(val))
    #     for tmp in val:
    #         print(len(tmp))
    # print(res_spl_calc)
    
    # print(num_resample)




                                # asd = hfile[ens][P][irrep][lv]["mean"][()]
                                # print(type(hfile[ens][P][irrep][lv]["mean"][()]))
                                # exit()
                                # res_calc.append()

                            # beta, m0, mpi, mrho, ld = infile[1:,infile[0] == ens+P+irrep]
            #             beta = float(beta)
            #             m0 = float(m0)
            #             mpi = float(mpi)
            #             mrho = float(mrho)
            #             ld = bool(ld)
            #             info["beta"] = beta
            #             info["m0"] = m0
            #             info["mpi"] = mpi
            #             info["mrho"] = mrho
            #             NL = hfile[ens]["lattice"][()][3]
            #             info["NL"] = NL
            #             # print(beta,m0,NL,dvec)
            #             for i in range(num_lv):
            #                 # print("E=",i)
            #                 E = hfile[ens][P][irrep]["E%i"%i][()][0]
            #                 E_m = hfile[ens][P][irrep]["Delta_E%i"%i][()][0]
            #                 E_p = hfile[ens][P][irrep]["Delta_E%i"%i][()][0]
            #                 res, res_sampled, info_tmp = result_sampled(NL, E, E_m, E_p, dvec, mpi, irrep, ld, resampling=resampling, num_resample=num_resample)
            #                 for key, val in info_tmp.items():
            #                     info[key] = val
            #                 save_to_hdf(res, res_sampled, info, ens, P, irrep, i, corrfitname)


if __name__ == "__main__":
    fit_one_phaseshift("_evp_deriv_false",6.9,-0.92)


    # filename = "close_res_new"
    # m0 = -0.863

    # res, res_spl = read_from_hdf(filename)
    # res, res_spl = get_fits(res, res_spl,m0=m0)
    # save_to_hdf(res, res_spl, filename+"_fit")