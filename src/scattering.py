import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm
from tqdm import tqdm
import h5py
import sys
import os

def save_to_hdf(res,res_sample, filename):
    with h5py.File("rho_pipi_scattering_analysis/output/hdf5/"+filename+".hdf5","w") as hfile:
    # with h5py.File("output/hdf5/"+filename+".hdf5","w") as hfile:
        for key, val in res.items():
            hfile.create_dataset("orig_"+key, data = val)
        for key, val in res_sample.items():
            hfile.create_dataset("sample_"+key, data = val)

def read_from_hdf(filename):
    res, res_tmp = [{},{}]
    with h5py.File("output/hdf5/"+filename+".hdf5","r") as hfile:
        for key in hfile.keys():
            if key[:5] == "orig_":
                res[key[5:]] = hfile[key][()]
            if key[:7] == "sample_":
                res_tmp[key[7:]] = hfile[key][()]
    return res, res_tmp

def result_sampled(info,N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,mpi,irrep,num_resample=50, resampling = "gauss"):
    res = {}
    res_sample = {}
    res["resampling"] = resampling
    res["num_resample"] = num_resample
    res["mpi"] = mpi
    res["N_L"] = N_L
    # print(type(N_L), type(mpi))
    # print(N_L, mpi)
    res["L_prime"] = N_L*mpi
    res["dvec"] = dvec
    res["d"] = [np.sqrt(np.dot(dvec[i],dvec[i])) for i in range(len(dvec))]
    res["d2"] = [np.dot(dvec[i],dvec[i]) for i in range(len(dvec))]

    for key, val in info.items():
        res[key] = val
    for key, val in get_rizz(E_pipi,N_L,dvec,mpi,irrep).items():
        res[key] = val
    for key in res.keys():
        res_sample[key] = []

    if resampling == "gauss":
        for i in tqdm(range(num_resample)):
            if i < num_resample//2:
                E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
            else:
                E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    elif resampling == "lin":
        for E_pipi_tmp in tqdm(np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_resample)):
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    return res, res_sample


def cot_delta_000_T1(q2,N_L,mpi):
    return wlm(0,0,0,0,0,mpi,mpi,q2,int(N_L))

def cot_delta_001_A1(q2,N_L,mpi):                                      # A1 (A2 in Luka)
    p = [0,0,1]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = 2*wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first + second
def cot_delta_001_E(q2,N_L,mpi):                                      # E (A2 in Luka)
    p = [0,0,1]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first - second

def cot_delta_002_A1(q2,N_L,mpi):
    p = [0,0,2]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = 2*wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first + second
def cot_delta_002_E(q2,N_L,mpi):
    p = [0,0,2]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first - second

def cot_delta_110_A1(q2,N_L,mpi):                                                          # A1=B2 (B3 in Luka)
    p = [1,1,0]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    third = complex(0,1) * np.sqrt(6) * wlm(2,2,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first - second - third
def cot_delta_110_B1(q2,N_L,mpi):                                                          # B1 (B1 or B2 in Luka? Ambiguity in sign)
    p = [1,1,0]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = 2*wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first + second # sign of third term depends on convention. Please check!!
# def cot_delta_110_B2(q2,N_L,mpi):                                                          # A1=B2 (B3 in Luka) # skipped for now
#     p = [1,1,0]
#     first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
#     second = wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
#     third = complex(0,1) * np.sqrt(6) * wlm(2,2,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
#     return first - second + third

def cot_delta_111_A1(q2,N_L,mpi):                                                              # A1 (A2 in Luka) (formula in 1206.4141v2 different [old]) 
    p = [1,1,1]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = complex(0,1)*np.sqrt(8/3)*wlm(2,2,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    term = wlm(2,1,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    third = np.sqrt(8/3)*(term.real+term.imag)
    return first - second - third
def cot_delta_111_E(q2,N_L,mpi):                                                              # E (E in Luka) 
    p = [1,1,1]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = complex(0,1)*np.sqrt(6)*wlm(2,2,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first + second


def cot_delta_mom(dvec, irrep):
    if list(dvec) == [0,0,0]:
        if irrep == "T1":
            return cot_delta_000_T1
    elif list(dvec) == [0,0,1] or  list(dvec) == [0,0,2]:
        if irrep == "A1":
            return cot_delta_001_A1
        elif irrep == "E":
            return cot_delta_001_E
    elif list(dvec) == [1,1,0]:
        if irrep == "A1":
            return cot_delta_110_A1
        elif irrep == "B1":
            return cot_delta_110_B1
    elif list(dvec) == [0,1,1]:
        if irrep == "A1":
            return cot_delta_110_A1
        elif irrep == "E":
            return cot_delta_110_B1
    elif list(dvec) == [1,1,1]:
        if irrep == "A1":
            return cot_delta_111_A1
        elif irrep == "B1":
            return cot_delta_111_E
    else:
        print("wrong momentum or irrep")
        exit()

def Ecm_prime(E_prime, P_prime):
    return np.sqrt(E_prime**2-P_prime**2)

def pstar_prime(Ecm_prime):
    return np.sqrt(Ecm_prime**2/4-1)

def Ecm_lat_disp(E, Pvec):
    return np.arccosh(np.cosh(E)-2*(np.sin(Pvec[0]/2)**2+np.sin(Pvec[1]/2)**2+np.sin(Pvec[2]/2)**2))

def pstar_lat_disp(Ecm, mpi):
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(Ecm/2)-np.cosh(mpi))))

def get_rizz(E_pipis, N_Ls, dvecs, mpi, irrep):
    result = {}
    key_list = ["En","En_prime","E_cm","E_cm_prime","E_cm_ld","E_cm_ld_prime","s","s_prime","s_ld","s_ld_prime","pstar","pstar_prime","pstar_ld","pstar_ld_prime","p2star","p2star_prime","p2star_ld","p2star_ld_prime","q","q_ld","q2","q2_ld","cot_PS","cot_PS_ld","tan_PS","tan_PS_ld","PS","PS_ld", "p3cotPS", "p3cotPS_prime", "p3cotPS_ld", "p3cotPS_ld_prime", "p3cotPS_Ecm", "p3cotPS_Ecm_prime", "p3cotPS_Ecm_ld", "p3cotPS_Ecm_ld_prime", "sigma", "sigma_prime", "sigma_ld", "sigma_ld_prime"]

    for key in key_list:
        result[key] = []

    for i in range(len(E_pipis)):
        result["En"].append(E_pipis[i])
        Pvec = 2*np.pi*dvecs[i]/N_Ls[i]
        P_prime = 2*np.pi*np.sqrt(np.dot(dvecs[i],dvecs[i]))/(N_Ls[i]*mpi)
        En_prime = E_pipis[i]/mpi
        result["En_prime"].append(En_prime)
        if En_prime**2 - P_prime**2 < 4:
            for key in key_list[2:]:
                result[key].append(0)
        else:
            tmp = {}
            tmp["E_cm_prime"] = Ecm_prime(E_pipis[i]/mpi,P_prime)
            tmp["E_cm"] = tmp["E_cm_prime"]*mpi
            tmp["s_prime"] = tmp["E_cm_prime"]**2
            tmp["s"] = tmp["E_cm"]**2
            tmp["pstar_prime"] = np.sqrt(tmp["s_prime"]/4-1)
            tmp["p2star_prime"] = tmp["pstar_prime"]**2
            tmp["pstar"] = tmp["pstar_prime"]*mpi
            tmp["p2star"] = tmp["p2star_prime"]*mpi**2
            q2 = tmp["p2star"]*(N_Ls[i]/(2*np.pi))**2
            tmp["q"] = np.sqrt(q2)
            tmp["q2"] = q2
            cot_PS = cot_delta_mom(dvecs[i],irrep)(q2, N_Ls[i],mpi)
            tmp["cot_PS"] = cot_PS
            tmp["tan_PS"] = 1/cot_PS
            PS = 360*np.arctan(1/cot_PS)/(2*np.pi)
            tmp["PS"] = complex(PS.real%180,PS.imag%180)
            tmp["p3cotPS"] = tmp["pstar"]**3*cot_PS
            tmp["p3cotPS_prime"] = tmp["pstar_prime"]**3*cot_PS
            tmp["p3cotPS_Ecm"] = tmp["pstar"]**3/tmp["E_cm"]*cot_PS
            tmp["p3cotPS_Ecm_prime"] = tmp["pstar_prime"]**3/tmp["E_cm_prime"]*cot_PS
            tmp["sigma"] = 4*np.pi*3/(cot_PS**2+1)/tmp["p2star"]
            tmp["sigma_prime"] = 4*np.pi*3/(cot_PS**2+1)/tmp["p2star_prime"]

            tmp["E_cm_ld"] = Ecm_lat_disp(E_pipis[i],Pvec)
            tmp["E_cm_ld_prime"] = tmp["E_cm_ld"]/mpi
            tmp["s_ld"] = tmp["E_cm_ld"]**2
            tmp["s_ld_prime"] = tmp["E_cm_ld_prime"]**2
            tmp["pstar_ld"] = pstar_lat_disp(tmp["E_cm_ld"],mpi).real
            tmp["pstar_ld_prime"] = tmp["pstar_ld"]/mpi
            tmp["p2star_ld"] = tmp["pstar_ld"]**2
            tmp["p2star_ld_prime"] = tmp["pstar_ld_prime"]**2
            q2_ld = tmp["p2star_ld"]*(N_Ls[i]/(2*np.pi))**2
            tmp["q2_ld"] = q2_ld
            tmp["q_ld"] = np.sqrt(q2_ld)
            cot_PS_ld = cot_delta_mom(dvecs[i],irrep)(q2_ld, N_Ls[i],mpi).real
            tmp["cot_PS_ld"] = cot_PS_ld
            tmp["tan_PS_ld"] = 1/cot_PS_ld
            PS_ld = 360*np.arctan(1/cot_PS_ld)/(2*np.pi)
            tmp["PS_ld"] = complex(PS_ld.real%180,PS_ld.imag%180)
            tmp["p3cotPS_ld"] = tmp["pstar_ld"]**3*cot_PS_ld
            tmp["p3cotPS_ld_prime"] = tmp["pstar_ld_prime"]**3*cot_PS_ld
            tmp["p3cotPS_Ecm_ld"] = tmp["pstar_ld"]**3*cot_PS_ld/tmp["E_cm_ld"]
            tmp["p3cotPS_Ecm_ld_prime"] = tmp["pstar_ld_prime"]**3*cot_PS_ld/tmp["E_cm_ld_prime"]
            tmp["sigma_ld"] = 12*np.pi/(cot_PS_ld**2+1)/tmp["p2star_ld"]
            tmp["sigma_ld_prime"] = 12*np.pi/(cot_PS_ld**2+1)/tmp["p2star_ld_prime"]


            for key, val in tmp.items():
                result[key].append(tmp[key])
    return result

def read_hdf5_fitresults(pref, NTs, NLs, beta, m0, num_lvl=2):
    NL_arr, dvec_arr, en_arr, en_m_arr, en_p_arr = [[],[],[],[],[]]
    with h5py.File("output/data/isospin1_fitresults"+pref+".hdf5","r") as hfile:
        for i in range(len(NTs)):
            NT=NTs[i]
            for NL in NLs[i]:
                hfile_str = "Lt%iLs%ibeta"%(NT, NL)+("%f"%beta).rstrip("0")+"m"+("%f"%m0).rstrip("0")+"/"
                for key in hfile[hfile_str].keys():
                    for i in range(num_lvl):
                        NL_arr.append(float(NL))
                        dvec_arr.append([int(key[2]),int(key[4]),int(key[6])])
                        en_arr.append(hfile[hfile_str+key+"/E%i"%i][()][0])
                        en_m_arr.append(hfile[hfile_str+key+"/Delta_E%i"%i][()][0])
                        en_p_arr.append(hfile[hfile_str+key+"/Delta_E%i"%i][()][0])
    return np.asarray(NL_arr), np.asarray(dvec_arr), np.asarray(en_arr), np.asarray(en_m_arr), np.asarray(en_p_arr)
    
def mpi_m0(m0):
    if m0 == -0.92:
        return 0.38649
    elif m0 == -0.863:
        return 0.20590
    elif m0 == -0.867:
        return 0.14810

def mrho_m0(m0):
    if m0 == -0.92:
        return 0.5494
    elif m0 == -0.863:
        return 0.3773
    elif m0 == -0.867:
        return 0.3530

def calc_PS_Fabian(fitresname,NTs,NLs,beta,m0,pref = "",resampling="lin",num_resample=5):
    NL_arr, dvec_arr, en_arr, en_m_arr, en_p_arr = read_hdf5_fitresults(fitresname,NTs,NLs,beta,m0,num_lvl=2)
    mpi = mpi_m0(m0)
    mrho = mrho_m0(m0)
    info={}
    info["beta"],info["m_1"],info["m_2"], info["mrho"], info["mpi"], info["en_lv"] = [beta,m0,m0,mrho,mpi,2]

if __name__ == "__main__":
    calc_PS_Fabian("fit_non_res", "non_res", resampling="lin", num_resample=100)