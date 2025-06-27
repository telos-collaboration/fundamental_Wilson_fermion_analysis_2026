import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm
from tqdm import tqdm
import h5py
import sys
import os

# def save_to_hdf_old(res,res_sample, info, filename):
#     directory = "../output/scattering/"
#     fpathname = directory+filename+".hdf5"
#     if not os.path.exists(directory):
#         os.makedirs(directory)
#     with h5py.File(fpathname,"w") as hfile:
#     # with h5py.File("output/hdf5/"+filename+".hdf5","w") as hfile:
#         for key, val in res.items():
#             hfile.create_dataset("mean/"+key, data = val)
#         for key, val in res_sample.items():
#             hfile.create_dataset("sample/"+key, data = val)
#         for key, val in info.items():
#             hfile.create_dataset("info/"+key, data = val)

def save_to_hdf(res,res_sample, info, ens, P, irrep, lv, filename):
    group = ens+"/"+P+"/"+irrep+"/"+"lv"+"%i/"%lv
    # directory = "../output/scattering/"
    # fpathname = directory+filename+".hdf5"
    # if not os.path.exists(directory):
    #     os.makedirs(directory)
    # print("../output/scattering/isospin1_scattering"+filename+".hdf5")
    with h5py.File("../output/scattering/isospin1_scattering"+filename+".hdf5","a") as hfile:
        for key, val in res.items():
            hfile.create_dataset(group+"mean/"+key, data = val)
        for key, val in res_sample.items():
            hfile.create_dataset(group+"sample/"+key, data = val)
        for key, val in info.items():
            hfile.create_dataset(group+"info/"+key, data = val)

# def read_from_hdf(filename):                  # not needed for now
#     res, res_tmp, info = [{},{},{}]
#     with h5py.File("../output/scattering/"+filename+".hdf5","r") as hfile:
#         for key in hfile.keys():
#             if key[:5] == "mean/":
#                 res[key[5:]] = hfile[key][()]
#             if key[:7] == "sample/":
#                 res_tmp[key[7:]] = hfile[key][()]
#             if key[:5] == "info/":
#                 info[key[5:]] = hfile[key][()]
#     return res, res_tmp, info

def result_sampled(N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,mpi,irrep,ld,resampling="gauss",num_resample=50):
    res = {}
    info = {}
    res_sample = {}
    info["resampling"] = resampling
    info["num_resample"] = num_resample
    # info["mpi"] = mpi
    # res["N_L"] = N_L
    # print(type(N_L), type(mpi))
    # print(N_L, mpi)
    # print(N_L,mpi)
    info["L_prime"] = N_L*mpi
    info["dvec"] = dvec
    info["d"] = [np.sqrt(np.dot(dvec[i],dvec[i])) for i in range(len(dvec))]
    info["d2"] = [np.dot(dvec[i],dvec[i]) for i in range(len(dvec))]

    # for key, val in info.items():
    #     res[key] = val
    for key, val in get_rizz(E_pipi,N_L,dvec,mpi,irrep,ld).items():
        res[key] = val
    for key in res.keys():
        res_sample[key] = []

    # print("res:  ", res.keys())
    # print(res_sample.keys())
    # exit()

    if resampling == "gauss":
        for i in tqdm(range(num_resample)):
            if i < num_resample//2:
                E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
            else:
                E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep,ld)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    elif resampling == "lin":
        for E_pipi_tmp in tqdm(np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_resample)):
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep,ld)
            # print("res_tmp:  ", res_tmp.keys())
            for key, val in res_tmp.items():
                # print(key)
                res_sample[key].append(val)
    return res, res_sample, info


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

# def Ecm_prime(E_prime, P_prime):
#     return np.sqrt(E_prime**2-P_prime**2)

def Ecm(ld):
    if ld:
        return Ecm_ld
    else:
        return Ecm_cont

def Ecm_cont(E, Pvec):
    return np.sqrt(E**2-np.dot(Pvec,Pvec))

def Ecm_ld(E, Pvec):
    return np.arccosh(np.cosh(E)-2*(np.sin(Pvec[0]/2)**2+np.sin(Pvec[1]/2)**2+np.sin(Pvec[2]/2)**2))

# def pstar(E, P):
#     return np.sqrt(E**2-P**2)
def pstar(ld):
    if ld:
        return pstar_ld
    else:
        return pstar_cont

def pstar_cont(Ecm, mpi):
    return np.sqrt(Ecm**2/4-mpi**2)

def pstar_ld(Ecm, mpi):
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(Ecm/2)-np.cosh(mpi))))

# def pstar_prime(Ecm_prime):
#     return np.sqrt(Ecm_prime**2/4-1)

# def Ecm_lat_disp(E, Pvec):
#     return np.arccosh(np.cosh(E)-2*(np.sin(Pvec[0]/2)**2+np.sin(Pvec[1]/2)**2+np.sin(Pvec[2]/2)**2))

# def pstar_lat_disp(Ecm, mpi):
#     return 2*np.arcsin(np.sqrt(0.5*(np.cosh(Ecm/2)-np.cosh(mpi))))

def get_rizz(E_pipi, N_L, dvec, mpi, irrep,ld):
    res = {}
    # key_list = ["En","En_prime","E_cm","E_cm_prime","E_cm_ld","E_cm_ld_prime","s","s_prime","s_ld","s_ld_prime","pstar","pstar_prime","pstar_ld","pstar_ld_prime","p2star","p2star_prime","p2star_ld","p2star_ld_prime","q","q_ld","q2","q2_ld","cot_PS","cot_PS_ld","tan_PS","tan_PS_ld","PS","PS_ld", "p3cotPS", "p3cotPS_prime", "p3cotPS_ld", "p3cotPS_ld_prime", "p3cotPS_Ecm", "p3cotPS_Ecm_prime", "p3cotPS_Ecm_ld", "p3cotPS_Ecm_ld_prime", "sigma", "sigma_prime", "sigma_ld", "sigma_ld_prime"]
    key_list = ["En_prime","E_cm","E_cm_prime","s","s_prime","pstar","pstar_prime","p2star","p2star_prime","q","q2","cot_PS","tan_PS","PS", "p3cotPS", "p3cotPS_prime", "p3cotPS_Ecm", "p3cotPS_Ecm_prime", "sigma", "sigma_prime"]

    # for key in key_list:
        # result[key] = []

    # for i in range(len(E_pipis)):
    # res["En"].append(E_pipis[i])

    # for x in [dvec,N_L]:
    #     print(type(x),x)

    Pvec = [2*np.pi*x/N_L for x in dvec]
    # for x in [Pvec,]:
    #     print(type(x),x)
    P_prime = 2*np.pi*np.sqrt(np.dot(dvec,dvec))/(N_L*mpi)
    En_prime = E_pipi/mpi
    res["En_prime"] = En_prime
    # print("\t\t", En_prime**2, P_prime**2, En_prime**2 - P_prime**2)
    if En_prime**2 - P_prime**2 < 4:
        for key in key_list:
            res[key] = 0
    else:
        # res = {}
        res["E_cm"] = Ecm(ld)(E_pipi,Pvec)
        res["E_cm_prime"] = res["E_cm"]/mpi
        res["s_prime"] = res["E_cm_prime"]**2
        res["s"] = res["E_cm"]**2
        res["pstar"] = pstar(ld)(res["E_cm"],mpi)    # np.sqrt(res["s_prime"]/4-1)
        res["pstar_prime"] = res["pstar"]/mpi
        res["p2star"] = res["pstar"]**2
        res["p2star_prime"] = res["pstar_prime"]**2
        q2 = res["p2star"]*(N_L/(2*np.pi))**2
        res["q"] = np.sqrt(q2)
        res["q2"] = q2
        cot_PS = cot_delta_mom(dvec,irrep)(q2, N_L,mpi)
        res["cot_PS"] = cot_PS
        res["tan_PS"] = 1/cot_PS
        PS = 360*np.arctan(1/cot_PS)/(2*np.pi)
        res["PS"] = complex(PS.real%180,PS.imag%180)
        res["p3cotPS"] = res["pstar"]**3*cot_PS
        res["p3cotPS_prime"] = res["pstar_prime"]**3*cot_PS
        res["p3cotPS_Ecm"] = res["pstar"]**3/res["E_cm"]*cot_PS
        res["p3cotPS_Ecm_prime"] = res["pstar_prime"]**3/res["E_cm_prime"]*cot_PS
        res["sigma"] = 4*np.pi*3/(cot_PS**2+1)/res["p2star"]
        res["sigma_prime"] = 4*np.pi*3/(cot_PS**2+1)/res["p2star_prime"]

            # tmp["E_cm_ld"] = Ecm_lat_disp(E_pipis[i],Pvec)
            # tmp["E_cm_ld_prime"] = tmp["E_cm_ld"]/mpi
            # tmp["s_ld"] = tmp["E_cm_ld"]**2
            # tmp["s_ld_prime"] = tmp["E_cm_ld_prime"]**2
            # tmp["pstar_ld"] = pstar_lat_disp(tmp["E_cm_ld"],mpi).real
            # tmp["pstar_ld_prime"] = tmp["pstar_ld"]/mpi
            # tmp["p2star_ld"] = tmp["pstar_ld"]**2
            # tmp["p2star_ld_prime"] = tmp["pstar_ld_prime"]**2
            # q2_ld = tmp["p2star_ld"]*(N_Ls[i]/(2*np.pi))**2
            # tmp["q2_ld"] = q2_ld
            # tmp["q_ld"] = np.sqrt(q2_ld)
            # cot_PS_ld = cot_delta_mom(dvecs[i],irrep)(q2_ld, N_Ls[i],mpi).real
            # tmp["cot_PS_ld"] = cot_PS_ld
            # tmp["tan_PS_ld"] = 1/cot_PS_ld
            # PS_ld = 360*np.arctan(1/cot_PS_ld)/(2*np.pi)
            # tmp["PS_ld"] = complex(PS_ld.real%180,PS_ld.imag%180)
            # tmp["p3cotPS_ld"] = tmp["pstar_ld"]**3*cot_PS_ld
            # tmp["p3cotPS_ld_prime"] = tmp["pstar_ld_prime"]**3*cot_PS_ld
            # tmp["p3cotPS_Ecm_ld"] = tmp["pstar_ld"]**3*cot_PS_ld/tmp["E_cm_ld"]
            # tmp["p3cotPS_Ecm_ld_prime"] = tmp["pstar_ld_prime"]**3*cot_PS_ld/tmp["E_cm_ld_prime"]
            # tmp["sigma_ld"] = 12*np.pi/(cot_PS_ld**2+1)/tmp["p2star_ld"]
            # tmp["sigma_ld_prime"] = 12*np.pi/(cot_PS_ld**2+1)/tmp["p2star_ld_prime"]


            # for key, val in tmp.items():
            #     result[key].append(tmp[key])
    return res

# def read_hdf5_fitresults_old(pref, NTs, NLs, beta, m0, num_lvl=2):
#     NL_arr, dvec_arr, en_arr, en_m_arr, en_p_arr = [[],[],[],[],[]]
#     with h5py.File("output/data/isospin1_fitresults"+pref+".hdf5","r") as hfile:
#         for i in range(len(NTs)):
#             NT=NTs[i]
#             for NL in NLs[i]:
#                 hfile_str = "Lt%iLs%ibeta"%(NT, NL)+("%f"%beta).rstrip("0")+"m"+("%f"%m0).rstrip("0")+"/"
#                 for key in hfile[hfile_str].keys():
#                     for i in range(num_lvl):
#                         NL_arr.append(float(NL))
#                         dvec_arr.append([int(key[2]),int(key[4]),int(key[6])])
#                         en_arr.append(hfile[hfile_str+key+"/E%i"%i][()][0])
#                         en_m_arr.append(hfile[hfile_str+key+"/Delta_E%i"%i][()][0])
#                         en_p_arr.append(hfile[hfile_str+key+"/Delta_E%i"%i][()][0])
#     return np.asarray(NL_arr), np.asarray(dvec_arr), np.asarray(en_arr), np.asarray(en_m_arr), np.asarray(en_p_arr)

    #     for i in range(len(NTs)):
    #         NT=NTs[i]
    #         for NL in NLs[i]:
    #             hfile_str = "Lt%iLs%ibeta"%(NT, NL)+("%f"%beta).rstrip("0")+"m"+("%f"%m0).rstrip("0")+"/"
    #             for key in hfile[hfile_str].keys():
    #                 for i in range(num_lvl):
    #                     NL_arr.append(float(NL))
    #                     dvec_arr.append([int(key[2]),int(key[4]),int(key[6])])
    #                     en_arr.append(hfile[hfile_str+key+"/E%i"%i][()][0])
    #                     en_m_arr.append(hfile[hfile_str+key+"/Delta_E%i"%i][()][0])
    #                     en_p_arr.append(hfile[hfile_str+key+"/Delta_E%i"%i][()][0])
    # return np.asarray(NL_arr), np.asarray(dvec_arr), np.asarray(en_arr), np.asarray(en_m_arr), np.asarray(en_p_arr)
    
# def mpi_m0(m0):
#     if m0 == -0.92:
#         return 0.38649
#     elif m0 == -0.863:
#         return 0.20590
#     elif m0 == -0.867:
#         return 0.14810

# def mrho_m0(m0):
#     if m0 == -0.92:
#         return 0.5494
#     elif m0 == -0.863:
#         return 0.3773
#     elif m0 == -0.867:
#         return 0.3530

def calc_all_phaseshifts(corrfitname,pref = "std",resampling="lin",num_resample=5,num_lv=2):
    info = {}
    infile = np.transpose(np.genfromtxt("../input/scattering_input.csv",delimiter=";",skip_header=1,dtype=str))
    # infile[3] = [float(infile[3,i]) for i in range(len(infile[0]))]
    # infile[4] = [ bool(infile[4,i]) for i in range(len(infile[0]))]
    # exit()
    with h5py.File("../output/data/isospin1_fitresults"+corrfitname+".hdf5","r") as hfile:
        for ens in hfile:
            # print(ens)
            for P in hfile[ens]:
                if P[0] == "p":
                    # print(P)
                    dvec = [int(P[2]),int(P[4]),int(P[6])]
                    for irrep in hfile[ens][P]:
                        beta, m0, mpi, mrho, ld = infile[1:,infile[0] == ens+P+irrep]
                        beta = float(beta)
                        m0 = float(m0)
                        mpi = float(mpi)
                        mrho = float(mrho)
                        ld = bool(ld)
                        info["beta"] = beta
                        info["m0"] = m0
                        info["mpi"] = mpi
                        info["mrho"] = mrho
                        # for x in [beta, m0, mpi, ld]:
                        #     print(type(x),x)
                        NL = hfile[ens]["lattice"][()][3]
                        # print(beta,m0,NL,dvec)
                        for i in range(num_lv):
                            # print("E=",i)
                            E = hfile[ens][P][irrep]["E%i"%i][()][0]
                            E_m = hfile[ens][P][irrep]["Delta_E%i"%i][()][0]
                            E_p = hfile[ens][P][irrep]["Delta_E%i"%i][()][0]
                            res, res_sampled, info_tmp = result_sampled(NL, E, E_m, E_p, dvec, mpi, irrep, ld, resampling=resampling, num_resample=num_resample)
                            for key, val in info_tmp.items():
                                info[key] = val
                            save_to_hdf(res, res_sampled, info, ens, P, irrep, i, corrfitname)
                            # print(res)
                            # for key, val in res.items():
                            #     print(key, val)
                        # E0 = hfile[ens][P][irrep]["E0"][()][0]
                        # E0_m = hfile[ens][P][irrep]["Delta_E0"][()][0]
                        # E0_p = hfile[ens][P][irrep]["Delta_E0"][()][0]
                        # E1 = hfile[ens][P][irrep]["E1"][()][0]
                        # E1_m = hfile[ens][P][irrep]["Delta_E1"][()][0]
                        # E1_p = hfile[ens][P][irrep]["Delta_E1"][()][0]
                        # res, res_sampled = result_sampled(NL, E1, E1_m, E1_p, dvec, mpi, irrep, ld, resampling=resampling, num_resample=num_resample)
                        # print(NL, E0, E0_m, E0_p, E1, E1_m, E1_p, dvec)



    # NL_arr, dvec_arr, en_arr, en_m_arr, en_p_arr = 
    # read_hdf5_fitresults(corrfitname,num_lvl=2)
    # mpi = mpi_m0(m0)
    # mrho = mrho_m0(m0)
    # info={}
    # info["beta"],info["m_1"],info["m_2"], info["mrho"], info["mpi"], info["en_lv"] = [beta,m0,m0,mrho,mpi,2]

if __name__ == "__main__":
    calc_all_phaseshifts("_evp_deriv_false", "non_res", resampling="lin", num_resample=3)