import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm
from tqdm import tqdm
import h5py


def save_to_hdf(res,res_sample, filename):
    with h5py.File("output/hdf5/"+filename+".hdf5","w") as hfile:
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

def result_sampled(info,N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,mpi,num_resample=50, resampling = "gauss"):
    res = {}
    res_sample = {}
    res["resampling"] = resampling
    res["num_resample"] = num_resample
    res["mpi"] = mpi
    res["N_L"] = N_L
    res["L_prime"] = N_L*mpi
    res["dvec"] = dvec
    res["d"] = [np.sqrt(np.dot(dvec[i],dvec[i])) for i in range(len(dvec))]
    res["d2"] = [np.dot(dvec[i],dvec[i]) for i in range(len(dvec))]

    for key, val in info.items():
        res[key] = val
    for key, val in get_rizz(E_pipi,N_L,dvec,mpi).items():
        res[key] = val
    for key in res.keys():
        res_sample[key] = []

    if resampling == "gauss":
        for i in tqdm(range(num_resample)):
            if i < num_resample//2:
                E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
            else:
                E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
    elif resampling == "lin":
        for E_pipi_tmp in tqdm(np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_resample)):
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    return res, res_sample


def cot_delta_000(q2,N_L):
    return wlm(0,0,0,0,0,1,1,q2,int(N_L))
def cot_delta_001(q2,N_L):
    return wlm(0,0,0,0,1,1,1,q2,int(N_L)) + 2*wlm(2,0,0,0,1,1,1,q2,int(N_L))
def cot_delta_110(q2,N_L):
    return wlm(0,0,1,1,0,1,1,q2,int(N_L)) - wlm(2,0,1,1,0,1,1,q2,int(N_L)) + np.sqrt(3/2) * complex(0,1) * (wlm(2,2,1,1,0,1,1,q2,int(N_L)) - wlm(2,-2,1,1,0,1,1,q2,int(N_L)))

def cot_delta_mom(dvec):
    if list(dvec) == [0,0,0]:
        return cot_delta_000
    if list(dvec) == [0,0,1]:
        return cot_delta_001
    elif list(dvec) == [1,1,0]:
        return cot_delta_110
    else:
        print("wrong momentum")
        exit()

def Ecm_prime(E_prime, P_prime):
    return np.sqrt(E_prime**2-P_prime**2)

def pstar_prime(Ecm_prime):
    return np.sqrt(Ecm_prime**2/4-1)

def Ecm_lat_disp(E, Pvec):
    return np.arccosh(np.cosh(E)-2*(np.sin(Pvec[0]/2)**2+np.sin(Pvec[1]/2)**2+np.sin(Pvec[2]/2)**2))

def pstar_lat_disp(Ecm, mpi):
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(Ecm/2)-np.cosh(mpi))))

def get_rizz(E_pipis, N_Ls, dvecs, mpi):
    result = {}
    key_list = ["En","En_prime","E_cm","E_cm_prime","E_cm_ld","E_cm_ld_prime","s","s_prime","s_ld","s_ld_prime","pstar","pstar_prime","pstar_ld","pstar_ld_prime","p2star","p2star_prime","p2star_ld","p2star_ld_prime","q","q_ld","q2","q2_ld","cot_PS","cot_PS_ld","tan_PS","tan_PS_ld","PS","PS_ld", "p3cotPS", "p3cotPS_prime", "p3cotPS_ld", "p3cotPS_ld_prime", "p3cotPS_Ecm", "p3cotPS_Ecm_prime", "p3cotPS_Ecm_ld", "p3cotPS_Ecm_ld_prime"]

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
            tmp["q2"] = q2
            cot_PS = cot_delta_mom(dvecs[i])(q2, N_Ls[i])
            tmp["cot_PS"] = cot_PS
            tmp["tan_PS"] = 1/cot_PS
            PS = 360*np.arctan(1/cot_PS)/(2*np.pi)
            tmp["PS"] = complex(PS.real%180,PS.imag%180)
            tmp["p3cotPS"] = tmp["pstar"]**3*cot_PS
            tmp["p3cotPS_prime"] = tmp["pstar_prime"]**3*cot_PS
            tmp["p3cotPS_Ecm"] = tmp["pstar"]**3/tmp["E_cm"]*cot_PS
            tmp["p3cotPS_Ecm_prime"] = tmp["pstar_prime"]**3/tmp["E_cm_prime"]*cot_PS

            tmp["E_cm_ld"] = Ecm_lat_disp(E_pipis[i],Pvec)
            tmp["E_cm_ld_prime"] = tmp["E_cm_ld"]/mpi
            tmp["s_ld"] = tmp["E_cm_ld"]**2
            tmp["s_ld_prime"] = tmp["E_cm_ld_prime"]**2
            tmp["pstar_ld"] = pstar_lat_disp(tmp["E_cm_ld"],mpi)
            tmp["pstar_ld_prime"] = tmp["pstar_ld"]/mpi
            tmp["p2star_ld"] = tmp["pstar_ld"]**2
            tmp["p2star_ld_prime"] = tmp["pstar_ld_prime"]**2
            q2_ld = tmp["p2star_ld"]*(N_Ls[i]/(2*np.pi))**2
            tmp["q2_ld"] = q2_ld
            cot_PS_ld = cot_delta_mom(dvecs[i])(q2_ld, N_Ls[i])
            tmp["cot_PS_ld"] = cot_PS_ld
            tmp["tan_PS_ld"] = 1/cot_PS_ld
            PS_ld = 360*np.arctan(1/cot_PS_ld)/(2*np.pi)
            tmp["PS_ld"] = complex(PS_ld.real%180,PS_ld.imag%180)
            tmp["p3cotPS_ld"] = tmp["pstar_ld"]**3*cot_PS_ld
            tmp["p3cotPS_ld_prime"] = tmp["pstar_ld_prime"]**3*cot_PS_ld
            tmp["p3cotPS_Ecm_ld"] = tmp["pstar_ld"]**3/tmp["E_cm_ld"]*cot_PS_ld
            tmp["p3cotPS_Ecm_ld_prime"] = tmp["pstar_ld_prime"]**3/tmp["E_cm_ld_prime"]*cot_PS_ld

            for key, val in tmp.items():
                result[key].append(tmp[key])
    return result


def calc_PS(name):
    info={}
    data = np.transpose(np.genfromtxt("data/%s.dat"%name))
    NL_arr = data[0]
    dvec_arr = np.transpose(data[1:4])
    en_lv_arr = data[4]
    en_arr = data[5]
    en_m_arr = data[6]
    en_p_arr = data[7]
    mpi = data[8][0]
    mrho = data[9][0]
    info["beta"],info["m_1"],info["m_2"], info["mrho"], info["en_lv"] = [6.9,-0.92,-0.92,mrho,en_lv_arr]

    res, res_sampled = result_sampled(info,NL_arr, en_arr, en_m_arr, en_p_arr, dvec_arr, mpi, resampling="lin")
    save_to_hdf(res, res_sampled, name)

if __name__ == "__main__":
    calc_PS("PS_69_092")
    # calc_PS("Plymouth")
    # calc_PS("Lang_Prelovsek")