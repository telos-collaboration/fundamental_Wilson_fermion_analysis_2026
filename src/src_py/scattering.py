import numpy as np
import h5py
import sys
from pylink_wlm import wlm_func_c as wlm
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

num_resample=None

def save_to_hdf(res,res_sample, info, ens, P, irrep, lv, outfile):
    group = ens+"/"+P+"/"+irrep+"/"+"lv"+"%i/"%lv
    with h5py.File(outfile,"a") as hfile:
        for key, val in res.items():
            hfile.create_dataset(group+"mean/"+key, data = val)
        for key, val in res_sample.items():
            hfile.create_dataset(group+"sample/"+key, data = val)
        for key, val in info.items():
            hfile.create_dataset(group+"info/"+key, data = val)

def result_sampled(N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,mpi,irrep,ld=False,resampling="gauss",num_resample=50):
    res = {}
    info = {}
    res_sample = {}
    info["resampling"] = resampling
    info["num_resample"] = num_resample
    info["L_prime"] = N_L*mpi
    info["dvec"] = dvec
    info["d"] = [np.sqrt(np.dot(dvec[i],dvec[i])) for i in range(len(dvec))]
    info["d2"] = [np.dot(dvec[i],dvec[i]) for i in range(len(dvec))]

    for key, val in get_rizz(E_pipi,N_L,dvec,mpi,irrep,ld).items():
        res[key] = val
    for key in res.keys():
        res_sample[key] = []

    if resampling == "gauss":
        for i in tqdm(range(num_resample)):
            if i < num_resample//2:                                     # I insert a bias here!! Check if we actually want that!
                E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
            else:
                E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep,ld)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    elif resampling == "lin":
        for E_pipi_tmp in tqdm(np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_resample)):
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep,ld)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    return res, res_sample, info

def result_sampled_parallel(N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,mpi,irrep,ld):
    # ld=args[6] == "True"
    resampling=args[5]
    num_res= int(args[4])
    res = {}
    info = {}
    res_sample = {}
    info["resampling"] = resampling
    info["num_resample"] = num_res
    info["L_prime"] = N_L*mpi
    info["dvec"] = dvec
    info["d"] = [np.sqrt(np.dot(dvec[i],dvec[i])) for i in range(len(dvec))]
    info["d2"] = [np.dot(dvec[i],dvec[i]) for i in range(len(dvec))]


    for key, val in get_rizz(E_pipi,N_L,dvec,mpi,irrep,ld).items():
        res[key] = val
    for key in res.keys():
        res_sample[key] = []

    if resampling == "gauss":
        for i in range(num_res):
            if i < num_res//2:
                E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
            else:
                E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep,ld)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    elif resampling == "lin":
        for E_pipi_tmp in np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_res):
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi,irrep,ld)
            for key, val in res_tmp.items():
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
    # p = [0,1,1]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = 2*wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first + second # sign of third term depends on convention. Please check!!
def cot_delta_110_B2(q2,N_L,mpi):                                                       # A1=B2 (B3 in Luka)
    p = [1,1,0]
    first = wlm(0,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    second = wlm(2,0,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    third = complex(0,1) * np.sqrt(6) * wlm(2,2,p[0],p[1],p[2],mpi,mpi,q2,int(N_L))
    return first - second + third

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
        elif irrep == "B2":
            return cot_delta_110_B2
    elif list(dvec) == [0,1,1]:
        if irrep == "A1":
            return cot_delta_110_A1
        elif irrep == "B1":
            return cot_delta_110_B1
    elif list(dvec) == [1,1,1]:
        if irrep == "A1":
            return cot_delta_111_A1
        elif irrep == "E":
            return cot_delta_111_E
    else:
        print("wrong momentum or irrep")
        raise ValueError("wrong momentum or irrep in cot_delta_mom")

def Ecm(ld):
    if ld:
        return Ecm_ld
    else:
        return Ecm_cont

def Ecm_cont(E, Pvec):
    return np.sqrt(E**2-np.dot(Pvec,Pvec))

def Ecm_ld(E, Pvec):
    return np.arccosh(np.cosh(E)-2*(np.sin(Pvec[0]/2)**2+np.sin(Pvec[1]/2)**2+np.sin(Pvec[2]/2)**2))

def pstar(ld):
    if ld:
        return pstar_ld
    else:
        return pstar_cont

def pstar_cont(Ecm, mpi):
    return np.sqrt(Ecm**2/4-mpi**2)

def pstar_ld(Ecm, mpi):
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(Ecm/2)-np.cosh(mpi))))

def get_rizz(E_pipi, N_L, dvec, mpi, irrep,ld=False):
    res = {}
    key_list = ["aEn","En_prime","E_cm","E_cm_prime","s","s_prime","pstar","pstar_prime","p2star","p2star_prime","q","q2","cot_PS","tan_PS","PS", "p3cotPS", "p3cotPS_prime", "p3cotPS_Ecm", "p3cotPS_Ecm_prime", "sigma", "sigma_prime"]

    Pvec = [2*np.pi*x/N_L for x in dvec]
    En_prime = E_pipi/mpi
    if E_pipi**2 >= np.dot(Pvec,Pvec):
        E_cm_prime = Ecm(ld)(E_pipi,Pvec)/mpi
    else:
        E_cm_prime = 0
    if E_cm_prime < 2 or E_cm_prime > 4:
        for key in key_list:
            res[key] = float(0)
    else:
        res["E_cm_prime"] = E_cm_prime
        res["E_cm"] = E_cm_prime*mpi
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
        res["sigma"] = 12*np.pi*res["p2star"]**2/(res["p2star"]**3+res["p3cotPS"]**2)
        res["sigma_prime"] = res["sigma"]*mpi**2

    res["aEn"] = E_pipi
    res["En_prime"] = En_prime
    res["E_cm_prime"] = E_cm_prime
    res["E_cm"] = E_cm_prime*mpi
    res["s_prime"] = res["E_cm_prime"]**2
    res["s"] = res["E_cm"]**2
    res["dvec"] = "%i%i%i"%(dvec[0],dvec[1],dvec[2])
    res["N_L"] = N_L
    return res

def calc_all_phaseshifts(input_file, fitresults, h5file, resampling="lin",num_resample=10):
    infile = np.transpose(np.genfromtxt(input_file,delimiter=";",skip_header=1,dtype=str))
    with h5py.File(fitresults,"r") as hfile:
        for ens in hfile:
            for P in hfile[ens]:
                if P[0] == "p":
                    dvec = [int(P[2]),int(P[4]),int(P[6])]
                    #NOTE: Somethins is broken here with the extra irreps
                    for irrep in hfile[ens][P]:
                        if irrep != "pi":
                            beta, m0, mpi, mrho, ld, num_lv = infile[1:,infile[0] == ens+P+irrep] 
                            beta = float(beta[0])
                            m0 = float(m0[0])
                            mpi = float(mpi[0])
                            mrho = float(mrho[0])
                            ld = ld[0] == "True"
                            num_lv = int(num_lv[0])
                            info = {}
                            info["beta"] = beta
                            info["m0"] = m0
                            info["mpi"] = mpi
                            info["mrho"] = mrho
                            NL = hfile[ens]["lattice"][()][3]
                            info["NL"] = NL
                            for i in range(num_lv):
                                E = hfile[ens][P][irrep]["E"][()][i]
                                E_m = hfile[ens][P][irrep]["Delta_E"][()][i]
                                E_p = hfile[ens][P][irrep]["Delta_E"][()][i]
                                res, res_sampled, info_tmp = result_sampled(NL, E, E_m, E_p, dvec, mpi, irrep, ld, resampling=resampling, num_resample=num_resample)
                                for key, val in info_tmp.items():
                                    info[key] = val
                                save_to_hdf(res, res_sampled, info, ens, P, irrep, i, h5file)

def unpack_and_run(args):
    return result_sampled_parallel(*args)

def calc_all_phaseshifts_parallel(input_file, fitresults, h5file):
    info=[]
    NLs, Es, E_ms, E_ps, dvecs, mpis, irreps, enss, Ps, lvs, lds = [[],[],[],[],[],[],[],[],[],[],[]]
    infile = np.transpose(np.genfromtxt(input_file,delimiter=";",skip_header=1,dtype=str))
    with h5py.File(fitresults,"r") as hfile:
        for ens in hfile:
            for P in hfile[ens]:
                if P[0] == "p":
                    dvec = [int(P[2]),int(P[4]),int(P[6])]
                    #NOTE: Somethins is broken here with the extra irreps
                    for irrep in hfile[ens][P]:
                        if irrep != "pi":
                            # print(ens,P,irrep)
                            beta, m0, mpi, mrho, ld, num_lv = infile[1:,infile[0] == ens+P+irrep] 
                            beta = float(beta[0])
                            m0 = float(m0[0])
                            mpi = float(mpi[0])
                            mrho = float(mrho[0])
                            ld = ld[0] == "True"
                            num_lv = int(num_lv[0])
                            for i in range(num_lv):
                                info.append({})
                                info[-1]["beta"] = beta
                                info[-1]["m0"] = m0
                                info[-1]["mpi"] = mpi
                                info[-1]["ld"] = ld
                                info[-1]["mrho"] = mrho
                                NL = hfile[ens]["lattice"][()][3]
                                info[-1]["NL"] = NL
                                E = hfile[ens][P][irrep]["E"][()][i]
                                E_m = hfile[ens][P][irrep]["Delta_E"][()][i]
                                E_p = hfile[ens][P][irrep]["Delta_E"][()][i]
                                NLs.append(NL)
                                Es.append(E)
                                E_ms.append(E_m)
                                E_ps.append(E_p)
                                dvecs.append(dvec)
                                mpis.append(mpi)
                                irreps.append(irrep)
                                enss.append(ens)
                                Ps.append(P)
                                lvs.append(i)
                                lds.append(ld)
        
    params = list(zip(NLs, Es, E_ms, E_ps, dvecs, mpis, irreps,lds))
    
    results = [None] * len(params)    
    num_cores = os.cpu_count()
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = {executor.submit(unpack_and_run, p): i for i, p in enumerate(params)}
        for f in tqdm(as_completed(futures), total=len(futures)):
            idx = futures[f]
            results[idx] = f.result()
    for i in range(len(results)):
        res, res_sampled, info_tmp = results[i]
        for key, val in info_tmp.items():
            info[i][key] = val
        save_to_hdf(res, res_sampled, info[i], enss[i], Ps[i], irreps[i], lvs[i], h5file)

if __name__ == "__main__":
    # avod hard-coding of names outside of main
    args = sys.argv
    input_file = args[1]
    fitresults = args[2]
    h5fileout  = args[3]
    num_res= int(args[4])
    resampling=args[5]

    # calc_all_phaseshifts(input_file, fitresults, h5fileout,resampling,num_res)
    calc_all_phaseshifts_parallel(input_file, fitresults, h5fileout)