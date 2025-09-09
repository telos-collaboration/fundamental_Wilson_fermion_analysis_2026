import numpy as np
import h5py
import sys
from scipy.optimize import curve_fit
# from tqdm import tqdm
import fit_models as fm

import warnings
warnings.simplefilter("error")


def curve_fit_try(func, x, y, num_res):
    try:
        popt, pcov = curve_fit(func, x, y)
        return popt
    except:
        return np.zeros(num_res)

def get_fits(data,err=True):
    res = {}
    for model in fm.all_models:
        xaxis = data[model.xaxis]
        yaxis = data[model.yaxis]
        tmp = model.fit(xaxis,yaxis,err)
        for key, val in tmp.items():
            res[key] = val
    return res

# def get_fits_sample_again(data_spl,err=True):         # this did not work
#     res = {}
#     for model in fm.all_models:
#         xtmp = data_spl[model.xaxis]
#         ytmp = data_spl[model.yaxis]
#         xaxis = []
#         yaxis = []
#         for i in range(len(xtmp)):
#             for j in range(len(xtmp[i])):
#                 xaxis.append(xtmp[i,j])
#                 yaxis.append(ytmp[i,j])
#         print(len(xaxis))
#         # print(len(xaxis[0]))
#         # exit()
#         tmp = model.fit(xaxis,yaxis,err)
#         for key, val in tmp.items():
#             res[key] = val
#     return res

def get_fits_spl(data_spl,res,err=False):
    res_spl = {}
    for key in res:
        res_spl[key] = []
    for model in fm.all_models:
        xaxis_spl = data_spl[model.xaxis]
        yaxis_spl = data_spl[model.yaxis]
        # print(len(xaxis_spl),len(xaxis_spl[0]))
        # exit()
        for i in range(len(xaxis_spl)):
            tmp = model.fit(xaxis_spl[i],yaxis_spl[i],err)
            for key, val in tmp.items():
                res_spl[key].append(val)
    return res_spl
    
def genfromtxt_skip_empty(filename, **kwargs):
    with open(filename) as f:   
        # keep only non-empty lines that don't start with "#"
        lines = [line for line in f if line.strip() and not line.lstrip().startswith("#")]
        return np.genfromtxt(lines, **kwargs)

def fit_one_phaseshift(h5file_out, input_file, beta, m0):
    res_scat = {}
    res_spl_scat = {}

    infile = np.transpose(genfromtxt_skip_empty(input_file,delimiter=";",skip_header=1,dtype=str))

    with h5py.File(h5file_out,"a") as hfile:
        for ens in hfile:
            if str(beta) in ens and str(m0) in ens:
                for P in hfile[ens]:
                    if P[0] == "p":
                        for irrep in hfile[ens][P]:
                            for lv in hfile[ens][P][irrep]:
                                if lv[:2] == "lv":
                                    # print(ens+P+irrep+lv)
                                    fit_in = infile[1][infile[0] == ens+P+irrep+lv]
                                    fit = None
                                    if len(fit_in) == 0:
                                        fit = "False"
                                    else:
                                        fit = fit_in[0]
                                    if fit == "True":
                                        hfile[ens][P][irrep][lv]["fit"] = True
                                        p2star_prime = hfile[ens][P][irrep][lv]["mean"]["p2star_prime"][()]
                                        # print(p2star_prime)
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
                                            raise ValueError("Energy not in elastic window at: %s"%(ens+P+irrep+lv))
                                    elif fit == "False":
                                        hfile[ens][P][irrep][lv]["fit"] = False
                                        pass
                                    else:
                                        raise RuntimeError("Wrong assignment in 'fit_scatter_input.csv' at: %s"%(ens+P+irrep+lv))

        if res_scat == {} or len(res_scat["p2star_prime"]) < 3:
            raise RuntimeError("Less than 3 energy levels selected in 'fit_scatter_input.csv' for beta=%f and m0=%f"%(beta, m0))

        for key in res_scat:
            res_scat[key] = np.asarray(res_scat[key])
            res_spl_scat[key] = np.transpose(np.asarray(res_spl_scat[key]))

        fit_beta_m = "fit_b%f_m%f"%(beta,m0)

        print(len(res_scat["p2star_prime"]))

        res_fit = get_fits(res_scat,err=True)
        res_spl_fit = get_fits_spl(res_spl_scat,res_fit,err=False)

        for key, val in res_fit.items():
            mean_group = hfile.require_group(fit_beta_m+"/mean")
            mean_group.create_dataset(key, data=val)
        for key, val in res_spl_fit.items():
            spl_group = hfile.require_group(fit_beta_m+"/sample")
            spl_group.create_dataset(key, data=val)

def fit_all_phase_shifts(h5file, input_file):
    print("Fitting phase shifts...")
    fit_one_phaseshift(h5file, input_file,6.9,-0.92)
    # fit_one_phaseshift(h5file, input_file,7.05,-0.863)
    fit_one_phaseshift(h5file, input_file,7.05,-0.867)
    print("Done!")


if __name__ == "__main__":

    args = sys.argv
    h5file = args[1]
    input_file = args[2]

    fit_all_phase_shifts(h5file,input_file)