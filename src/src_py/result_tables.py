import numpy as np
import math
import sys
import h5py

def nth(num):
    if num <= 10:
        return 1
    elif num <= 500:
        return num//10
    else:
        return num//100

num_perc = math.erf(1/np.sqrt(2))

def sig_dig(x):
    return -int(math.floor(np.log10(abs(x))))

def round_to_sig(x, y=0):
    if y == 0:
        y = x
    return round(x, sig_dig(y))

def get_str(val, errp, errm=None, sd = None):
    if sd == None:
        sd = sig_dig(min(errp, errm))
    if errm == None:
        errm = errp
    res = ""
    if sd == 0:
        res = "$%i^{+%i}_{-%i}$"%(int(val*10**sd), int(round(errp*10**sd)), int(round(errm*10**sd)))
    else: 
        res = "$(%i^{+%i}_{-%i})\\times 10^{%i}$"%(int(val*10**sd), int(round(errp*10**sd)), int(round(errm*10**sd)), -sd)
    # res = "$%i^{+%i}_{-%i}$"%(int(val*10**sd), int(round(errp*10**sd)), int(round(errm*10**sd)))
    return res

# V=134
# E=13

# print(V, E)
# get_str(V,E)

# exit()

def get_data(h5file_scatter_fit, beta, m0, fit):                # wont work with current get_data()
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
                                    fit_this = False
                                    if fit:
                                        fit_this = hfile[ens][P][irrep][lv]["fit"][()]
                                    if fit_this or not fit:
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

def s_p2(p2):
    return 4*(1+p2)

def p2_s(s):
    return s/4-1

def delta_x(x,p2):
    return np.arctan(p2**(3/2)/x)*360/(2*np.pi)

def delta_res_x(x,p2):
    return np.arctan(p2**(3/2)/(2*np.sqrt(1+p2)*x))*360/(2*np.pi)

# def delete_zeros(value, error):
#     if value == 0:
#         return "--"
#     else:
#         return format_with_error(value, error)
    
# def format_with_error(value: float, error: float, two_digits=False) -> str:
#     """
#     Format a value with uncertainty in LaTeX style: x = value(error).
    
#     Parameters:
#         value (float): Central value
#         error (float): Uncertainty
#         two_digits (bool): If True, show two significant digits in the error
    
#     Returns:
#         str: formatted string, e.g. "132.124(56)"
#     """
#     if error <= 0:
#         raise ValueError("Uncertainty must be positive")

#     # Determine number of significant digits for the error
#     digits = 2 if two_digits else 1

#     # Get the order of magnitude of the error
#     exponent = int(math.floor(math.log10(error)))
    
#     # Round error to `digits` significant figures
#     rounded_error = round(error, -exponent + (digits - 1))

#     # Round value to the same decimal place
#     decimal_places = -exponent + (digits - 1)
#     rounded_value = round(value, decimal_places+1)%180

#     # Scale the error for the parentheses
#     error_str = f"{int(rounded_error * 10**(-exponent + (digits - 1)))}"

#     # Format with correct number of decimals
#     value_str = f"{rounded_value:.{max(decimal_places,0)}f}"

#     return f"{value_str}({error_str})"

# print(format_with_error(140.1234123,3.21312))
# exit()

def irrep_str(irrep):
    if irrep == "T1":
        return "$T_1$"
    elif irrep == "A1":
        return "$A_1$"
    elif irrep == "B1":
        return "$B_1$"
    elif irrep == "E":
        return "$E$"
    else:
        raise ValueError

def result_tables(h5file, beta, m0, TBLDIR):
    info, info_nf, fit_param_mean, fit_param_spl, scat_mean, scat_spl, scat_nf_mean, scat_nf_spl = get_data(h5file, beta, m0, False)

    lvs = np.asarray(info["lv"])
    for i in range(len(lvs)):
        lvs[i] = lvs[i] + 1
    irreps = np.asarray(scat_mean["irrep"])
    N_Ls = np.asarray(scat_mean["N_L"])
    dvecs = scat_mean["dvec"]
    dvecs = [[int(x.decode("utf-8")[0]),int(x.decode("utf-8")[1]),int(x.decode("utf-8")[2])] for x in dvecs]
    d2s = np.asarray([np.dot(d,d) for d in dvecs])

    E_m = np.asarray(scat_mean["aEn"])
    E_s = np.asarray(scat_spl["aEn"])
    s_m = np.asarray(scat_mean["E_cm"])
    s_s = np.asarray(scat_spl["E_cm"])
    PS_m = np.asarray(scat_mean["PS"])
    PS_s = np.asarray(scat_spl["PS"])


    length = len(E_s[0])

    print(PS_s.shape)

    def mean_err(arr):
        mean = np.asarray([sorted(arr[i])[length//2-1] for i in range(len(arr))])
        err = np.asarray([(abs(sorted(arr[i])[math.ceil(length*(1+num_perc)/2)] - mean[i]) + 
                           abs(sorted(arr[i])[math.floor(length*(1-num_perc)/2)] - mean[i]))/4 for i in range(len(arr))])
        # fac=10**20
        # err_tmp = [round(err[i]*fac,0)/fac for i in range(len(err))]
        # for i in range(len(err)):
        #     print(err[i], err_tmp[i])
        return mean, err

    # E_err = [abs(sorted(E_s[i])[math.ceil(length*(1+num_perc)/2)]-sorted(E_s[i])[math.floor(length*(1-num_perc)/2)])/2 for i in range(len(E_s))]
    # s_err = [abs(sorted(s_s[i])[math.ceil(length*(1+num_perc)/2)]-sorted(s_s[i])[math.floor(length*(1-num_perc)/2)])/2 for i in range(len(E_s))]
    # PS_err = [abs(sorted(PS_s[i])[math.ceil(length*(1+num_perc)/2)]-sorted(PS_s[i])[math.floor(length*(1-num_perc)/2)])/2 for i in range(len(E_s))]

    E_m, E_err = mean_err(E_s)
    s_m, s_err = mean_err(s_s)
    PS_m, PS_err = mean_err(PS_s)

    start_str = "\\begin{tabular}{|l|l|l|c|l|l|c|c|}\n \\hline \n"
    end_str = "\\end{tabular}"
    first_line = "\t\t$N_L$ & $\\vec{P}$ & $\Lambda$ & n & $aE_n^{\\Lambda,\\,\\vec{P}}$ & $a\sqrt{s_n^{\Lambda,\\,\\vec{P}}}$ & $\delta_1$   & Incl. \\\\ \\hline \\hline\n"
    latex_seperator =  " & "

    with open(TBLDIR+"latex_scattering_b%1.2f_m%1.3f_table.txt"%(beta,m0),"w") as f:
        f.write(start_str)
        f.write(first_line)
        # f.write(after_first_line_str)
        for i in range(len(E_m)):
            line_str = "\t\t"
            line_str += "%i"%N_Ls[i] + latex_seperator
            line_str += "%i"%d2s[i] + latex_seperator
            line_str += irrep_str(irreps[i]) + latex_seperator
            line_str += "%i"%lvs[i] + latex_seperator
            line_str += "%f\t%f"%(E_m[i],E_err[i]) + latex_seperator
            line_str += "%f\t%f"%(s_m[i],s_err[i]) + latex_seperator
            line_str += "%f\t%f"%(PS_m[i],PS_err[i]) + latex_seperator
            # line_str += format_with_error(E_m[i],E_err[i]) + latex_seperator
            # line_str += format_with_error(s_m[i],s_err[i]) + latex_seperator
            # line_str += delete_zeros(PS_m[i],PS_err[i]) + latex_seperator
            line_str += "Yes"

            f.write(line_str)
            f.write("\\\\")
            if i < len(E_m)-1:
                if N_Ls[i+1] > N_Ls[i]:
                    f.write("  \\hline") 
            f.write("\n")
        f.write("\t\\hline\n")
        f.write(end_str)

if __name__ == "__main__":

    args = sys.argv
    TBLDIR = args[1]
    h5file  = args[2]

    result_tables(h5file, 6.9, -0.92, TBLDIR)
    result_tables(h5file, 7.05, -0.863, TBLDIR)
    result_tables(h5file, 7.05, -0.867, TBLDIR)