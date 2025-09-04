import numpy as np
import warnings
# from typing import Optional, Tuple, Callable
from scipy.optimize import curve_fit

def safe_curve_fit(func, xdata, ydata, p0=None, name="",err=False, **kwargs):
    """
    Run scipy.optimize.curve_fit safely.

    Parameters
    ----------
    func : callable
        The model function f(x, *params).
    xdata : array_like
        X data to fit.
    ydata : array_like
        Y data to fit.
    p0 : array_like, optional
        Initial guess for parameters.
    name : str
        A string you can pass to identify this fit.
    kwargs : dict
        Extra arguments passed to curve_fit.

    Returns
    -------
    (popt, pcov) : tuple
        Optimal parameters and covariance if fit succeeded.
    None : if fit failed due to warning or error.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("error", category=RuntimeWarning)
        try:
            popt, pcov = curve_fit(func, xdata, ydata, p0=p0, **kwargs)
            # print(f"Fit '{name}' succeeded.")
            return popt, pcov
        except (RuntimeWarning, RuntimeError, ValueError) as e:
            if err:
                print(f"⚠️ Fit '{name}' failed: {e}")
            return None

def cot(x):
    return 1/np.tan(x)

class fitting_model:
    def __init__(self, name, model, param_names):
        self.name = name
        self.model = model
        self.param_names = param_names
        self.num_params = len(param_names)
    def fit(self, xdata, ydata,err=False):
        try:
            popt = safe_curve_fit(self.model, xdata, ydata, np.ones(self.num_params),self.name,err)[0]
            # popt = curve_fit(self.model, xdata, ydata, np.ones(self.num_params))[0]
            res = {}
            for i in range(self.num_params):
                res[self.param_names[i]] = popt[i]
            return res
        except:
            return np.zeros(self.num_params)

def mixed_model(NR_model, R_model):
    def mixed(p2,*args):
        args_NR = args[:NR_model.num_params]
        args_R = args[NR_model.num_params:]
        delta_NR = np.arctan(p2**(3/2)/NR_model.model(p2, *args_NR))
        delta_R = np.arctan(p2**(3/2)/R_model.model(p2, *args_R))
        return p2**(3/2)*cot(delta_NR+delta_R)
    mixed_name = NR_model.name+"__"+R_model.name
    mixed_param_list = [*NR_model.param_names,*R_model.param_names]
    return fitting_model(mixed_name, mixed, [mixed_name+"__"+par_nam for par_nam in mixed_param_list])

### non-resonant models

def ERE_0(p2, a1_0):
    return 0 if a1_0 == 0 else -1/a1_0**3+0*p2
def ERE_1(p2, a1_1, r1_1):
    return 0 if a1_1 == 0 or r1_1 == 0 else  -1/a1_1**3+p2/(2*r1_1)
def ERE_2(p2, a1_2, r1_2, c1_2):
    return 0 if a1_2 == 0 or r1_2 == 0 else  -1/a1_2**3+p2/(2*r1_2)+c1_2*p2**2
def NR_I(p2, A_0):
    return p2**(3/2)*cot(A_0)
def NR_II(p2, A_1, B_1):
    return p2**(3/2)*cot(A_1+B_1*p2)

### resonant models

def ECM_p2(p2):
    return 2*np.sqrt(1+p2)

def BW_I(p2, m_R_I, gVPP2_I):
    ECM = ECM_p2(p2)
    return 0 if gVPP2_I == 0 else 6*np.pi*(m_R_I**2-ECM**2)*ECM/gVPP2_I
def BW_II(p2, m_R_II, gVPP2_II, r0_II): 
    ECM = ECM_p2(p2)
    kR2 = m_R_II**2/4-1
    return 0 if gVPP2_II == 0 else 6*np.pi*(m_R_II**2-ECM**2)*(1+p2*r0_II)/(1+kR2*r0_II)*ECM/gVPP2_II

ERE_0_model = fitting_model("ERE_0",ERE_0,["a1_0"])
ERE_1_model = fitting_model("ERE_1",ERE_1,["a1_1", "r1_1"])
ERE_2_model = fitting_model("ERE_2",ERE_2,["a1_2", "r1_2", "c1_2"])
NR_I_model = fitting_model("NR_I",NR_I,["A_0"])
NR_II_model = fitting_model("NR_II",NR_II,["A_1", "B_1"])
BW_I_model = fitting_model("BW_I",BW_I,["m_R_I", "gVPP2_I"])
BW_II_model = fitting_model("BW_II",BW_II,["m_R_II", "gVPP2_II","r0_II"])

non_res_models = [ERE_0_model,ERE_1_model,ERE_2_model,NR_I_model,NR_II_model]
res_models = [BW_I_model,BW_II_model]
mixed_models = []
for non_res_model in non_res_models:
    for res_model in res_models:
        mixed_models.append(mixed_model(non_res_model,res_model))

all_models = [*non_res_models,*res_models,*mixed_models]

if __name__ == "__main__":
    p2_arr = np.linspace(0,2)
    p3cot_PS_arr = np.ones(len(p2_arr))

    ERE_1_model = fitting_model("ERE_0",ERE_1,["a1_1", "r1_1"])
    BW_II_model = fitting_model("BW_II",BW_II,["m_R_II", "gVPP2_II","r0_II"])
    ERE_1__BW_II__model = mixed_model(ERE_1_model,BW_II_model)

    print(ERE_1_model.name)
    print(BW_II_model.name)
    print(ERE_1__BW_II__model.name)

    print(ERE_1_model.param_names)
    print(BW_II_model.param_names)
    print(ERE_1__BW_II__model.param_names)

    ERE_1_vals = [0.4532,0.123]
    BW_II_vals = [0.213,0.543,0.012]
    mixed_vals = [*ERE_1_vals,*BW_II_vals]

    # print(mixed_vals)
    # exit()

    y_ERE_1 = [ERE_1_model.model(p2,*ERE_1_vals) for p2 in p2_arr]
    y_BW_II = [BW_II_model.model(p2,*BW_II_vals) for p2 in p2_arr]

    def add_PS(p2, p3cot1, p3cot2):
        return p2**(3/2)*cot(np.arctan(p2**(3/2)/p3cot1)+np.arctan(p2**(3/2)/p3cot2))

    y_correct_mix = [add_PS(p2,ERE_1_model.model(p2,*ERE_1_vals),BW_II_model.model(p2,*BW_II_vals)) for p2 in p2_arr]
    y_test_mix = [ERE_1__BW_II__model.model(p2, *mixed_vals) for p2 in p2_arr]

    for i in range(len(y_correct_mix)):
        print(y_correct_mix[i], y_test_mix[i])

    j=0
    for model in all_models:
        print(model.name, j, "\t\t\t", model.fit(p2_arr,p3cot_PS_arr))
        j = j+1
