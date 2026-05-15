import numpy as np
import warnings
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
            return popt, pcov
        except (RuntimeWarning, RuntimeError, ValueError) as e:
            if err:
                print(f"⚠️ Fit '{name}' failed: {e}")
            return None

def cot(x):
    return 1/np.tan(x)

class fitting_model:
    def __init__(self, name, model, param_names, xaxis, yaxis, p0 = None):
        self.name = name
        self.model = model
        self.param_names = param_names
        self.num_params = len(param_names)
        if p0 == None:
            p0 = np.ones(self.num_params)
        self.p0 = p0
        self.xaxis = xaxis
        self.yaxis = yaxis
    def fit(self,xaxis,yaxis,err=False):
        try:
            popt = safe_curve_fit(self.model, xaxis, yaxis, self.p0,self.name,err,maxfev=4000)[0]
            res = {}
            for i in range(self.num_params):
                res[self.param_names[i]] = popt[i]
            return res
        except:
            res = {}
            for i in range(self.num_params):
                res[self.param_names[i]] = 0
            return res

def mixed_model(NR_model, R_model):
    def mixed(p2,*args): # This is adding the p3cotPS
        args_NR = args[:NR_model.num_params]
        args_R = args[NR_model.num_params:]
        return NR_model.model(p2, *args_NR) + R_model.model(p2, *args_R)
    mixed_name = NR_model.name+"__"+R_model.name
    mixed_param_list = [*NR_model.param_names,*R_model.param_names]
    return fitting_model(mixed_name, mixed, [mixed_name+"__"+par_nam for par_nam in mixed_param_list])

### non-resonant models

def ERE_0(p2, a1_0):
    return a1_0+0*p2
def ERE_1(p2, a1_1, r1_1):
    return a1_1+p2*r1_1
def ERE_2(p2, a1_2, r1_2, c1_2):
    return a1_2+p2*r1_2+p2**2*c1_2
def NR_I(s, A_0):
    return A_0+0*s
def NR_II(s, A_1, B_1):
    return A_1+B_1*s

### resonant models

def BW_I(s, m_R_I, gVPP2_I):
    return (m_R_I**2-s)*gVPP2_I

def BW_II(s, m_R_II, gVPP2_II, r02_II): 
    kR2 = m_R_II**2/4-1
    p2 = s/4-1
    return (m_R_II-s)*gVPP2_II*(1+p2*r02_II)/(1+kR2*r02_II)

def BW_I_PS(s, exp_mR2_BW_I_PS, gVPP2_BW_I_PS): 
    return ((180/(np.pi))*np.arctan((gVPP2_BW_I_PS*(s/4-1)**(3/2))/(6*np.pi*np.sqrt(s)*(4+np.exp(exp_mR2_BW_I_PS)-s))))%180

def BW_II_PS(s, exp_MR2_BW_II_PS, gVPP2_BW_II_PS, r02_BW_II_PS): 
    return ((180/(np.pi))*np.arctan((gVPP2_BW_II_PS*(s/4-1)**(3/2)*(4+exp_MR2_BW_II_PS*r02_BW_II_PS))/(6*np.pi*np.sqrt(s)*(4+np.exp(exp_MR2_BW_II_PS)-s)*(4+r02_BW_II_PS*(s-4)))))%180

def BW_I_NR_I_PS(s, exp_mR2_BW_I_NR_I_PS, gVPP2_BW_I_NR_I_PS, A_BW_I_NR_I): 
    return ((180/(np.pi))*np.arctan((gVPP2_BW_I_NR_I_PS*(s/4-1)**(3/2))/(6*np.pi*np.sqrt(s)*(4+np.exp(exp_mR2_BW_I_NR_I_PS)-s)))+A_BW_I_NR_I)%180

def BW_I_NR_II_PS(s, exp_mR2_BW_I_NR_II_PS, gVPP2_BW_I_NR_II_PS, A_BW_I_NR_II, B_BW_I_NR_II): 
    return ((180/(np.pi))*np.arctan((gVPP2_BW_I_NR_II_PS*(s/4-1)**(3/2))/(6*np.pi*np.sqrt(s)*(4+np.exp(exp_mR2_BW_I_NR_II_PS)-s)))+A_BW_I_NR_II+B_BW_I_NR_II*s)%180

def BW_II_NR_I_PS(s, exp_MR2_BW_II_NR_I_PS, gVPP2_BW_II_NR_I_PS, r02_BW_II_NR_I_PS, A_BW_II_NR_I): 
    return ((180/(np.pi))*np.arctan((gVPP2_BW_II_NR_I_PS*(s/4-1)**(3/2)*(4+exp_MR2_BW_II_NR_I_PS*r02_BW_II_NR_I_PS))/(6*np.pi*np.sqrt(s)*(4+np.exp(exp_MR2_BW_II_NR_I_PS)-s)*(4+r02_BW_II_NR_I_PS*(s-4))))+A_BW_II_NR_I)%180

ERE_0_model = fitting_model("ERE_0",ERE_0,["a1_0",],"p2star_prime", "p3cotPS_prime")
ERE_1_model = fitting_model("ERE_1",ERE_1,["a1_1", "r1_1"],"p2star_prime", "p3cotPS_prime")
ERE_2_model = fitting_model("ERE_2",ERE_2,["a1_2", "r1_2", "c1_2"],"p2star_prime", "p3cotPS_prime")
NR_I_model = fitting_model("NR_I",NR_I,["A_0",],"s_prime", "PS")
NR_II_model = fitting_model("NR_II",NR_II,["A_1", "B_1"],"s_prime", "PS")
BW_I_model = fitting_model("BW_I",BW_I,["m_R_I", "gVPP2_I"],"s_prime", "p3cotPS_Ecm_prime")
BW_II_model = fitting_model("BW_II",BW_II,["m_R_II", "gVPP2_II","r02_II"],"s_prime", "p3cotPS_Ecm_prime")
BW_I_PS_model = fitting_model("BW_I_PS",BW_I_PS,["exp_mR2_BW_I_PS", "gVPP2_BW_I_PS"],"s_prime", "PS")
BW_II_PS_model = fitting_model("BW_II_PS",BW_II_PS,["exp_MR2_BW_II_PS","gVPP2_BW_II_PS","r02_BW_II_PS"],"s_prime", "PS")
BW_I_NR_I_PS_model = fitting_model("BW_I_NR_I_PS",BW_I_NR_I_PS,["exp_mR2_BW_I_NR_I_PS","gVPP2_BW_I_NR_I_PS","A_BW_I_NR_I"],"s_prime", "PS")
BW_I_NR_II_PS_model = fitting_model("BW_I_NR_II_PS",BW_I_NR_II_PS,["exp_mR2_BW_I_NR_II_PS","gVPP2_BW_I_NR_II_PS","A_BW_I_NR_II","B_BW_I_NR_II"],"s_prime", "PS")
BW_II_NR_I_PS_model = fitting_model("BW_II_NR_I_PS",BW_II_NR_I_PS,["exp_MR2_BW_II_NR_I_PS","gVPP2_BW_II_NR_I_PS","r02_BW_II_NR_I_PS","A_BW_II_NR_I"],"s_prime", "PS")

non_res_models = [ERE_0_model,ERE_1_model,ERE_2_model,NR_I_model,NR_II_model]
res_models = [BW_I_model,BW_II_model,BW_I_PS_model,BW_II_PS_model,BW_I_NR_I_PS_model,BW_I_NR_II_PS_model,BW_II_NR_I_PS_model]

all_models = [*non_res_models,*res_models]

