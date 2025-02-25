import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import csv
from tqdm import tqdm

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=T, tmin=tmin, tmax=tmax, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def first_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p['dE'])
    a = p['a']
    chi2 = fit.chi2     
    dof = fit.dof
    return E, a, chi2, dof

def print_fit_param(fit):
    E, a, chi2, dof = first_fit_parameters(fit) 
    print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
    print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))
    print('chi2/dof = ', chi2/dof, '\n')


def fit_correlator_without_bootstrap(avg,T,tmin,tmax,Nmax,antisymmetric,plotname="test",plotdir="./output/plots/",plotting=False,printing=False):
    T = - abs(T) if antisymmetric else abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax))
    p0 = None
    for N in range(1,Nmax+1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print('nterm =', N, 30 * '=')
            print(fit)

    E, a, chi2, dof = first_fit_parameters(fit) 
    if plotting:
        os.makedirs(plotdir+plotname, exist_ok=True)
        fit.show_plots(view='ratio',save=plotdir+plotname+'/ratio.pdf')
        fit.show_plots(view='log'  ,save=plotdir+plotname+'/data.pdf')
    return E, a, chi2, dof

def fit_all_files(infile,outfile,betas, m0s, Ls, Ts, groups, tmins, tmaxs,binsize=1):

    fid = h5py.File(infile,'r')
    
    print("Fitting correlators.... ")
    for i in tqdm(range(0,len(groups)),disable=True):

        beta = betas[i] 
        m = m0s[i] 
        L = Ls[i]
        T = Ts[i] 
        group = groups[i] 
        tmin = tmins[i]
        tmax = tmaxs[i]

        # read the data from the hdf5 file
        ev       = get_hdf5_value(fid,group+"/eigvals")
        Delta_ev = get_hdf5_value(fid,group+"/Delta_eigvals")
        cov_ev   = get_hdf5_value(fid,group+"/cov_eigvals")

        # Rescale data such that eig(t=0)=1
        # Note: There is an issue when t_min < t0
        #       C(t0) ihas no variance and destabilises the fit
        eig2 = dict(Gab=gv.gvar(ev[:,0]/ev[0,0],Delta_ev[:,0]/ev[0,0]))
        eig1 = dict(Gab=gv.gvar(ev[:,1]/ev[0,1],Delta_ev[:,1]/ev[0,1]))
        # Use full covariance matrix estimator
        eig2 = dict(Gab=gv.gvar(ev[:,0]/ev[0,0],cov_ev[:,:,0]/ev[0,0]/ev[0,0]))
        eig1 = dict(Gab=gv.gvar(ev[:,1]/ev[0,1],cov_ev[:,:,1]/ev[0,1]/ev[0,1]))

        plotname = "beta{}_m{}_L{}_T{}".format(beta,m,L,T)
        antisymmetric = True
        printing=False
        plotting=False
        plotdir = "./output/plots/"
        Nmax = 10

        E1, a1, chi2_1, dof1 = fit_correlator_without_bootstrap(eig1,T,tmin,tmax,Nmax,antisymmetric,plotname,plotdir,plotting,printing)
        E2, a2, chi2_2, dof2 = fit_correlator_without_bootstrap(eig2,T,tmin,tmax,Nmax,antisymmetric,plotname,plotdir,plotting,printing)
        
        f = h5py.File(outfile, "a")
        f.create_dataset(group+"/tmin", data=tmin)
        f.create_dataset(group+"/tmax", data=tmax)
        f.create_dataset(group+"/antisymmetric", data=antisymmetric)
        f.create_dataset(group+"/E0", data=[E_i.mean for E_i in E1])
        f.create_dataset(group+"/E1", data=[E_i.mean for E_i in E2])
        f.create_dataset(group+"/Delta_E0", data=[E_i.sdev for E_i in E1])
        f.create_dataset(group+"/Delta_E1", data=[E_i.sdev for E_i in E2])
        f.create_dataset(group+"/a0", data=[a_i.mean for a_i in a1])
        f.create_dataset(group+"/a1", data=[a_i.mean for a_i in a2])
        f.create_dataset(group+"/Delta_a0", data=[a_i.sdev for a_i in a1])
        f.create_dataset(group+"/Delta_a1", data=[a_i.sdev for a_i in a2])
        f.create_dataset(group+"/chi2_0", data=chi2_1)
        f.create_dataset(group+"/chi2_1", data=chi2_2)
        f.create_dataset(group+"/dof0", data=dof1)
        f.create_dataset(group+"/dof1", data=dof2)
        f.close()
        print(group)
        print(E1[0],chi2_1/dof1)
        print(E2[0],chi2_2/dof2)


def read_filelist_fitparam(parameterfile):
    reader = csv.reader(open(parameterfile))
    # create list that contain the fitting information
    beta  = []
    m0  = []
    L = []
    T = []
    tmins  = []
    tmaxs  = []
    groups = []
    # skip line containing headers
    next(reader, None)
    for row in reader:
        beta.append(row[0])
        m0.append(row[1])
        L.append(int(row[2]))
        T.append(int(row[3]))
        groups.append(row[4])
        tmins.append(int(row[5]))
        tmaxs.append(int(row[6]))

    return beta, m0, L, T, groups, tmins, tmaxs

parameterfile  = './input/pipi_fitintervals.csv'
betas, m0s, Ls, Ts, groups, tmins, tmaxs = read_filelist_fitparam(parameterfile)

infile  = './data/isospin1_eigenvalues_t0_3_deriv.hdf5'
outfile = './data/isospin1_fitresults_t0_3_deriv.hdf5'

if os.path.exists(outfile):
    os.remove(outfile)

fit_all_files(infile,outfile,betas, m0s, Ls, Ts, groups, tmins, tmaxs)