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


def fit_correlator_without_bootstrap(avg,T,tmin,tmax,Nmax,plotname="test",plotdir="./output/plots/",plotting=True,printing=False,antisymmetric=False):
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
    if True:
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

        print(ev.shape)
        Nops = ev.shape[1]
        eig1 = dict(Gab=gv.gvar(ev[:,Nops-1],Delta_ev[:,Nops-1]))
        eig2 = dict(Gab=gv.gvar(ev[:,Nops-2],Delta_ev[:,Nops-2]))

        # data needed for the pion decay constant
        # (normalisation of correlator is important here!)
        plotname = "beta{}_m{}_L{}_T{}".format(beta,m,L,T)

        antisymmetric = False
        plotdir = "./output/plots/"
        Nmax = 10

        #fit_correlator_without_bootstrap(avg,T,tmin,tmax,Nmax,tp,plotting=False,printing=False):
        E1, a1, chi2_1, dof1 = fit_correlator_without_bootstrap(eig1,T,tmin,tmax,Nmax,plotname,plotdir,antisymmetric)
        E2, a2, chi2_2, dof2 = fit_correlator_without_bootstrap(eig2,T,tmin,tmax,Nmax,plotname,plotdir,antisymmetric)
        print(group)
        print(E1[0])
        print(E2[0])


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

infile  = './data/isospin1_eigenvalues_t0_8.hdf5'
outfile = './data/isospin1_test.hdf5'
fit_all_files(infile,outfile,betas, m0s, Ls, Ts, groups, tmins, tmaxs)