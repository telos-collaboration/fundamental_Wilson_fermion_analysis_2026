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

def make_models(T,tmin,tmax,t0):
    """ 
    Create corrfitter model for G(t).

    Here I provide a dedicated range of timeslices to be used in fitting
    Using simply (tmin,tmax) is not sufficient if t0 is in (tmin,tmax) 
    because at t=t0 the correlation matrix has a vanishing variance which 
    destabilises the fit.

    Therefore, I create a list 'tfit' of timslices to be included, such that
    t0 is excluded. 

    Furthermore, if T is provided (i.e. it is not 'None') then the symmetry 
    of the correlator is included specifically and the largest t needed is T/2
    
    """
    tfit = range(tmin,tmax+1)
    tfit = list(filter(lambda x: x != t0,tfit))
    if T is not None:
        tfit = list(filter(lambda x: x < abs(T)//2,tfit))

    return [cf.Corr2(datatag='Gab', tp=T, tfit=tfit, a='a', b='a', dE='dE')]

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


def fit_correlator_without_bootstrap(avg,T,tmin,tmax,t0,Nmax,antisymmetric,plotname="test",plotdir="./output/plots/",plotting=False,printing=False):
    T = - abs(T) if antisymmetric else abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax,t0))
    p0 = None
    for N in range(1,Nmax+1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print_fit_param(fit)

    E, a, chi2, dof = first_fit_parameters(fit) 
    if plotting:
        os.makedirs(plotdir+plotname, exist_ok=True)
        fit.show_plots(view='ratio',save=plotdir+plotname+'/ratio.pdf')
        fit.show_plots(view='log'  ,save=plotdir+plotname+'/data.pdf')
    return E, a, chi2, dof

def fit_all_files(infile,outfile,parameterfile):

    fid = h5py.File(infile,'r')
    reader = csv.reader(open(parameterfile))
    next(reader, None) # skip line containing headers
    
    print("Fitting correlators.... ")
    for row in tqdm((reader),disable=False):

        group, tmin, tmax = row[0], int(row[1]), int(row[2]) 

        if group not in fid:
            continue 

        # read the data from the hdf5 file
        ev       = get_hdf5_value(fid,group+"/eigvals")
        Delta_ev = get_hdf5_value(fid,group+"/Delta_eigvals")
        cov_ev   = get_hdf5_value(fid,group+"/cov_eigvals")
        t0       = get_hdf5_value(fid,group+"/t0") - 1 # offset for 1-indexing in julia 
        T        = ev.shape[0]
        antisymmetric = get_hdf5_value(fid,group+"/deriv") 

        # Rescale data such that eig(t=0)=1
        # Note: There was an issue when t_min < t0
        #       C(t0) has no variance and destabilises the fit
        #         t0  is now excluded from the fir
        # Use full covariance matrix estimator
        eig2 = dict(Gab=gv.gvar(ev[:,0]/ev[0,0],cov_ev[:,:,0]/ev[0,0]/ev[0,0]))
        eig1 = dict(Gab=gv.gvar(ev[:,1]/ev[0,1],cov_ev[:,:,1]/ev[0,1]/ev[0,1]))

        plotname = group
        plotdir  = "./output/plots/"
        printing = False
        plotting = False
        Nmax = 10

        E1, a1, chi2_1, dof1 = fit_correlator_without_bootstrap(eig1,T,tmin,tmax,t0,Nmax,antisymmetric,plotname,plotdir,plotting,printing)
        E2, a2, chi2_2, dof2 = fit_correlator_without_bootstrap(eig2,T,tmin,tmax,t0,Nmax,antisymmetric,plotname,plotdir,plotting,printing)
        
        f = h5py.File(outfile, "a")
        f.create_dataset(group+"/tmin", data=tmin)
        f.create_dataset(group+"/tmax", data=tmax)
        f.create_dataset(group+"/Nmax", data=Nmax)
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


infile         = './data/isospin1_eigenvalues_t0_8_deriv.hdf5'
outfile        = './data/isospin1_fitresults_t0_8_deriv.hdf5'
parameterfile  = './input/pipi_fitintervals.csv'

if os.path.exists(outfile):
    os.remove(outfile)

fit_all_files(infile,outfile,parameterfile)