import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import os
import os.path as op
import matplotlib.pyplot as plt
import csv
import sys
from tqdm import tqdm

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax):
    """  Create corrfitter model for G(t). """
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


def fit_correlator_without_bootstrap(avg,T,tmin,tmax,Nmax,antisymmetric,plotname="test",plotdir="./tmp/plots/",plotting=False,printing=False):
    T = - abs(T) if antisymmetric else abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax))
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
    lines  = sum(1 for row in csv.reader(open(parameterfile))) - 1
    reader = csv.reader(open(parameterfile))
    next(reader, None) # skip line containing headers
    
    print("Fitting of correlators started...")
    for row in tqdm((reader), total=lines , desc="Fit eigenvalues", disable=True):

        ensemble, p, irrep, tmin, tmax, Nmax = row[0], row[1], row[2], int(row[3]), int(row[4]), int(row[5])
        print(f"{ensemble},{irrep},{p}")

        lattice = get_hdf5_value(fid,op.join(ensemble,"lattice"))
        T = lattice[0]

        plotname = op.join(ensemble,irrep)
        plotdir  = "./data_assets/plots/"
        printing = False
        plotting = False
        
        f = h5py.File(outfile, "a")
        

        corr_label = "Cpi" if irrep == "pi" else f"{irrep}/C"
        covm_label = "cov_Cpi" if irrep == "pi" else f"{irrep}/cov_C"

        C = get_hdf5_value(fid,op.join(ensemble,p,corr_label))
        cov = get_hdf5_value(fid,op.join(ensemble,p,covm_label))
        var = gv.gvar(C[:],cov[:,:])
        cor = dict(Gab=var/var[0])
            
        E, a, chi2, dof = fit_correlator_without_bootstrap(cor,T,tmin,tmax,Nmax,False,plotname,plotdir,plotting,printing)
        
        f.create_dataset(op.join(ensemble,p,irrep,"E"), data=[E[0].mean])
        f.create_dataset(op.join(ensemble,p,irrep,"Delta_E"), data=[E[0].sdev])
        f.create_dataset(op.join(ensemble,p,irrep,"a"), data=[a[0].mean])
        f.create_dataset(op.join(ensemble,p,irrep,"Delta_a"), data=[a[0].sdev])
        f.create_dataset(op.join(ensemble,p,irrep,"chi2"), data=chi2)
        f.create_dataset(op.join(ensemble,p,irrep,"dof"), data=dof)
        f.create_dataset(op.join(ensemble,p,irrep,"Nmax"), data=Nmax)
        f.create_dataset(op.join(ensemble,p,irrep,"tmin"), data=tmin)
        f.create_dataset(op.join(ensemble,p,irrep,"tmax"), data=tmax)

        if not op.join(ensemble,"lattice") in f.keys():
            f.create_dataset(op.join(ensemble,"lattice"), data=lattice)

        f.close()
    print("Fitting of correlators finished.")

args = sys.argv
infile         = args[1]
outfile        = args[2]
parameterfile  = args[3]

fit_all_files(infile,outfile,parameterfile)