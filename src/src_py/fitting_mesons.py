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

        ensemble, p, irrep, tmin1, tmax1 = row[0], row[1], row[2], int(row[3]), int(row[4])
        Nmax = int(row[12])
        print(f"{ensemble},{irrep},{p}")
        # read the data from the hdf5 file
        lattice = get_hdf5_value(fid,op.join(ensemble,"lattice"))
        T = lattice[0]

        plotname = op.join(ensemble,irrep)
        plotdir  = "./data_assets/plots/"
        printing = False
        plotting = False
        
        f = h5py.File(outfile, "a")
        
        # Fit pion correlator (if it exists)
        if "Cpi" in fid[op.join(ensemble,p)].keys():
            Cpi = get_hdf5_value(fid,op.join(ensemble,p,"Cpi"))
            cov_pi = get_hdf5_value(fid,op.join(ensemble,p,"cov_Cpi"))
            var_pi = gv.gvar(Cpi[:],cov_pi[:,:])
            cor_pi = dict(Gab=var_pi/var_pi[0])
            E, a, chi2, dof = fit_correlator_without_bootstrap(cor_pi,T,1,T//2-1,Nmax,False,plotname,plotdir,plotting,printing)
            f.create_dataset(op.join(ensemble,p,"pi/E"), data=[E[0].mean])
            f.create_dataset(op.join(ensemble,p,"pi/Delta_E"), data=[E[0].sdev])
            f.create_dataset(op.join(ensemble,p,"pi/a"), data=[a[0].mean])
            f.create_dataset(op.join(ensemble,p,"pi/Delta_a"), data=[a[0].sdev])
            f.create_dataset(op.join(ensemble,p,"pi/chi2"), data=chi2)
            f.create_dataset(op.join(ensemble,p,"pi/dof"), data=dof)
            f.create_dataset(op.join(ensemble,p,"pi/Nmax"), data=Nmax)
            f.create_dataset(op.join(ensemble,p,"pi/tmin1"), data=1)
            f.create_dataset(op.join(ensemble,p,"pi/tmax1"), data=T//2-1)

        # If it hasn't been done yet: Fit zero-momentum pi and rho correlator
        p0 = "p(0,0,0)"
        if op.join(ensemble,p0,"pi/E") not in f.keys():
            Cpi = get_hdf5_value(fid,op.join(ensemble,p0,"Cpi"))
            cov_pi = get_hdf5_value(fid,op.join(ensemble,p0,"cov_Cpi"))
            var_pi = gv.gvar(Cpi[:],cov_pi[:,:])
            cor_pi = dict(Gab=var_pi/var_pi[0])
            E, a, chi2, dof = fit_correlator_without_bootstrap(cor_pi,T,1,T//2-1,Nmax,False,plotname,plotdir,plotting,printing)
            f.create_dataset(op.join(ensemble,p0,"pi/E"), data=[E[0].mean])
            f.create_dataset(op.join(ensemble,p0,"pi/Delta_E"), data=[E[0].sdev])
            f.create_dataset(op.join(ensemble,p0,"pi/a"), data=[a[0].mean])
            f.create_dataset(op.join(ensemble,p0,"pi/Delta_a"), data=[a[0].sdev])
            f.create_dataset(op.join(ensemble,p0,"pi/chi2"), data=chi2)
            f.create_dataset(op.join(ensemble,p0,"pi/dof"), data=dof)
            f.create_dataset(op.join(ensemble,p0,"pi/Nmax"), data=Nmax)
            f.create_dataset(op.join(ensemble,p0,"pi/tmin1"), data=1)
            f.create_dataset(op.join(ensemble,p0,"pi/tmax1"), data=T//2-1)

        if op.join(ensemble,p0,"T1/E") not in f.keys():
            C_T1 = get_hdf5_value(fid,op.join(ensemble,p0,"T1/C"))
            cov_T1 = get_hdf5_value(fid,op.join(ensemble,p0,"T1/cov_C"))
            var_T1 = gv.gvar(C_T1[:],cov_T1[:,:])
            cor_T1 = dict(Gab=var_T1/var_T1[0])
            E, a, chi2, dof = fit_correlator_without_bootstrap(cor_T1,T,1,T//2-1,Nmax,False,plotname,plotdir,plotting,printing)
            f.create_dataset(op.join(ensemble,p0,"T1/E"), data=[E[0].mean])
            f.create_dataset(op.join(ensemble,p0,"T1/Delta_E"), data=[E[0].sdev])
            f.create_dataset(op.join(ensemble,p0,"T1/a"), data=[a[0].mean])
            f.create_dataset(op.join(ensemble,p0,"T1/Delta_a"), data=[a[0].sdev])
            f.create_dataset(op.join(ensemble,p0,"T1/chi2"), data=chi2)
            f.create_dataset(op.join(ensemble,p0,"T1/dof"), data=dof)
            f.create_dataset(op.join(ensemble,p0,"T1/Nmax"), data=Nmax)
            f.create_dataset(op.join(ensemble,p0,"T1/tmin1"), data=1)
            f.create_dataset(op.join(ensemble,p0,"T1/tmax1"), data=T//2-1)


        # Fit rho correlator in E irrep (if it exists)
        if "E" in fid[op.join(ensemble,p)].keys():
            C_E = get_hdf5_value(fid,op.join(ensemble,p,"E/C"))
            cov_E = get_hdf5_value(fid,op.join(ensemble,p,"E/cov_C"))
            var_E = gv.gvar(C_E,cov_E)
            cor_E = dict(Gab=var_E/var_E[0])
            E, a, chi2, dof = fit_correlator_without_bootstrap(cor_E,T,1,T//2-1,Nmax,False,plotname,plotdir,plotting,printing)
            f.create_dataset(op.join(ensemble,p,"E/E"), data=[E[0].mean])
            f.create_dataset(op.join(ensemble,p,"E/Delta_E"), data=[E[0].sdev])
            f.create_dataset(op.join(ensemble,p,"E/a"), data=[a[0].mean])
            f.create_dataset(op.join(ensemble,p,"E/Delta_a"), data=[a[0].sdev])
            f.create_dataset(op.join(ensemble,p,"E/chi2"), data=chi2)
            f.create_dataset(op.join(ensemble,p,"E/dof"), data=dof)
            f.create_dataset(op.join(ensemble,p,"E/Nmax"), data=Nmax)
            f.create_dataset(op.join(ensemble,p,"E/tmin1"), data=1)
            f.create_dataset(op.join(ensemble,p,"E/tmax1"), data=T//2-1)


        # Fit rho correlator in B1 irrep (if it exists)
        if "B1" in fid[op.join(ensemble,p)].keys():
            C_B1   = get_hdf5_value(fid,op.join(ensemble,p,"B1/C"))
            cov_B1 = get_hdf5_value(fid,op.join(ensemble,p,"B1/cov_C"))
            var_B1 = gv.gvar(C_B1,cov_B1)
            cor_B1 = dict(Gab=var_B1/var_B1[0])
            E, a, chi2, dof = fit_correlator_without_bootstrap(cor_B1,T,1,T//2-1,Nmax,False,plotname,plotdir,plotting,printing)
            f.create_dataset(op.join(ensemble,p,"B1/E"), data=[E[0].mean])
            f.create_dataset(op.join(ensemble,p,"B1/Delta_E"), data=[E[0].sdev])
            f.create_dataset(op.join(ensemble,p,"B1/a"), data=[a[0].mean])
            f.create_dataset(op.join(ensemble,p,"B1/Delta_a"), data=[a[0].sdev])
            f.create_dataset(op.join(ensemble,p,"B1/chi2"), data=chi2)
            f.create_dataset(op.join(ensemble,p,"B1/dof"), data=dof)
            f.create_dataset(op.join(ensemble,p,"B1/Nmax"), data=Nmax)
            f.create_dataset(op.join(ensemble,p,"B1/tmin1"), data=1)
            f.create_dataset(op.join(ensemble,p,"B1/tmax1"), data=T//2-1)

        if not op.join(ensemble,"lattice") in f.keys():
            f.create_dataset(op.join(ensemble,"lattice"), data=lattice)

        f.close()
    print("Fitting of correlators finished.")

args = sys.argv
infile         = args[1]
outfile        = args[2]
parameterfile  = args[3]

fit_all_files(infile,outfile,parameterfile)