import os.path as op
import gvar as gv
import os
import h5py
import csv
import sys
from tqdm import tqdm
from corrfitter_utils import get_hdf5_value, fit_correlator_without_bootstrap

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
            
        E, a, chi2, dof, fit, Dfit = fit_correlator_without_bootstrap(cor,T,tmin,tmax,Nmax,False,plotname,plotdir,plotting,printing)
        
        f.create_dataset(op.join(ensemble,p,irrep,"E"), data=[E[0].mean])
        f.create_dataset(op.join(ensemble,p,irrep,"Delta_E"), data=[E[0].sdev])
        f.create_dataset(op.join(ensemble,p,irrep,"a"), data=[a[0].mean])
        f.create_dataset(op.join(ensemble,p,irrep,"Delta_a"), data=[a[0].sdev])
        f.create_dataset(op.join(ensemble,p,irrep,"chi2"), data=chi2)
        f.create_dataset(op.join(ensemble,p,irrep,"dof"), data=dof)
        f.create_dataset(op.join(ensemble,p,irrep,"Nmax"), data=Nmax)
        f.create_dataset(op.join(ensemble,p,irrep,"tmin"), data=tmin)
        f.create_dataset(op.join(ensemble,p,irrep,"tmax"), data=tmax)
        f.create_dataset(op.join(ensemble,p,irrep,"fit"), data=fit)
        f.create_dataset(op.join(ensemble,p,irrep,"Delta_fit"), data=Dfit)

        if not op.join(ensemble,"lattice") in f.keys():
            f.create_dataset(op.join(ensemble,"lattice"), data=lattice)

        f.close()
    print("Fitting of correlators finished.")

args = sys.argv
infile         = args[1]
outfile        = args[2]
parameterfile  = args[3]

fit_all_files(infile,outfile,parameterfile)