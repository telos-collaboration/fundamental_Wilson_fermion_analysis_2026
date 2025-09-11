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

        ensemble, p, irrep, tmin1, tmax1, tmin2, tmax2, use3x3 = row[0], row[1], row[2], int(row[3]), int(row[4]), int(row[5]), int(row[6]), row[11]
        use3x3 = use3x3 == "true"
        Nmax   = int(row[12])
        fit_id = row[13]

        print(f"{ensemble},{irrep},{p}")
        # read the data from the hdf5 file
        lattice = get_hdf5_value(fid,op.join(ensemble,"lattice"))
        if use3x3 and "eigvals_3x3" in fid[op.join(ensemble,p,irrep,fit_id)].keys():
            ev = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"eigvals_3x3"))
            cov_ev = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"cov_eigvals_3x3"))
            std_ev = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"Delta_eigvals_3x3"))
        else:
            ev = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"eigvals"))
            cov_ev = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"cov_eigvals"))
            std_ev = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"Delta_eigvals"))

        T = ev.shape[0]
        antisymmetric = get_hdf5_value(fid,op.join(ensemble,p,irrep,fit_id,"deriv")) 

        # Rescale data such that eig(t=0)=1 and use full covariance matrix estimator
        var1 = gv.gvar(ev[:,0],cov_ev[:,:,0]/1)
        var2 = gv.gvar(ev[:,1],cov_ev[:,:,1]/1)
        
        eig1 = dict(Gab=var1/var1[0])
        eig2 = dict(Gab=var2/var2[0])
        plotname = op.join(ensemble,irrep)
        plotdir  = "./data_assets/plots/"
        printing = False
        plotting = False

        E1, a1, chi2_1, dof1, fit1, Dfit1, tfit = fit_correlator_without_bootstrap(eig1,T,tmin1,tmax1,Nmax,antisymmetric,plotname,plotdir,plotting,printing)
        E2, a2, chi2_2, dof2, fit2, Dfit2, tfit = fit_correlator_without_bootstrap(eig2,T,tmin2,tmax2,Nmax,antisymmetric,plotname,plotdir,plotting,printing)

        E_A1 = [E1[0].mean,E2[0].mean]
        a_A1 = [a1[0].mean,a2[0].mean]
        Delta_E_A1 = [E1[0].sdev,E2[0].sdev]
        Delta_a_A1 = [a1[0].sdev,a2[0].sdev]
        
        f = h5py.File(outfile, "a")
        if not op.join(ensemble,"lattice") in f.keys():
            f.create_dataset(op.join(ensemble,"lattice"), data=lattice)

        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"tmin1"), data=tmin1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"tmax1"), data=tmax1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"tmin2"), data=tmin2)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"tmax2"), data=tmax2)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"Nmax"), data=Nmax)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"antisymmetric"), data=antisymmetric)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"E"), data=E_A1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"a"), data=a_A1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"Delta_E"), data=Delta_E_A1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"Delta_a"), data=Delta_a_A1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"Delta_a1"), data=[a_i.sdev for a_i in a2])
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"chi2_0"), data=chi2_1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"chi2_1"), data=chi2_2)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"dof0"), data=dof1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"dof1"), data=dof2)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"tfit"), data=tfit)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"fit0"), data=fit1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"fit1"), data=fit2)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"Delta_fit0"), data=Dfit1)
        f.create_dataset(op.join(ensemble,fit_id,p,irrep,"Delta_fit1"), data=Dfit2)
        f.close()
    print("Fitting of correlators finished.")

args = sys.argv
infile         = args[1]
outfile        = args[2]
parameterfile  = args[3]

if op.exists(outfile):
    os.remove(outfile)

fit_all_files(infile,outfile,parameterfile)