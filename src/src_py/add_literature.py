import h5py
import numpy as np

if __name__ == "__main__":
    literature = np.transpose(np.genfromtxt("assets/tables/literature.csv", str, delimiter=";",skip_header=1))
    # print(literature.shape)

    NL = literature[0].astype(int)
    p1 = literature[1].astype(int)
    p2 = literature[2].astype(int)
    p3 = literature[3].astype(int)
    lv = literature[4].astype(int)
    E = literature[5].astype(float)
    Delta_E = literature[6].astype(float)
    irrep = literature[7]
    mpi = literature[8].astype(float)
    mrho = literature[9].astype(float)
    name = literature[10]

    with h5py.File("data_assets/literature.hdf5","w") as hfile:
        for i in range(len(NL)):
            ens = name[i]+str(i)
            lattice_group = hfile.require_group(ens)
            lattice_group.create_dataset("lattice", data=[NL[i],NL[i],NL[i],NL[i]])
            en_lv_group = hfile.require_group(ens+"/p(%i,%i,%i)/%s"%(p1[i],p2[i],p3[i],irrep[i]))
            en_lv_group.create_dataset("E", data=[E[i],])
            en_lv_group.create_dataset("Delta_E", data=[Delta_E[i],])
    
    with open("metadata/scattering_input_literature.csv", "w") as f:
        f.write("ensPirrep;beta;m0;mpi;ld;num_lv\n")
        for i in range(len(NL)):
            pstr = "p(%i,%i,%i)"%(p1[i],p2[i],p3[i])
            f.write("%s%i%s%s;0.0;0.0;%f;%f;False;1\n"%(name[i],i,pstr,irrep[i],mpi[i],mrho[i]))

    done = True

    if done:
        with h5py.File("data_assets/literaure_scattering.hdf5","r") as hfile:
            for ens in hfile:
                for p in hfile[ens]:
                    if p[0] == "p":
                        for irrep in hfile[ens][p]:
                            for lv in hfile[ens][p][irrep]:
                                if lv[:2] == "lv":
                                    E = hfile[ens][p][irrep][lv]["mean"]["aEn"][()]
                                    E_cm_prime = hfile[ens][p][irrep][lv]["mean"]["E_cm_prime"][()]
                                    PS = hfile[ens][p][irrep][lv]["mean"]["PS"][()]
                                    print("%s\nE:\t%f\nEcmp:\t%f\ndelta:\t%f\n\n"%(ens,E,E_cm_prime,PS))




    
        




# "/Lt24Ls14beta6.9m-0.92/p(0,1,1)/A1/E"





# scat_fit_mean.setdefault("irrep",[]).append(irrep)

# mean_group = hfile.require_group(fit_beta_m+"/mean")
# mean_group.create_dataset(key, data=val)