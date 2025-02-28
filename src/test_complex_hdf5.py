import h5py
import numpy as np

# zz = complex(123, 3223)
# zznp = np.complex(123, 3223)


# print(type(zz))
# print(type(zznp))


# with h5py.File("test_complex.hdf5", "w") as f:
#     f.create_dataset("complex_numer", data=zz)

# with h5py.File("test_complex.hdf5", "r") as f:
#     print(f["complex_numer"][()])
#     print(type(f["complex_numer"][()]))

NT = [32,32,24,24,32,32,32,32,32,32]
NL = [16,16,14,14,24,24,16,16,24,24]
ds = [[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[1,1,0],[1,1,0],[1,1,0],[1,1,0]]
lvs = [1,2,1,2,1,2,1,2,1,2]
mpi = 0.38649
mrho = 0.5494

with h5py.File("data/fitresults_Feb26.hdf5", "r") as f:
    f.visit(print)
    with open("data/fitresults_Feb26.dat","w") as out:
        for i in range(len(NT)):
            hdf5str = "Lt%iLs%ibeta6.9m1-0.92m2-0.92/p(%i,%i,%i)/"%(NT[i],NL[i],ds[i][0],ds[i][1],ds[i][2])
            out.write("%i %i %i %i %i %f %f %f %f %f\n"%(NL[i],ds[i][0],ds[i][1],ds[i][2],lvs[i],f[hdf5str+"E%i"%(lvs[i]-1)][()][0],f[hdf5str+"Delta_E%i"%(lvs[i]-1)][()][0],f[hdf5str+"Delta_E%i"%(lvs[i]-1)][()][0],mpi,mrho))
        






