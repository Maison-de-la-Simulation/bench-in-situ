import zarr
import h5py


z2 = zarr.open('sum_over_xy.zarr', mode='r')
print(z2)

hf = h5py.File("deisa_sumXY.h5", "w")
hf.create_dataset("ekin_sum_over_XY", data=z2)
hf.close()

z2 = zarr.open('slice.zarr', mode='r')
print(z2)

z3 = z2[:,0,:,:]
print(z3)

hf = h5py.File("deisa_slices.h5", "w")
hf.create_dataset("slices", data=z3)
hf.close()



z2 = zarr.open('fourier_amplitudes.zarr', mode='r')
print(z2)

hf = h5py.File("deisa_fourier.h5", "w")
hf.create_dataset("slices", data=z2)
hf.close()