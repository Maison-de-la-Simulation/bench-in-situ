###################################################################################################
# Copyright (c) 2020-2022 Centre national de la recherche scientifique (CNRS)
# Copyright (c) 2020-2022 Commissariat a l'énergie atomique et aux énergies alternatives (CEA)
# Copyright (c) 2020-2022 Institut national de recherche en informatique et en automatique (Inria)
# Copyright (c) 2020-2022 Université Paris-Saclay
# Copyright (c) 2020-2022 Université de Versailles Saint-Quentin-en-Yvelines
#
# SPDX-License-Identifier: MIT
#
###################################################################################################

from deisa import Deisa
import os
import sys
import h5py
import dask.array as da
import yaml

# Initialize Deisa

scheduler_info = sys.argv[1] if len(sys.argv)>1 else "scheduler.json"
# scheduler_info = "scheduler.json"

nb_workers = 1
deisa = Deisa(scheduler_info, nb_workers, use_ucx=os.environ.get("DASK_DISTRIBUTED__COMM__UCX__INFINIBAND", False))

print("getting client")
client = deisa.get_client()
# Get client
print("getting deisa array")
arrays = deisa.get_deisa_arrays()

# Select data
gt = arrays["global_t"][...]
mx = len(gt[0,0,0,0,:])
my = len(gt[0,0,0,:,0])
mz = len(gt[0,0,:,0,0])

mt = len(gt[:,0,0,0,0])

assert(isinstance(mx, int))
assert(isinstance(my, int))
assert(isinstance(mz, int))
print("X-dim =", mx, flush=True)
print("Y-dim =", my, flush=True)
print("Z-dim =", mz, flush=True)
z_pos = int(mz/3)
print("getting slice at z =", z_pos, flush=True)

t_stride = 100

slice = arrays["global_t"][0: mt: t_stride, :, z_pos, :, :]
# gt[time, var, z, y, x]
# slice[time, var, y, x]

# Check contract
arrays.check_contract()

# Construct a lazy task graph
id = 0
iu = 2
iv = 3
iw = 4


ekin_deisa = (
    0.5
    * slice[:, id, :, :]
    * (  slice[:, iu, :, :] * slice[:, iu, :, :]
    + slice[:, iv, :, :] * slice[:, iv, :, :]
    + slice[:, iw, :, :] * slice[:, iw, :, :]
    )
    / (mz * mz)
)
# ekin_deisa[time, y, x]

sum_over_xy = ekin_deisa.sum(axis=(1,2))


ekin_deisa_rechunked = ekin_deisa.rechunk({0: 1, 1: -1, 2: -1}) #no chunking along dim 0, 1, and 2
# npix = ekin_deisa_rechunked.shape[1]
ekin_fft2 = da.fft.fft2(ekin_deisa_rechunked) # fft over the last two axes
fourier_amplitudes = da.absolute(ekin_fft2) **2
# fourier_amplitudes = fourier_amplitudes.reshape(mt/t_stride, mx*my)
# kfreq = da.fft.fftfreq(npix) * npix
# kfreq2D = da.meshgrid(kfreq, kfreq)
# knrm = da.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)
# knrm = knrm.flatten()
# kbins = da.arange(0.5, npix // 2 + 1, 1.0)
# kvals = 0.5 * (kbins[1:] + kbins[:-1])

sum_over_xy = sum_over_xy.to_zarr("sum_over_xy.zarr", overwrite=True, compute=False)
slice = slice.to_zarr("slice.zarr", overwrite=True, compute=False)
fourier_amplitudes = fourier_amplitudes.to_zarr("fourier_amplitudes.zarr", overwrite=True, compute=False)

# s1 = client.persist(ekin_deisa)
s2 = client.persist(sum_over_xy)
s3 = client.persist(slice)
s4 = client.persist(fourier_amplitudes)
# Submit the task graph to the scheduler

# Sign contract
arrays.validate_contract()

# client.compute(s1).result()
client.compute(s2).result()
client.compute(s3).result()
client.compute(s4).result()



#for t in range(0, mt, t_stride):
#    hf = h5py.File("deisa_output_" + str(t) + ".h5", "w")
#    hf.create_dataset("U", data=s3[t/t_stride,iu,:,:])
#    hf.create_dataset("V", data=s3[t/t_stride,iv,:,:])
#    hf.create_dataset("W", data=s3[t/t_stride,iw,:,:])
#    hf.create_dataset("fft2", data=s4[t/t_stride,:])
#    hf.close()

#hf = h5py.File("deisa_sumXY.h5", "w")
#hf.create_dataset("ekin_sum_over_XY", data=s2)
#hf.close()

print("Done ", flush=True)
deisa.wait_for_last_bridge_and_shutdown()
# client.shutdown()
