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

os.environ["DASK_DISTRIBUTED__COMM__UCX__INFINIBAND"] = "True"

# Initialize Deisa

# scheduler_info = sys.argv[1] if sys.argv[1] else "scheduler.json"
scheduler_info = "scheduler.json"

nb_workers = 1
deisa = Deisa(scheduler_info, nb_workers)

print("getting client")
client = deisa.get_client()
# Get client
print("getting deisa array")
arrays = deisa.get_deisa_arrays()

# Select data
gt = arrays["global_t"][...]
mz = len(gt[0,0,0,:,0])
assert(isinstance(mz, int))
print("Z-dim = ", mz, flush=True)

print("getting slice")
slice = arrays["global_t"][:, :, :, 5, :]

# Check contract
arrays.check_contract()

# Construct a lazy task graph
id = 0
iu = 2
iv = 3
iw = 4

ekin_deisa = (
    0.5
    * slice[:, :, :, id]
    * (  slice[:, :, :, iu] * slice[:, :, :, iu]
    + slice[:, :, :, iv] * slice[:, :, :, iv]
    + slice[:, :, :, iw] * slice[:, :, :, iw]
    )
    / (mz * mz)
)

sum_over_xy = ekin_deisa.sum(axis=(1,2))



# ekin_deisa_rechunked = ekin_deisa.rechunk(({0: -1, 1: -1, 2: -1})) #no chanking along dim 0 and dim 1
# npix = ekin_deisa_rechunked.shape[1]
# fourier_image = da.fft.fftn(ekin_deisa_rechunked)
# fourier_amplitudes = da.absolute(fourier_image) ** 2
# kfreq = da.fft.fftfreq(npix) * npix
# kfreq2D = da.meshgrid(kfreq, kfreq)
# knrm = da.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)
# knrm = knrm.flatten()
# fourier_amplitudes2 = fourier_amplitudes.flatten()
# kbins = da.arange(0.5, npix // 2 + 1, 1.0)
# kvals = 0.5 * (kbins[1:] + kbins[:-1])

s1 = client.persist(ekin_deisa)
s2 = client.persist(sum_over_xy)
# Submit the task graph to the scheduler

# s1, s2, s3, s4= client.persist([knrm, fourier_amplitudes, kbins, kvals])
# s1, s3, s4 = client.persist([knrm, kbins, kvals])
# s2 = client.persist(fourier_amplitudes)

# Sign contract
arrays.validate_contract()

client.compute(s1).result()
client.compute(s2).result()
# client.compute(s3).result()
# client.compute(s4).result()


# for t in range(100):
#     hf = h5py.File("fft_from_deisa_" + str(t) + ".h5", "w")
#     hf.create_dataset("knrm", data=s2[t])
#     # hf.create_dataset("knrm", data=s1)
#     # hf.create_dataset("fourier_amplitudes", data=fourier_amplitudes)
#     # hf.create_dataset("kbins", data=s3)
#     # hf.create_dataset("kvals", data=s4)

#     hf.close()

print("Done ", flush=True)
deisa.wait_for_last_bridge_and_shutdown()
# client.shutdown()
