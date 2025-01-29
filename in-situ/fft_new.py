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
from dask.distributed import performance_report

import os
import sys
#import h5py
import dask
import dask.array as da

os.environ["DASK_DISTRIBUTED__COMM__UCX__INFINIBAND"] = "True"

# Initialize Deisa

scheduler_info = sys.argv[1] if sys.argv[1] else "scheduler.json"
nb_workers = 1
deisa = Deisa(scheduler_file_name=scheduler_info, nb_workers=nb_workers)

# Get client
client = deisa.get_client()
arrays = deisa.get_deisa_arrays()

# Select data
#iz_middle_gloc = mz * nz / 2 - 1
t = 3  # Temps auquel on veut le sdp

global_t = arrays["global_t"][:,:,:,:,:]
slice = global_t[t, 0, :, :, :]
print("global_t=" + str(global_t), flush=True)
print("slice=" + str(slice), flush=True)

#nz = 16
nz = len(slice[:,0,0])
assert(isinstance(nz, int))
#mz = 1
mz = len(global_t[0,:,0,0,0])
assert(isinstance(mz, int))
prefix = "deisa"
num_restart = 0

print("nz=" + str(nz))
print("mz=" + str(mz))

# Check contract
#arrays.check_contract()

# Construct a lazy task graph
id = 0
iu = 2
iv = 3
iw = 4

ekin_deisa = (
    0.5
    * slice[:, :, id]
    * (
        slice[:, :, iu] * slice[:, :, iu]
        + slice[:, :, iv] * slice[:, :, iv]
        + slice[:, :, iw] * slice[:, :, iw]
    )
    / (mz * nz * mz * nz)
)

ekin_deisa_rechunked = ekin_deisa.rechunk({0: -1, 1: -1})

npix = ekin_deisa_rechunked.shape[0]
print("npix=", npix)
fourier_image = da.fft.fftn(ekin_deisa_rechunked)
fourier_amplitudes = da.absolute(fourier_image) ** 2
kfreq = da.fft.fftfreq(npix) * npix
kfreq2D = da.meshgrid(kfreq, kfreq)
knrm = da.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)
knrm = knrm.flatten()
fourier_amplitudes = fourier_amplitudes.flatten()
kbins = da.arange(0.5, npix // 2 + 1, 1.0)
kvals = 0.5 * (kbins[1:] + kbins[:-1])


# Submit the task graph to the scheduler

s1, s2, s3, s4 = client.persist([knrm, fourier_amplitudes, kbins, kvals])
# Sign contract
#arrays.validate_contract()

with performance_report(filename="dask-report.html"), dask.config.set(array_optimize=None):
    res2 = client.compute(s2).result()
    res1 = client.compute(s1).result()
    res3 = client.compute(s3).result()
    res4 = client.compute(s4).result()

    print("res1=" + str(res1.sum()))
    print("res2=" + str(res2.sum()))
    print("res3=" + str(res3))
    print("res4=" + str(res4))


#hf = h5py.File("fft_from_deisa_" + prefix + "_" + str(num_restart) + ".h5", "w")
#hf.create_dataset("knrm", data=knrm)
#hf.create_dataset("fourier_amplitudes", data=fourier_amplitudes)
#hf.create_dataset("kbins", data=kbins)
#hf.create_dataset("kvals", data=kvals)
#hf.close()

print("Done", flush=True)
deisa.wait_for_last_bridge_and_shutdown()
