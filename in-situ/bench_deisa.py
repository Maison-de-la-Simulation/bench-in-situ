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
import dask
import dask.array as da
from dask.distributed import performance_report

# Initialize Deisa
scheduler_info = sys.argv[1] if len(sys.argv)>1 else "scheduler.json"
# scheduler_info = "scheduler.json"

deisa = Deisa(scheduler_file_name=scheduler_info, 
              nb_workers=os.environ.get("DASK_NB_WORKERS", 1),
              use_ucx=os.environ.get("DASK_DISTRIBUTED__COMM__UCX__INFINIBAND", False))

print("getting client")
client = deisa.get_client()
# Get client
print("getting deisa array")
arrays = deisa.get_deisa_arrays()

# Select data
gt = arrays["global_t"][:,:,:,:,:]
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

t_stride = 1

slice = arrays["global_t"][0: mt: t_stride, :, z_pos, :, :]

# Check contract
#arrays.check_contract()

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

# output task graph
#sum_over_xy.visualize(filename="sum_over_xy")
#slice.visualize(filename="slice")
#fourier_amplitudes.visualize(filename="fourier_amplitudes")


# s1 = client.persist(ekin_deisa)
s2 = client.persist(sum_over_xy)
s3 = client.persist(slice)
s4 = client.persist(fourier_amplitudes)
# Submit the task graph to the scheduler

# Sign contract
#arrays.validate_contract()

with performance_report(filename="dask-report.html"), dask.config.set(array_optimize=None):
    res2 = client.compute(s2).result()
#    res1 = client.compute(s1).result()
    res3 = client.compute(s3).result()
    res4 = client.compute(s4).result()

 #   print("res1=" + str(res1.sum()))
    print("res2=" + str(res2.sum()))
    print("res3=" + str(res3.sum()))
    print("res4=" + str(res4.sum()))

print("Done ", flush=True)
deisa.wait_for_last_bridge_and_shutdown()
# client.shutdown()
