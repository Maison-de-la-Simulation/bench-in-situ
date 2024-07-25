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

from distributed import performance_report

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

z_pos = 4
print("getting slice at z =", z_pos, flush=True)

t_stride = 100
mt = 500
mz=16

with performance_report(filename="bench_perf_report.html"):
    slice = arrays["global_t"][0: mt: t_stride, :, z_pos, :, :]


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

    sum_over_xy = ekin_deisa.sum(axis=(1,2))


    ekin_deisa_rechunked = ekin_deisa.rechunk({0: 1, 1: -1, 2: -1}) #no chunking along dim 0, 1, and 2
    # npix = ekin_deisa_rechunked.shape[1]
    ekin_fft2 = da.fft.fft2(ekin_deisa_rechunked) # fft over the last two axes
    fourier_amplitudes = da.absolute(ekin_fft2) **2


    sum_over_xy = sum_over_xy.to_zarr("sum_over_xy.zarr", overwrite=True, compute=False)
    slice = slice.to_zarr("slice.zarr", overwrite=True, compute=False)
    fourier_amplitudes = fourier_amplitudes.to_zarr("fourier_amplitudes.zarr", overwrite=True, compute=False)

    print(sum_over_xy)
    print(slice)
    print("===============================")
    
    # Sign contract
    arrays.validate_contract()
    # s1,s2,s3 = client.compute(sum_over_xy, slice, fourier_amplitudes)
    sum_over_xy.compute()
    slice.compute()
    fourier_amplitudes.compute()



print("Done ", flush=True)
deisa.wait_for_last_bridge_and_shutdown()
# client.shutdown()
