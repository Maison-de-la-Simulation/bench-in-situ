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
import time
from distributed.diagnostics import MemorySampler
from pprint import pformat
import matplotlib.pyplot as plt

# Initialize Deisa
if len(sys.argv) < 3:
    raise Exception("Number of dask workers not set. Usage: python3 bench_deisa.py <n_dask_workers> <scheduler_file_name>")
else:
    nb_workers = int(sys.argv[1])
    scheduler_file_name=str(sys.argv[2])
    print(f"parameters: dask workers - {nb_workers}, schedueler_file - {scheduler_file_name}", flush=True)

deisa = Deisa(
    scheduler_file_name=scheduler_file_name,
    nb_workers=nb_workers,
    use_ucx=os.environ.get("DASK_DISTRIBUTED__COMM__UCX__INFINIBAND", False),
)

print("getting client")
client = deisa.get_client()
# Get client
print("getting deisa array")
arrays = deisa.get_deisa_arrays()

# Select data
gt = arrays["global_t"][:, :, :, :, :]
mx = len(gt[0, 0, 0, 0, :])
my = len(gt[0, 0, 0, :, 0])
mz = len(gt[0, 0, :, 0, 0])

mt = len(gt[:, 0, 0, 0, 0])

assert isinstance(mx, int)
assert isinstance(my, int)
assert isinstance(mz, int)
print("X-dim =", mx, flush=True)
print("Y-dim =", my, flush=True)
print("Z-dim =", mz, flush=True)
z_pos = int(mz / 3)
print("getting slice at z =", z_pos, flush=True)

t_stride = 1

slice = arrays["global_t"][0:mt:t_stride, :, z_pos, :, :]

# Check contract
# arrays.check_contract()

# Construct a lazy task graph
id = 0
iu = 2
iv = 3
iw = 4
ms = MemorySampler()

with performance_report(filename="dask-report.html"), dask.config.set(
    array_optimize=None
), ms.sample("collection 1"):
    ekin_deisa = (
        0.5
        * slice[:, id, :, :]
        * (
            slice[:, iu, :, :] * slice[:, iu, :, :]
            + slice[:, iv, :, :] * slice[:, iv, :, :]
            + slice[:, iw, :, :] * slice[:, iw, :, :]
        )
        / (mz * mz)
    )
    # better to lessen the memory used by .persist methods
    # also in general, we do not need to persist in chains of computation, in this specific case,
    # I think we do sum_over_xy deletes ekin_persisted, causing later computations to be recomputed
    # thus making it fail -- also, not super clear why rechunking is needed.
    ekin_persisted = client.persist(ekin_deisa)

    sum_over_xy = ekin_persisted.sum(axis=(1, 2))

    ekin_deisa_rechunked = ekin_persisted.rechunk(
        {0: 1, 1: -1, 2: -1}
    )  # no chunking along dim 0, 1, and 2

    # npix = ekin_deisa_rechunked.shape[1]
    ekin_fft2 = da.fft.fft2(ekin_deisa_rechunked)  # fft over the last two axes
    fourier_amplitudes = da.absolute(ekin_fft2) ** 2
    # fourier_amplitudes = fourier_amplitudes.reshape(mt/t_stride, mx*my)
    # kfreq = da.fft.fftfreq(npix) * npix
    # kfreq2D = da.meshgrid(kfreq, kfreq)
    # knrm = da.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)
    # knrm = knrm.flatten()
    # kbins = da.arange(0.5, npix // 2 + 1, 1.0)
    # kvals = 0.5 * (kbins[1:] + kbins[:-1])

    # output task graph
    # sum_over_xy.visualize(filename="sum_over_xy")
    # slice.visualize(filename="slice")
    # fourier_amplitudes.visualize(filename="fourier_amplitudes")

    ts = time.time()
    res2 = ekin_persisted.compute()
    te = time.time()
    print(f"time ekin: {te-ts}")

    ts = time.time()
    res3 = sum_over_xy.compute()
    te = time.time()
    print(f"time sum over xy: {te-ts}")

    ts = time.time()
    res4 = fourier_amplitudes.compute()
    te = time.time()
    print(f"time fourier: {te-ts}")


# diagnostics info
l1 = client.run(lambda dask_worker: dask_worker.transfer_outgoing_log)
l2 = client.run(lambda dask_worker: dask_worker.transfer_incoming_log)

with open("outgoing.txt", "w") as f1, open("incoming.txt", "w") as f2, open(
    "results.txt", "w"
) as f3:
    f1.write(pformat(l1))
    f2.write(pformat(l2))
    print(f"{res2=}\n{res3=}\n{res4=}", file=f3)

res = ms.plot(align=True)
if isinstance(res, plt.Axes):
    res = res.get_figure()

res.savefig("plot.png")

print("Done ", flush=True)
# deisa.wait_for_last_bridge_and_shutdown()
client.close()
