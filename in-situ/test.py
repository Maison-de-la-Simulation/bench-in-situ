import os

print(">>> LD_LIBRARY_PATH=" + str(os.environ.get("LD_LIBRARY_PATH", None)), flush=True)

# from mpi4py import MPI
from dask import delayed
from dask.distributed import performance_report
from deisa import Deisa

# os.environ["DASK_DISTRIBUTED__COMM__UCX__INFINIBAND"] = "True"

# Scheduler file name and configuration file
scheduler_info = "../build/scheduler.json"
config_file = "../deisa.yml"
# pdi_config_file = '../io.yml'
#
#
# comm = MPI.COMM_WORLD
# size = comm.Get_size()
# rank = comm.Get_rank()
#
# print("MPI size=" + str(size) + ", rank=" + str(rank) + ", Is_initialized=" + str(MPI.Is_initialized()), flush=True)
#
# with open(pdi_config_file, 'r') as f:
#     try:
#         config = yaml.safe_load(f)
#         print(config)
#
#     except yaml.YAMLError as e:
#         print("could not parse file " + pdi_config_file + ". Error=" + str(e))
#         exit(e)
#
# print("<<< pre init", flush=True)
# pdi.init(yaml.dump(config['pdi']))
# print("<<< post init", flush=True)

# Initialize Deisa
deisa = Deisa(scheduler_info, config_file)
print(">>>> Deisa", flush=True)

# Get client
client = deisa.get_client()

# either: Get data descriptor as a list of Deisa arrays object
arrays = deisa.get_deisa_arrays()
print("in-situ: arrays.arrays=", arrays.arrays, flush=True)
print("in-situ: arrays.contract=", arrays.contract, flush=True)
# or: Get data descriptor as a dict of Dask array
# arrays = Deisa.get_dask_arrays()

# # Select data
gt = arrays["global_t"][...]
print(">>>>>>>>>>>>>>> gt=" + str(gt), flush=True)

# Check contract
arrays.check_contract()
print("in-situ: arrays.arrays=", arrays.arrays, flush=True)
print("in-situ: arrays.contract=", arrays.contract, flush=True)


# def Derivee(F, dx):
#     """
#     First Derivative
#        Input: F        = function to be derivate
#               dx       = step of the variable for derivative
#        Output: dFdx = first derivative of F
#     """
#     c0 = 2. / 3.
#     dFdx = c0 / dx * (F[3: - 1] - F[1: - 3] - (F[4:] - F[:- 4]) / 8.)
#     return dFdx


def toggle_high_frequency(gt, ts):
    # HIGH frequency TODO
    data = gt[ts:, 0, 0]
    print("toggle_high_frequency ! ts=" + str(ts), flush=True)
    # deisa.set_return_control({'freq': 'high', 'ts': ts})
    # deisa.set_return_control((True, ts))
    # deisa.set_feedback_data(42)
    # feedback = 1
    # pdi.expose('feedback', feedback, pdi.OUT)

    return data.sum()


def toggle_low_frequency(gt, ts):
    # LOW frequency TODO
    data = gt[ts:, 0, 0]
    print("toggle_low_frequency ! ts=" + str(ts), flush=True)
    # deisa.set_return_control({'freq': 'low', 'ts': ts})
    # deisa.set_return_control((False, ts))
    # feedback = 0
    # pdi.expose('feedback', feedback, pdi.OUT)
    return data.sum()


res1 = client.persist(toggle_low_frequency(gt, 1), release=True)
res100 = client.persist(toggle_high_frequency(gt, 100), release=True)
res200 = client.persist(toggle_low_frequency(gt, 200), release=True)

arrays.validate_contract()
# del gt

# print("res= " + str(client.compute(res).result()), flush=True)
print("res1= " + str(client.compute(res1).result()), flush=True)
print("res100= " + str(client.compute(res100).result()), flush=True)
print("res200= " + str(client.compute(res200).result()), flush=True)

# # py-bokeh is needed if you wanna see the perf report
# with performance_report(filename="dask-report.html"):
#     print("!!!!!!!!!!!!!!!!!!", flush=True)
#
#     d = delayed(hello)(42)
#     client.compute(d)
#
#     # Construct a lazy task graph
#     # cpt = Derivee(gt, 1).mean()
#
#     # Submit the task graph to the scheduler
#     # s = client.persist(cpt, release=True)
#
#     # Sign contract
#     arrays.validate_contract()
#     # del gt
#     # Print the result, note that "s" is a future object, to get the result of the computation,
#     # we call `s.result()` to retreive it.
#     # print(client.compute(s).result(), flush=True)

print("Done", flush=True)
client.shutdown()
# pdi.finalize()
