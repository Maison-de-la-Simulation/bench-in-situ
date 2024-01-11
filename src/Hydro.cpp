#include "Hydro.hpp"
#include "DistributedMemorySession.hpp"
#include "HydroVersion.hpp"
#include "Print.hpp"

#include <cstdlib>
#include <iostream>
#if defined(MPI_SESSION)
#include <mpi.h>
#include <MpiCudaAwareSupport.hpp>
#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda.h>
#endif // defined(KOKKOS_ENABLE_CUDA)
#endif // defined(MPI_SESSION)
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>

namespace hydro
{

const static std::string code_name {"EulerMultiD"};

#if defined(MPI_SESSION) && defined(KOKKOS_ENABLE_CUDA) && defined(Euler_ENABLE_MPI_CUDA_AWARE)
void intel_opa_patch(std::ostream& os)
{
    // Attempts to handle GPU binding necessary for Intel Omni-path network.
    if(isMpiCudaAwareAvailable())
    {
        auto psm2_cuda_str = std::getenv("PSM2_CUDA");
        if (psm2_cuda_str && std::atoi(psm2_cuda_str) == 1)
        {
            auto local_rank_str = std::getenv("KOKKOS_DEVICE_ID"); // Kokkos variable to select a device
            if (local_rank_str)
            {
                auto device_id = std::atoi(local_rank_str);

                cudaError_t cuda_err = cudaSetDevice(device_id);
                if (cuda_err != cudaSuccess)
                {
                    throw std::runtime_error(cudaGetErrorString(cuda_err));
                }
            }
            else
            {
                // Round-robin binding like in Kokkos
                local_rank_str = std::getenv("OMPI_COMM_WORLD_LOCAL_RANK"); // Open MPI
                if (!local_rank_str) local_rank_str = std::getenv("PMI_RANK"); // Intel MPI
                if (!local_rank_str) local_rank_str = std::getenv("MV2_COMM_WORLD_LOCAL_RANK"); // MVAPICH2
                if (!local_rank_str) local_rank_str = std::getenv("MPI_LOCALRANKID"); // also MVAPICH2
                if (!local_rank_str) local_rank_str = std::getenv("SLURM_LOCALID"); // SLURM
                if (local_rank_str)
                {
                    int count = -1;
                    cudaError_t cuda_err = cudaGetDeviceCount(&count);
                    if (cuda_err != cudaSuccess)
                    {
                        throw std::runtime_error(cudaGetErrorString(cuda_err));
                    }

                    auto local_rank = std::atoi(local_rank_str);
                    auto device_id = local_rank % count;

                    cuda_err = cudaSetDevice(device_id);
                    if (cuda_err != cudaSuccess)
                    {
                        throw std::runtime_error(cudaGetErrorString(cuda_err));
                    }

                    // Set KOKKOS_DEVICE_ID to be coherent inside initialization of Kokkos
                    setenv("KOKKOS_DEVICE_ID", std::to_string(device_id).c_str(), 0);
                    os << "[WARNING] PSM2 in a CUDA-aware build has been detected."
                       << "[WARNING] GPU binding has been done before MPI_Init.\n";
                }
                else
                {
                    throw std::runtime_error("PSM2 in a CUDA-aware build has been detected but no GPU binding could be found"\
                                             ", please set KOKKOS_DEVICE_ID environment variable.");
                }
            }
        }
    }
}
#else
void intel_opa_patch(std::ostream&)
{
}
#endif // defined(MPI_SESSION) && defined(KOKKOS_ENABLE_CUDA) && defined(Euler_ENABLE_MPI_CUDA_AWARE)

void initialize(int& argc, char**& argv)
{
    std::ostringstream oss;
    intel_opa_patch(oss);
    Session::initialize(argc, argv);
    Kokkos::initialize(argc, argv);
    oss << "Initializing " << code_name << '\n';
    Print() << oss.str();
}

void print_configuration(std::ostream& os)
{
    os << std::string(80, '#') << '\n';
    os << "Git information\n";
    os << std::string(80, '#') << '\n';

    os << "Build " << version::git_build_string << '\n';
    os << "Branch " << version::git_branch << '\n';

    os << std::string(80, '#') << '\n';
    os << "Kokkos configuration\n";
    os << std::string(80, '#') << '\n';

    Kokkos::print_configuration(os);

#if defined(MPI_SESSION) && defined(KOKKOS_ENABLE_CUDA)
    const int nproc = Session::getNProc();
    int cudaDeviceId {-1};
    cudaGetDevice(&cudaDeviceId);
    std::vector<int> cudaDeviceIds(nproc);
    ::MPI_Gather(&cudaDeviceId, 1, MPI_INT, cudaDeviceIds.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    os << std::string(80, '#') << '\n';
    os << "GPU pinning\n";
    os << std::string(80, '#') << '\n';
    for (int rank=0; rank<nproc; ++rank)
    {
        os << "MPI rank " << rank << " / " << nproc;
        os << " [MPI_COMM_WORLD],";
        os << " pinned to GPU " << cudaDeviceIds[rank] << ".\n";
    }
    os << std::string(80, '#') << std::endl;
#if defined(Euler_ENABLE_MPI_CUDA_AWARE)
    if (isMpiCudaAwareAvailable())
    {
        os << "MPI CUDA-aware is enabled.\n";
    }
    else
    {
        os << "MPI CUDA-aware is disabled.\n";
    }
#else
    os << "MPI CUDA-aware is disabled.\n";
#endif
    os << std::string(80, '#') << '\n';
#endif // MPI_SESSION && KOKKOS_ENABLE_CUDA
}

void finalize()
{
    Print() << "Finalizing " << code_name << std::endl;
    Kokkos::finalize();
    Session::finalize();
}

}
