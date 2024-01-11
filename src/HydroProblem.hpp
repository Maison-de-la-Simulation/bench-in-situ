#pragma once

#include "CopyFromBufferKernel.hpp"
#include "CopyToBufferKernel.hpp"
#include "HydroBoundariesKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#if defined(MPI_SESSION)
#include <MpiCudaAwareSupport.hpp>
#endif

#include <array>
#include <memory>
#include <type_traits>

namespace hydro
{

struct Problem
{
    using IntVector   = IntVectorNd<three_d>;
    using Array       = Array3d;
    using DualArray   = DualArray3d;

    Problem(const std::shared_ptr<Params>& params);
    Problem(const Problem& x) = default;
    Problem(Problem&& x) = default;
    virtual ~Problem() = default;
    Problem& operator=(const Problem& x) = default;
    Problem& operator=(Problem&& x) = default;

    void make_boundaries(Array u, const UniformGrid& grid);

    virtual void initialize(Array u, const UniformGrid& grid) const = 0;
    virtual void make_boundaries_user(Array, const UniformGrid&, int, int);

    std::shared_ptr<Params> m_params;
    Kokkos::Array<int, 2*three_d> m_boundaryTypes;
#if defined(MPI_SESSION)
    distributed_memory_session::MpiCommCart<three_d> comm;
    std::array<DualArray, 2*three_d> buffers_recv;
    std::array<DualArray, 2*three_d> buffers_send;
#endif
};

inline
Problem::Problem(const std::shared_ptr<Params>& params)
    : m_params        {params}
    , m_boundaryTypes {params->mesh.boundaryTypes}
#if defined(MPI_SESSION)
    , comm {}
    , buffers_recv {}
    , buffers_send {}
#endif
{
#if defined(MPI_SESSION)
    std::array<int, three_d> dom2;
    for (int i=0; i<three_d; ++i)
    {
        dom2[i] = static_cast<int>(m_params->mesh.dom[i]);
    }
    comm = distributed_memory_session::MpiCommCart<three_d> {dom2, true, true};

    std::array<int, three_d> CartCoords {comm.getCoords(comm.rank())};
    int ghostWidth = 2;
    for (int idim=0; idim<three_d; ++idim)
    {
        if (CartCoords[idim]!=0)
        {
            m_boundaryTypes[IL+2*idim] = boundary_t::MPI_COMM;
        }
        if (CartCoords[idim]!=static_cast<int>(m_params->mesh.dom[idim]-1))
        {
            m_boundaryTypes[IR+2*idim] = boundary_t::MPI_COMM;
        }
        int size {1};
        for (int idim2=0; idim2<three_d; ++idim2)
        {
            size *= idim2==idim ? ghostWidth : m_params->mesh.nbCells[idim2] + 2*ghostWidth;
        }
        buffers_recv[IL+2*idim] = DualArray("recv"+std::to_string(IL+2*idim), size);
        buffers_recv[IR+2*idim] = DualArray("recv"+std::to_string(IR+2*idim), size);
        buffers_send[IL+2*idim] = DualArray("send"+std::to_string(IL+2*idim), size);
        buffers_send[IR+2*idim] = DualArray("send"+std::to_string(IR+2*idim), size);
    }
#endif
}

inline
void Problem::make_boundaries(Array u, const UniformGrid& grid)
{
#if defined(MPI_SESSION)

    using mem_device = Kokkos::DefaultExecutionSpace::memory_space;
    using mem_host = Kokkos::HostSpace;
    using mem_device_from_mem_host = Kokkos::SpaceAccessibility<mem_host, mem_device>;

#if defined(KOKKOS_ENABLE_CUDA) && defined(Euler_ENABLE_MPI_CUDA_AWARE)
    constexpr bool is_device_cuda = std::is_same<Kokkos::DefaultExecutionSpace,
                                                 Kokkos::Cuda::execution_space>::value;
    constexpr bool is_mpi_device_aware = ((is_device_cuda && isMpiCudaAwareAvailable()) ||
                                          mem_device_from_mem_host::accessible);
#else
    constexpr bool is_mpi_device_aware = mem_device_from_mem_host::accessible;
#endif
    using t_view = std::conditional_t<is_mpi_device_aware,
                                      typename DualArray::t_dev,
                                      typename DualArray::t_host>;

    for (int idim=0; idim<three_d; ++idim)
    {
        auto& buffer_L_send = buffers_send[IL+2*idim];
        auto& buffer_R_send = buffers_send[IR+2*idim];

        auto& buffer_L_recv = buffers_recv[IL+2*idim];
        auto& buffer_R_recv = buffers_recv[IR+2*idim];

        ExecuteCopyToBuffer(u, buffer_L_send.view_device(), grid.curve(), idim, IL);
        ExecuteCopyToBuffer(u, buffer_R_send.view_device(), grid.curve(), idim, IR);
        Kokkos::fence();

        if (!is_mpi_device_aware)
        {
            buffer_L_send.modify_device();
            buffer_R_send.modify_device();
            buffer_L_send.sync_host();
            buffer_R_send.sync_host();
        }

        const int rank {grid.comm.rank()};
        const int min  {grid.comm.shift(static_cast<int>(idim), -1, rank)};
        const int max  {grid.comm.shift(static_cast<int>(idim), +1, rank)};

        const int type {grid.comm.template dataType <Real> ()};
        grid.comm.sendrecv(buffer_L_send.template view<t_view>().data(),
                           static_cast<int>(buffer_L_send.template view<t_view>().span()),
                           type, min, 111+static_cast<int>(idim),
                           buffer_R_recv.template view<t_view>().data(),
                           static_cast<int>(buffer_R_recv.template view<t_view>().span()),
                           type, max, 111+static_cast<int>(idim));

        grid.comm.sendrecv(buffer_R_send.template view<t_view>().data(),
                           static_cast<int>(buffer_R_send.template view<t_view>().span()),
                           type, max, 111+static_cast<int>(idim),
                           buffer_L_recv.template view<t_view>().data(),
                           static_cast<int>(buffer_L_recv.template view<t_view>().span()),
                           type, min, 111+static_cast<int>(idim));

        if (!is_mpi_device_aware)
        {
            buffer_L_recv.modify_host();
            buffer_R_recv.modify_host();
            buffer_L_recv.sync_device();
            buffer_R_recv.sync_device();
        }

        if (m_boundaryTypes[IL+2*idim] == boundary_t::MPI_COMM ||
            m_boundaryTypes[IL+2*idim] == boundary_t::PERIODIC)
        {
            ExecuteCopyFromBuffer(u, buffer_L_recv.view_device(), grid.curve(), idim, IL);
        }
        else
        {
            if (m_boundaryTypes[IL+2*idim] == boundary_t::USER_DEFINED)
            {
                this->make_boundaries_user(u, grid, idim, IL);
            }
            else
            {
                ExecuteBoundaries(u, grid, idim, IL, m_boundaryTypes[IL+2*idim]);
            }
        }
        if (m_boundaryTypes[IR+2*idim] == boundary_t::MPI_COMM ||
            m_boundaryTypes[IR+2*idim] == boundary_t::PERIODIC)
        {
            ExecuteCopyFromBuffer(u, buffer_R_recv.view_device(), grid.curve(), idim, IR);
        }
        else
        {
            if (m_boundaryTypes[IR+2*idim] == boundary_t::USER_DEFINED)
            {
                this->make_boundaries_user(u, grid, idim, IR);
            }
            else
            {
                ExecuteBoundaries(u, grid, idim, IR, m_boundaryTypes[IR+2*idim]);
            }
        }

        Kokkos::fence();
        grid.comm.synchronize();
    }
#else
    for (int idim=0; idim<three_d; ++idim)
    {
        for (int iside=0; iside<2; ++iside)
        {
            if (m_boundaryTypes[iside+2*idim] == boundary_t::USER_DEFINED)
            {
                this->make_boundaries_user(u, grid, idim, iside);
            }
            else
            {
                ExecuteBoundaries(u, grid, idim, iside, m_boundaryTypes[iside+2*idim]);
            }
        }
    }
#endif
}

inline
void Problem::make_boundaries_user(Array, const UniformGrid&, int, int)
{
    throw std::runtime_error("You haven't defined a user boundary condition.");
}

}
