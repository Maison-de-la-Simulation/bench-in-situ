#pragma once

#include "HydroConstants.hpp"
#include "HydroFillingCurve.hpp"
#include "HydroTypes.hpp"
#if defined(MPI_SESSION)
#include "MpiCommCart.hpp"
#endif
#include "Print.hpp"

namespace hydro
{
class UniformGrid : public StridedFillingCurve
{
    using Super      = StridedFillingCurve;
    using  IntVector =  IntVectorNd<three_d>;
    using RealVector = RealVector3d;

public:
    UniformGrid(const RealVector& low, const RealVector& up, const IntVector& nbCells,
                const IntVector& dom, const IntVector& ghostWidths);
    UniformGrid(const RealVector& low, const RealVector& up, const IntVector& nbCells,
                const IntVector& dom, Int ghostWidth);
    UniformGrid(const UniformGrid& x) = default;
    UniformGrid(UniformGrid&& x) = default;
    ~UniformGrid() = default;
    UniformGrid& operator=(const UniformGrid& x) = default;
    UniformGrid& operator=(UniformGrid&& x) = default;

    StridedFillingCurve& curve();

    StridedFillingCurve curve() const;

    KOKKOS_INLINE_FUNCTION
    RealVector getCellCenter(const IntVector& coords) const;

    KOKKOS_INLINE_FUNCTION
    RealVector getCellCenter(Int j) const;

    KOKKOS_INLINE_FUNCTION
    Real getdSdV(Int j, int iface) const;

    template <int idir, int iside>
    KOKKOS_INLINE_FUNCTION
    Real getdSdV(Int,
                 std::integral_constant<int, idir>,
                 std::integral_constant<int, iside>) const;

    KOKKOS_INLINE_FUNCTION
    Real dl(int idim) const;

    KOKKOS_INLINE_FUNCTION
    Real lo(int idim) const;

    KOKKOS_INLINE_FUNCTION
    Real hi(int idim) const;

    RealVector m_upGlobal;
    RealVector m_lowGlobal;
    RealVector m_up;
    RealVector m_low;
    IntVector m_dom;
    RealVector m_dl;
    RealVector m_invdl;
#if defined(MPI_SESSION)
    distributed_memory_session::MpiCommCart<three_d> comm;
#endif
};


inline
UniformGrid::UniformGrid(const RealVector& lowGlobal, const RealVector& upGlobal, const IntVector& nbCells,
                              const IntVector& dom, const IntVector& ghostWidths)
    : Super        {nbCells, ghostWidths}
    , m_upGlobal   {upGlobal}
    , m_lowGlobal  {lowGlobal}
    , m_up         {upGlobal}
    , m_low        {lowGlobal}
    , m_dom        {dom}
    , m_dl         {}
    , m_invdl      {}
#if defined(MPI_SESSION)
    , comm {}
#endif
{
#if defined(MPI_SESSION)
    std::array<int, three_d> dom2;
    for (int i=0; i<three_d; ++i)
    {
        dom2[i] = static_cast<int>(dom[i]);
    }
    comm = distributed_memory_session::MpiCommCart<three_d> {dom2, true, true};
    std::array<int, three_d> CartCoords {comm.getCoords(comm.rank())};
    for (int idim=0; idim<three_d; ++idim)
    {
        const Real length_loc {(m_up[idim] - m_low[idim]) / static_cast<Real>(dom[idim])};
        m_low[idim] = m_lowGlobal[idim] + length_loc * static_cast<Real>(CartCoords[idim]);
        m_up[idim] = m_low[idim] + length_loc;
    }
#endif
    for (int idim=0; idim<three_d; ++idim)
    {
        m_dl[idim]  = (m_up[idim] - m_low[idim]) / static_cast<Real>(nbCells[idim]);
        m_invdl[idim] = constants::one / m_dl[idim];
    }
}


inline
UniformGrid::UniformGrid(const RealVector& lowGlobal, const RealVector& upGlobal, const IntVector& nbCells,
                              const IntVector& dom, Int ghostWidth)
    : Super        {nbCells, ghostWidth}
    , m_upGlobal   {upGlobal}
    , m_lowGlobal  {lowGlobal}
    , m_up         {upGlobal}
    , m_low        {lowGlobal}
    , m_dom        {dom}
    , m_dl         {}
    , m_invdl      {}
#if defined(MPI_SESSION)
    , comm {}
#endif
{
#if defined(MPI_SESSION)
    std::array<int, three_d> dom2;
    for (int i=0; i<three_d; ++i)
    {
        dom2[i] = static_cast<int>(dom[i]);
    }
    comm = distributed_memory_session::MpiCommCart<three_d> {dom2, true, true};
    std::array<int, three_d> CartCoords {comm.getCoords(comm.rank())};
    for (int idim=0; idim<three_d; ++idim)
    {
        const Real length_loc {(m_up[idim] - m_low[idim]) / static_cast<Real>(dom[idim])};
        m_low[idim] = m_lowGlobal[idim] + length_loc * static_cast<Real>(CartCoords[idim]);
        m_up[idim] = m_low[idim] + length_loc;
    }
#endif
    for (int idim=0; idim<three_d; ++idim)
    {
        m_dl[idim]  = (m_up[idim] - m_low[idim]) / static_cast<Real>(nbCells[idim]);
        m_invdl[idim] = constants::one / m_dl[idim];
    }
}


inline
StridedFillingCurve& UniformGrid::curve()
{
    return *this;
}


inline
StridedFillingCurve UniformGrid::curve() const
{
    return *this;
}


KOKKOS_INLINE_FUNCTION
RealVector3d UniformGrid::getCellCenter(const IntVector& coords) const
{
    using constants::half;

    RealVector x {m_low};
    for (int idim=0; idim<three_d; ++idim)
    {
        x[idim] += (half+static_cast<Real>(coords[idim]-Super::m_ghostWidths[idim]))*m_dl[idim];
    }
    return x;
}


KOKKOS_INLINE_FUNCTION
RealVector3d UniformGrid::getCellCenter(Int j) const
{
    return getCellCenter(Super::indexToCoord(j));
}


KOKKOS_INLINE_FUNCTION
Real UniformGrid::getdSdV(Int, int iface) const
{
    return m_invdl[iface/2];
}


template <int idir, int iside>
KOKKOS_INLINE_FUNCTION
Real UniformGrid::getdSdV(Int,
                               std::integral_constant<int, idir>,
                               std::integral_constant<int, iside>) const
{
    return m_invdl[idir];
}


KOKKOS_INLINE_FUNCTION
Real UniformGrid::dl(int idim) const
{
    return m_dl[idim];
}


KOKKOS_INLINE_FUNCTION
Real UniformGrid::lo(int idim) const
{
    return m_low[idim];
}


KOKKOS_INLINE_FUNCTION
Real UniformGrid::hi(int idim) const
{
    return m_up[idim];
}

}
