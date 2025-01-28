#pragma once

#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"
#include "MHDSystem.hpp"
#include "linalg/Vector.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

// Should be used only inside a Kokkos::parallel_for
// This is temporary fix until Kokkos a Kokkos::TeamVectorRange with #pragma ivdep inside
template <class iType, class Team>
KOKKOS_INLINE_FUNCTION
static auto TeamVectorRange(const Team& team, const iType start, const iType end) noexcept
    -> decltype(Kokkos::ThreadVectorRange(team, start, end))
{
    const iType chunk = (end - start) / team.team_size();
    const iType start_vector = start + team.team_rank() * chunk;
    const iType end_vector = (team.team_rank()==team.team_size()-1) ? end : start_vector+chunk;

    return Kokkos::ThreadVectorRange(team, start_vector, end_vector);
}

struct BaseKernel
{
    using Euler           = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using VC              = typename Euler::VarCons;
    using VP              = typename Euler::VarPrim;
    static constexpr int nbvar = Euler::nbvar;

    using IntVector  = IntVectorNd<three_d>;

    using Array      = Array3d;
    using ConstArray = ConstArray3d;

    using ConsState  = typename Euler::ConsState;
    using PrimState  = typename Euler::PrimState;

    static Int computeLeagueSize(const IntVector& nbCells, int ghostdepth)
    {
        return (nbCells[IY]+2*ghostdepth)*(nbCells[IZ]+2*ghostdepth);
    }

    template <class kokkos_team>
    KOKKOS_FORCEINLINE_FUNCTION
    static IntVector teamToCoord(const kokkos_team& team, const UniformGrid& grid, int ghostdepth) noexcept
    {
        IntVector coord;

        coord[IZ] = team.league_rank() / (grid.m_nbCells[IY]+2*ghostdepth);
        coord[IY] = team.league_rank() - (grid.m_nbCells[IY]+2*ghostdepth)*coord[IZ];
        coord[IZ] += grid.m_ghostWidths[IZ] - ghostdepth;
        coord[IY] += grid.m_ghostWidths[IY] - ghostdepth;

        coord[IX] = 0;
        return coord;
    }

    template <class kokkos_team>
    KOKKOS_FORCEINLINE_FUNCTION static Int teamToLinearIndex(
            const kokkos_team& team,
            const UniformGrid& grid,
            int ghostdepth) noexcept
    {
        return grid.coordToIndex(teamToCoord(team, grid, ghostdepth));
    }

    KOKKOS_INLINE_FUNCTION
    void copy(const Array& array, Int j, Int j0) const
    {
        for (int ivar=0; ivar<nbvar; ++ivar)
        {
            array(j, ivar) = array(j0, ivar);
        }
    }

    KOKKOS_INLINE_FUNCTION
    void copy(const Array& dst, Int j, const ConstArray& src, Int j0) const
    {
        for (int ivar=0; ivar<nbvar; ++ivar)
        {
            dst(j, ivar) = src(j0, ivar);
        }
    }

    /// WARNING: If called inside a vectorizable loop, it MUST be called with a ConstArray (not Array)
    KOKKOS_FORCEINLINE_FUNCTION
    static ConsState getCons(const ConstArray& array, Int j) noexcept;

    /// WARNING: If called inside a vectorizable loop, it MUST be called with a ConstArray (not Array)
    KOKKOS_FORCEINLINE_FUNCTION
    static PrimState getPrim(const ConstArray& array, Int j) noexcept;

    /// WARNING: If called inside a vectorizable loop, it MUST be called with a ConstArray (not Array)
    KOKKOS_FORCEINLINE_FUNCTION
    static Vector<three_d, Real> getVelocity(const ConstArray& array, Int j) noexcept;

    /// WARNING: If called inside a vectorizable loop, it MUST be called with a ConstArray (not Array)
    KOKKOS_FORCEINLINE_FUNCTION
    static Vector<three_d, Real> getMagneticField(const ConstArray& array, Int j) noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    static void set(const Array& array, Int j, const PrimState& q_j) noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    static void set(const Array& array, Int j, const ConsState& q_j) noexcept;
};


KOKKOS_FORCEINLINE_FUNCTION
BaseKernel::ConsState BaseKernel::getCons(const ConstArray& array, Int j) noexcept
{
    ConsState u_j;
    u_j.d = array(j, VC::ID);
    u_j.e = array(j, VC::IE);
    u_j.m(IX) = array(j, VC::IMx);
    u_j.m(IY) = array(j, VC::IMy);
    u_j.m(IZ) = array(j, VC::IMz);
    u_j.B(IX) = array(j, VC::IBx);
    u_j.B(IY) = array(j, VC::IBy);
    u_j.B(IZ) = array(j, VC::IBz);
    u_j.dX = array(j, VC::I_dX);

    return u_j;
}


KOKKOS_FORCEINLINE_FUNCTION
BaseKernel::PrimState BaseKernel::getPrim(const ConstArray& array, Int j) noexcept
{
    PrimState q_j;
    q_j.d = array(j, VP::ID);
    q_j.p = array(j, VP::IP);
    q_j.v(IX) = array(j, VP::IUx);
    q_j.v(IY) = array(j, VP::IUy);
    q_j.v(IZ) = array(j, VP::IUz);
    q_j.B(IX) = array(j, VP::IBx);
    q_j.B(IY) = array(j, VP::IBy);
    q_j.B(IZ) = array(j, VP::IBz);
    q_j.X = array(j, VP::I_X);

    return q_j;
}


KOKKOS_FORCEINLINE_FUNCTION
Vector<three_d, Real> BaseKernel::getVelocity(const ConstArray& array, Int j) noexcept
{
    Vector<three_d, Real> v_j;
    v_j(IX) = array(j, VP::IUx);
    v_j(IY) = array(j, VP::IUy);
    v_j(IZ) = array(j, VP::IUz);
    return v_j;
}

KOKKOS_FORCEINLINE_FUNCTION
Vector<three_d, Real> BaseKernel::getMagneticField(const ConstArray& array, Int j) noexcept
{
    Vector<three_d, Real> B_j;
    B_j(IX) = array(j, VP::IBx);
    B_j(IY) = array(j, VP::IBy);
    B_j(IZ) = array(j, VP::IBz);
    return B_j;
}

KOKKOS_FORCEINLINE_FUNCTION
void BaseKernel::set(const Array& array, Int j, const PrimState& q_j) noexcept
{
    array(j, VP::ID) = q_j.d;
    array(j, VP::IP) = q_j.p;
    array(j, VP::IUx) = q_j.v(IX);
    array(j, VP::IUy) = q_j.v(IY);
    array(j, VP::IUz) = q_j.v(IZ);
    array(j, VP::IBx) = q_j.B(IX);
    array(j, VP::IBy) = q_j.B(IY);
    array(j, VP::IBz) = q_j.B(IZ);
    array(j, VP::I_X) = q_j.X;

}


KOKKOS_FORCEINLINE_FUNCTION
void BaseKernel::set(const Array& array, Int j, const ConsState& u_j) noexcept
{
    array(j, VC::ID) = u_j.d;
    array(j, VC::IE) = u_j.e;
    array(j, VC::IMx) = u_j.m(IX);
    array(j, VC::IMy) = u_j.m(IY);
    array(j, VC::IMz) = u_j.m(IZ);
    array(j, VC::IBx) = u_j.B(IX);
    array(j, VC::IBy) = u_j.B(IY);
    array(j, VC::IBz) = u_j.B(IZ);
    array(j, VC::I_dX) = u_j.dX;

}

}
