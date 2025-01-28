#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

namespace hydro
{

#if defined(FLOAT64_EXTENDED)
using Real = long double;
#elif defined(FLOAT64)
using Real = double;
#elif defined(FLOAT32)
using Real = float;
#else
static_assert(false, "You must define a floating point type")
#endif

#if defined(INT64)
using  Int = long int;
#elif defined(INT32)
using  Int = int;
#else
static_assert(false, "You must define an integer type")
#endif

enum dim_t : Int
{
    one_d=1,
    two_d=2,
    three_d=3
};

enum dir_t : Int
{
    IX=0,
    IY=1,
    IZ=2
};

enum side_t : int
{
    IL=0,
    IR=1
};

enum boundary_t : int
{
    UNKNOWN = -2,
    MPI_COMM = -1,
    USER_DEFINED = 0,
    REFLEXIVE = 1,
    NEUMANN = 2,
    PERIODIC = 3,
    ADHESIVE = 4
};

inline
boundary_t s2boundary(const std::string& name)
{
    boundary_t type;
    if (name == "mpi_comm")
    {
        type = boundary_t::MPI_COMM;
    }
    else if (name == "user_defined")
    {
        type = boundary_t::USER_DEFINED;
    }
    else if (name == "reflexive")
    {
        type = boundary_t::REFLEXIVE;
    }
    else if (name == "neumann")
    {
        type = boundary_t::NEUMANN;
    }
    else if (name == "periodic")
    {
        type = boundary_t::PERIODIC;
    }
    else if (name == "adhesive")
    {
        type = boundary_t::ADHESIVE;
    }
    else
    {
        type = boundary_t::UNKNOWN;
    }
    return type;
}

struct Tags
{
    template <int i>
    using Tag_t = std::integral_constant<int, i>;

    using Dim1_t = Tag_t<one_d>;
    using Dim2_t = Tag_t<two_d>;
    using Dim3_t = Tag_t<three_d>;
    static constexpr Dim1_t Dim1 = Dim1_t ();
    static constexpr Dim2_t Dim2 = Dim2_t ();
    static constexpr Dim3_t Dim3 = Dim3_t ();

    using DirX_t = Tag_t<IX>;
    using DirY_t = Tag_t<IY>;
    using DirZ_t = Tag_t<IZ>;
    static constexpr DirX_t DirX = DirX_t ();
    static constexpr DirY_t DirY = DirY_t ();
    static constexpr DirZ_t DirZ = DirZ_t ();

    using SideL_t = Tag_t<IL>;
    using SideR_t = Tag_t<IR>;
    static constexpr SideL_t SideL = SideL_t ();
    static constexpr SideR_t SideR = SideR_t ();
};

using Layout = Kokkos::LayoutLeft;

template <Int size>
using IntVectorNd = Kokkos::Array<Int, size>;

using RealVector3d = Kokkos::Array<Real, three_d>;

// Arrays with compile-time dimension
using Array3d = Kokkos::View<Real*[2*three_d+2+1], Layout>;

using DualArray3d = Kokkos::DualView<Real*[2*three_d+2+1], Layout>;

using ConstArray3d = Kokkos::View<const Real*[2*three_d+2+1], Layout>;

using HostArray3d = typename Array3d::HostMirror;

using HostConstArray3d = typename ConstArray3d::HostMirror;

// Arrays with run-time dimension
using ArrayDyn = Kokkos::View<Real**, Layout>;
using ConstArrayDyn = Kokkos::View<const Real**, Layout>;
using HostArrayDyn = ArrayDyn::HostMirror;
using HostConstArrayDyn = ConstArrayDyn::HostMirror;
}
