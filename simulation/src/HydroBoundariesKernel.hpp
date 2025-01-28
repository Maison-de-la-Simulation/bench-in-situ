#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroFillingCurve.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <exception>

namespace hydro
{

struct ReflexiveTag
{
};
struct AdhesiveTag
{
};
struct NeumannTag
{
};
struct PeriodicTag
{
};

class BoundariesKernel : public BaseKernel
{
    using Super = BaseKernel;

    using Euler = MHDSystem;
    using EquationOfState = typename Euler::EquationOfState;
    using ConsState = typename Euler::ConsState;
    using PrimState = typename Euler::PrimState;

    using IntVector = IntVectorNd< three_d>;

public:
    BoundariesKernel( const Array3d& u, const UniformGrid& grid, int idim, int side )
        : m_u( u )
        , m_grid( grid )
        , m_idim( idim )
        , m_side( side )
        , m_boundary_curve()
    {
        IntVector nbCells;
        for ( int idim2 = 0; idim2 < three_d; ++idim2 )
        {
            nbCells[ idim2 ] = idim2 == m_idim
                                   ? m_grid.m_ghostWidths[ idim2 ]
                                   : m_grid.m_nbCells[ idim2 ] + 2 * m_grid.m_ghostWidths[ idim2 ];
        }
        m_boundary_curve = StridedFillingCurve( nbCells, 0 );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReflexiveTag&, Int j ) const
    {
        IntVector coords = m_boundary_curve.indexToCoord( j );
        coords[ m_idim ] +=
            m_side == 0 ? 0 : m_grid.m_ghostWidths[ m_idim ] + m_grid.m_nbCells[ m_idim ];
        IntVector coords0 = coords;
        coords0[ m_idim ] = ( ( m_side == 1 ) ? 2 * m_grid.m_nbCells[ m_idim ] : 0 ) +
                            2 * m_grid.m_ghostWidths[ m_idim ] - 1 - coords[ m_idim ];
        const Int j0 = m_grid.coordToIndex( coords0 );
        ConsState u_j0 = Super::getCons( m_u, j0 );
        u_j0.m( m_idim ) = -u_j0.m( m_idim );

        Super::set( m_u, m_grid.coordToIndex( coords ), u_j0 );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const AdhesiveTag&, Int j ) const
    {
        IntVector coords = m_boundary_curve.indexToCoord( j );
        coords[ m_idim ] +=
            m_side == 0 ? 0 : m_grid.m_ghostWidths[ m_idim ] + m_grid.m_nbCells[ m_idim ];
        IntVector coords0 = coords;
        coords0[ m_idim ] = ( ( m_side == 1 ) ? 2 * m_grid.m_nbCells[ m_idim ] : 0 ) +
                            2 * m_grid.m_ghostWidths[ m_idim ] - 1 - coords[ m_idim ];
        const Int j0 = m_grid.coordToIndex( coords0 );
        ConsState u_j0 = Super::getCons( m_u, j0 );
        for ( int idim = 0; idim < three_d; ++idim )
        {
            u_j0.m( idim ) = -u_j0.m( idim );
        }

        Super::set( m_u, m_grid.coordToIndex( coords ), u_j0 );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const NeumannTag&, Int j ) const
    {
        IntVector coords = m_boundary_curve.indexToCoord( j );
        coords[ m_idim ] +=
            m_side == 0 ? 0 : m_grid.m_ghostWidths[ m_idim ] + m_grid.m_nbCells[ m_idim ];

        IntVector coords0 = coords;
        coords0[ m_idim ] = m_grid.m_ghostWidths[ m_idim ] +
                            ( ( m_side == 1 ) ? m_grid.m_nbCells[ m_idim ] - 1 : 0 );
        const Int j0 = m_grid.coordToIndex( coords0 );
        Super::copy( m_u, m_grid.coordToIndex( coords ), j0 );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const PeriodicTag&, Int j ) const
    {
        IntVector coords = m_boundary_curve.indexToCoord( j );
        coords[ m_idim ] +=
            m_side == 0 ? 0 : m_grid.m_ghostWidths[ m_idim ] + m_grid.m_nbCells[ m_idim ];

        IntVector coords0 = coords;
        coords0[ m_idim ] -= ( m_side == 1 ? 1 : -1 ) * m_grid.m_nbCells[ m_idim ];
        const Int j0 = m_grid.coordToIndex( coords0 );
        Super::copy( m_u, m_grid.coordToIndex( coords ), j0 );
    }

    Array3d m_u;
    UniformGrid m_grid;
    int m_idim;
    int m_side;
    StridedFillingCurve m_boundary_curve;
};

inline void ExecuteBoundaries( const Array3d& u, const UniformGrid& grid, int idim, int side,
                   int boundaryType )
{
    Int nbiter = 1;
    for ( int idim2 = 0; idim2 < three_d; ++idim2 )
    {
        nbiter *= idim2 == idim ? grid.m_ghostWidths[ idim2 ]
                                : grid.m_nbCells[ idim2 ] + 2 * grid.m_ghostWidths[ idim2 ];
    }

    if ( boundaryType == boundary_t::REFLEXIVE )
    {
        Kokkos::RangePolicy< int, ReflexiveTag > range( 0, nbiter );
        BoundariesKernel kernel( u, grid, idim, side );
        Kokkos::parallel_for( "Reflexive boundary kernel - RangePolicy", range, kernel );
    }
    else if ( boundaryType == boundary_t::NEUMANN )
    {
        Kokkos::RangePolicy< int, NeumannTag > range( 0, nbiter );
        BoundariesKernel kernel( u, grid, idim, side );
        Kokkos::parallel_for( "Neumann boundary kernel - RangePolicy", range, kernel );
    }
    else if ( boundaryType == boundary_t::PERIODIC )
    {
        Kokkos::RangePolicy< int, PeriodicTag > range( 0, nbiter );
        BoundariesKernel kernel( u, grid, idim, side );
        Kokkos::parallel_for( "Periodic boundary kernel - RangePolicy", range, kernel );
    }
    else if ( boundaryType == boundary_t::ADHESIVE )
    {
        Kokkos::RangePolicy< int, AdhesiveTag > range( 0, nbiter );
        BoundariesKernel kernel( u, grid, idim, side );
        Kokkos::parallel_for( "Periodic boundary kernel - RangePolicy", range, kernel );
    }
    else
    {
        throw std::runtime_error( "Unknown boundary condition" );
    }
}

} // namespace hydro
