#pragma once

#include "HydroFillingCurve.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

namespace hydro
{

class CopyFromBufferKernel
{
    using IntVector = IntVectorNd< three_d >;

public:
    CopyFromBufferKernel( const Array3d& u, const Array3d& buffer,
                          const StridedFillingCurve& curve, int idim, int side )
        : m_u( u )
        , m_buffer( buffer )
        , m_nvar( 0 )
        , m_curve( curve )
        , m_idim( idim )
        , m_side( side )
        , m_boundary_curve()
    {
        if ( m_buffer.extent( 1 ) != m_u.extent( 1 ) )
        {
            throw std::runtime_error( "Arrays don't have same size.\n" );
        }
        else
        {
            m_nvar = static_cast< int >( m_buffer.extent( 1 ) );
        }

        IntVector nbCells;
        for ( int idim2 = 0; idim2 < three_d; ++idim2 )
        {
            nbCells[ idim2 ] =
                idim2 == m_idim ? m_curve.m_ghostWidths[ idim2 ]
                                : m_curve.m_nbCells[ idim2 ] + 2 * m_curve.m_ghostWidths[ idim2 ];
        }
        m_boundary_curve = StridedFillingCurve( nbCells, 0 );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( Int j ) const
    {
        IntVector coords = m_boundary_curve.indexToCoord( j );
        coords[ m_idim ] +=
            ( m_side == 1 ? m_curve.m_nbCells[ m_idim ] + m_curve.m_ghostWidths[ m_idim ] : 0 );
        const Int j0 = m_curve.coordToIndex( coords );
        for ( int ivar = 0; ivar < m_nvar; ++ivar )
        {
            m_u( j0, ivar ) = m_buffer( j, ivar );
        }
    }

    Array3d m_u;
    ConstArray3d m_buffer;
    int m_nvar;
    StridedFillingCurve m_curve;
    int m_idim;
    int m_side;
    StridedFillingCurve m_boundary_curve;
};

inline
void ExecuteCopyFromBuffer( const Array3d& u, const Array3d& buffer,
                       const StridedFillingCurve& curve, int idim, int side )
{
    CopyFromBufferKernel kernel( u, buffer, curve, idim, side );
    Kokkos::RangePolicy< Int > range( 0, kernel.m_boundary_curve.nbCells() );
    Kokkos::parallel_for( "Copy from buffer kernel - RangePolicy", range, kernel );
}

} // namespace hydro
