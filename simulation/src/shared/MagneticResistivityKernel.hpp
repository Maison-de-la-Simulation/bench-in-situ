#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroConstants.hpp"
#include "MHDSystem.hpp"
#include "HydroProblem.hpp"
#include "HydroSolver.hpp"
#include "HydroTypes.hpp"
#include "HydroUniformGrid.hpp"

#include <Kokkos_Core.hpp>

namespace hydro
{

class MagneticResistivityKernel : public BaseKernel
{
public:
    using Super = BaseKernel;

    using MHD = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState = typename MHD::ConsState;
    using PrimState = typename MHD::PrimState;
    using RealVector = RealVector3d;

    using IntVector = IntVectorNd< three_d >;
    using TeamPolicy = Kokkos::TeamPolicy< Kokkos::IndexType< Int > >;
    using Team = typename TeamPolicy::member_type;
    static constexpr int ghostDepth = 0;

    MagneticResistivityKernel( const Params& params, const UniformGrid& grid,
                       const Array3d& u, const Array3d& q, Real dt )
        : m_u( u )
        , m_q( q )
        , m_grid( grid )
        , m_dt( dt )
        , m_D( params.hydro.resistivity_coefficient)
        , m_dl  {grid.m_dl}
        , m_splitted_diffusion_step( params.hydro.splitted_diffusion_step)
        , m_conservative_diffusion_update( params.hydro.conservative_diffusion_update)
        , m_x_start(grid.m_ghostWidths[ IX ] - ghostDepth)
        , m_x_end(grid.m_nbCells[ IX ] + grid.m_ghostWidths[ IX ] + ghostDepth)

    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const Int ix, const Int iy, const Int iz ) const
    {
    // Notation : La cellule centrale est q_j puis le reste est sous le format suivant :
    //            q_ _ _ _ les 3 derniers tirets représentent la coordonées i,j,k
    //            On passe à la cellule suivante en remplaçant le tiret avec p (plus) ou m (moins) dans la direction choisie
    //            Si on souhaite aller deux cases dans une direction, on mettra pp ou mm
    // Exemple : q_pp_ = u(i+1,j+1,k) ; q_pp__ = u(i+2,j,k) ; q_p_m = u(i+1,j,k-1) ; q_ppmmp = u(i+2,j-2,k+1)

    const Int j = m_grid.coordToIndex({ix,iy,iz});

    PrimState q_j = Super::getPrim( m_q, j );
    ConsState u_j = Super::getCons( m_u, j );

    Int j_p__ = j + m_grid.m_strides[ IX ];
    PrimState q_p__ = Super::getPrim( m_q, j_p__);

    Int j_m__ = j - m_grid.m_strides[ IX ];
    PrimState q_m__ = Super::getPrim( m_q, j_m__);

    Int j__m_ = j - m_grid.m_strides[ IY ];
    PrimState q__m_ = Super::getPrim( m_q, j__m_);

    Int j__p_ = j + m_grid.m_strides[ IY ];
    PrimState q__p_ = Super::getPrim( m_q, j__p_);

    Int j___m = j - m_grid.m_strides[ IZ ];
    PrimState q___m = Super::getPrim( m_q, j___m);

    Int j___p = j + m_grid.m_strides[ IZ ];
    PrimState q___p = Super::getPrim( m_q, j___p);

    Real Bx_flux_resistivite =  m_D * m_dt * ( (q_p__.B(IX) + q_m__.B(IX) - 2 * q_j.B(IX) ) / (m_dl[IX] * m_dl[IX]) +  (q__p_.B(IX) + q__m_.B(IX) - 2 * q_j.B(IX) ) / (m_dl[IY] * m_dl[IY]) + (q___p.B(IX) + q___m.B(IX) - 2 * q_j.B(IX) ) / (m_dl[IZ] * m_dl[IZ]) );
    Real By_flux_resistivite =  m_D * m_dt * ( (q_p__.B(IY) + q_m__.B(IY) - 2 * q_j.B(IY) ) / (m_dl[IX] * m_dl[IX]) +  (q__p_.B(IY) + q__m_.B(IY) - 2 * q_j.B(IY) ) / (m_dl[IY] * m_dl[IY]) + (q___p.B(IY) + q___m.B(IY) - 2 * q_j.B(IY) ) / (m_dl[IZ] * m_dl[IZ]) );
    Real Bz_flux_resistivite =  m_D * m_dt * ( (q_p__.B(IZ) + q_m__.B(IZ) - 2 * q_j.B(IZ) ) / (m_dl[IX] * m_dl[IX]) +  (q__p_.B(IZ) + q__m_.B(IZ) - 2 * q_j.B(IZ) ) / (m_dl[IY] * m_dl[IY]) + (q___p.B(IZ) + q___m.B(IZ) - 2 * q_j.B(IZ) ) / (m_dl[IZ] * m_dl[IZ]) );

    const Real emag0=0.5*dot(u_j.B, u_j.B);
                                    
    u_j.B(IX) += Bx_flux_resistivite;
    u_j.B(IY) += By_flux_resistivite;
    u_j.B(IZ) += Bz_flux_resistivite;  

    const Real emag1=0.5*dot(u_j.B, u_j.B);

    u_j.e += emag1-emag0;
                              
    Super::set( m_u, j, u_j );
                                  
    }
    

    const Array3d m_u;
    const ConstArray3d m_q;
    const UniformGrid m_grid;
    const Real m_dt;
    const Real m_D;
    const RealVector3d m_dl;
    const bool m_splitted_diffusion_step;
    const bool m_conservative_diffusion_update;

    Int m_x_start;
    Int m_x_end;


};

} // namespace hydro
