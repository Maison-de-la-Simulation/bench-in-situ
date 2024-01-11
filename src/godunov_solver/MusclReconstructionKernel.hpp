#pragma once

#include "HydroBaseKernel.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"

namespace hydro
{

class MusclReconstructionKernel : public BaseKernel
{
    using Super           = BaseKernel;

public:
    using MHD           = MHDSystem;
    using EquationOfState = typename MHD::EquationOfState;
    using ConsState           = typename MHD::ConsState;
    using PrimState           = typename MHD::PrimState;

    using RealVector = RealVector3d;
    using VC              = typename MHD::VarCons;
    using VP              = typename MHD::VarPrim;
    static constexpr int nbvar = MHD::nbvar;

    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::IndexType<Int>>;
    using Team       = typename TeamPolicy::member_type;
    static constexpr int ghostDepth = 1;

    MusclReconstructionKernel(const Params& params, const UniformGrid& grid,
                              Array3d q, const Kokkos::Array<Array3d, 2*three_d>& qr, Real dt)
        : m_gamma(params.thermo.gamma)
        , m_qr(qr)
        , m_q(q)
        , m_grid(grid)
        , m_dtdl()
        , m_slope_type(params.run.slope_type)
    {
        for (int idim=0; idim<three_d; ++idim)
        {
            m_dtdl[idim] = dt / grid.m_dl[idim];
        }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const Int ix, const Int iy, const Int iz ) const
    {
      const Int j = m_grid.coordToIndex({ix,iy,iz});
      const PrimState q_j = Super::getPrim(m_q, j);

      Kokkos::Array<PrimState, three_d> dq_j;

      PrimState qMinus = Super::getPrim(m_q, j - m_grid.m_strides[IX]);
      PrimState qPlus = Super::getPrim(m_q, j + m_grid.m_strides[IX]);
      // get hydro slopes dq
      dq_j[IX] = slope_unsplit_hydro_1d(q_j, qPlus, qMinus);

      qMinus= Super::getPrim(m_q, j - m_grid.m_strides[IY]);
      qPlus = Super::getPrim(m_q, j + m_grid.m_strides[IY]);
      // get hydro slopes dq
      dq_j[IY] = slope_unsplit_hydro_1d(q_j, qPlus, qMinus);

      qMinus = Super::getPrim(m_q, j - m_grid.m_strides[IZ]);
      qPlus  = Super::getPrim(m_q, j + m_grid.m_strides[IZ]);
      // get hydro slopes dq
      dq_j[IZ] = slope_unsplit_hydro_1d(q_j, qPlus, qMinus);

      //get un+1/2 with matrix
      const PrimState sq0_j = trace_unsplit_hydro_along_dir(q_j, dq_j);

      // get edge values
      PrimState qm_j;
      PrimState qp_j;

      qm_j.d = q_j.d - 0.5*dq_j[IX].d + 0.5*m_dtdl[IX]*sq0_j.d;
      qp_j.d = q_j.d + 0.5*dq_j[IX].d + 0.5*m_dtdl[IX]*sq0_j.d;
      qm_j.v = q_j.v - 0.5*dq_j[IX].v + 0.5*m_dtdl[IX]*sq0_j.v;
      qp_j.v = q_j.v + 0.5*dq_j[IX].v + 0.5*m_dtdl[IX]*sq0_j.v;
      qm_j.p = q_j.p - 0.5*dq_j[IX].p + 0.5*m_dtdl[IX]*sq0_j.p;
      qp_j.p = q_j.p + 0.5*dq_j[IX].p + 0.5*m_dtdl[IX]*sq0_j.p;
      qm_j.B = q_j.B - 0.5*dq_j[IX].B + 0.5*m_dtdl[IX]*sq0_j.B;
      qp_j.B = q_j.B + 0.5*dq_j[IX].B + 0.5*m_dtdl[IX]*sq0_j.B;

      qm_j.X = q_j.X - 0.5*dq_j[IX].X + 0.5*m_dtdl[IX]*sq0_j.X;
      qp_j.X = q_j.X + 0.5*dq_j[IX].X + 0.5*m_dtdl[IX]*sq0_j.X;

      Super::set(m_qr[0 + 2*IX], j, qm_j);
      Super::set(m_qr[1 + 2*IX], j, qp_j);

      qm_j.d = q_j.d - 0.5*dq_j[IY].d + 0.5*m_dtdl[IY]*sq0_j.d;
      qp_j.d = q_j.d + 0.5*dq_j[IY].d + 0.5*m_dtdl[IY]*sq0_j.d;
      qm_j.v = q_j.v - 0.5*dq_j[IY].v + 0.5*m_dtdl[IY]*sq0_j.v;
      qp_j.v = q_j.v + 0.5*dq_j[IY].v + 0.5*m_dtdl[IY]*sq0_j.v;
      qm_j.p = q_j.p - 0.5*dq_j[IY].p + 0.5*m_dtdl[IY]*sq0_j.p;
      qp_j.p = q_j.p + 0.5*dq_j[IY].p + 0.5*m_dtdl[IY]*sq0_j.p;
      qm_j.B = q_j.B - 0.5*dq_j[IY].B + 0.5*m_dtdl[IY]*sq0_j.B;
      qp_j.B = q_j.B + 0.5*dq_j[IY].B + 0.5*m_dtdl[IY]*sq0_j.B;

      qm_j.X = q_j.X - 0.5*dq_j[IY].X + 0.5*m_dtdl[IY]*sq0_j.X;
      qp_j.X = q_j.X + 0.5*dq_j[IY].X + 0.5*m_dtdl[IY]*sq0_j.X;

      Super::set(m_qr[0 + 2*IY], j, qm_j);
      Super::set(m_qr[1 + 2*IY], j, qp_j);

      qm_j.d = q_j.d - 0.5*dq_j[IZ].d + 0.5*m_dtdl[IZ]*sq0_j.d;
      qp_j.d = q_j.d + 0.5*dq_j[IZ].d + 0.5*m_dtdl[IZ]*sq0_j.d;
      qm_j.v = q_j.v - 0.5*dq_j[IZ].v + 0.5*m_dtdl[IZ]*sq0_j.v;
      qp_j.v = q_j.v + 0.5*dq_j[IZ].v + 0.5*m_dtdl[IZ]*sq0_j.v;
      qm_j.p = q_j.p - 0.5*dq_j[IZ].p + 0.5*m_dtdl[IZ]*sq0_j.p;
      qp_j.p = q_j.p + 0.5*dq_j[IZ].p + 0.5*m_dtdl[IZ]*sq0_j.p;
      qm_j.B = q_j.B - 0.5*dq_j[IZ].B + 0.5*m_dtdl[IZ]*sq0_j.B;
      qp_j.B = q_j.B + 0.5*dq_j[IZ].B + 0.5*m_dtdl[IZ]*sq0_j.B;

      qm_j.X = q_j.X - 0.5*dq_j[IZ].X + 0.5*m_dtdl[IZ]*sq0_j.X;
      qp_j.X = q_j.X + 0.5*dq_j[IZ].X + 0.5*m_dtdl[IZ]*sq0_j.X;

      Super::set(m_qr[0 + 2*IZ], j, qm_j);
      Super::set(m_qr[1 + 2*IZ], j, qp_j);
    };



    KOKKOS_INLINE_FUNCTION
    PrimState trace_unsplit_hydro_along_dir(const PrimState& q, const Kokkos::Array<PrimState, three_d>& dq) const noexcept
    {
        using namespace constants;
        PrimState sq0;

        Real div_u = (dq[IX].v(IX) + dq[IY].v(IY) + dq[IZ].v(IZ));

        Vector<three_d, Real> grad_d =  {dq[IX].d    , dq[IY].d    , dq[IZ].d};
        Vector<three_d, Real> grad_p =  {dq[IX].p    , dq[IY].p    , dq[IZ].p};
        Vector<three_d, Real> grad_vx = {dq[IX].v(IX), dq[IY].v(IX), dq[IZ].v(IX)};
        Vector<three_d, Real> grad_vy = {dq[IX].v(IY), dq[IY].v(IY), dq[IZ].v(IY)};
        Vector<three_d, Real> grad_vz = {dq[IX].v(IZ), dq[IY].v(IZ), dq[IZ].v(IZ)};
        Vector<three_d, Real> grad_Bx = {dq[IX].B(IX), dq[IY].B(IX), dq[IZ].B(IX)};
        Vector<three_d, Real> grad_By = {dq[IX].B(IY), dq[IY].B(IY), dq[IZ].B(IY)};
        Vector<three_d, Real> grad_Bz = {dq[IX].B(IZ), dq[IY].B(IZ), dq[IZ].B(IZ)};
        Vector<three_d, Real> grad_X  = {dq[IX].X,     dq[IY].X,      dq[IZ].X};

        Real r = (1.0/q.d);

        sq0.d=0.0; sq0.p=0.0; sq0.B(IX) =0.0;sq0.B(IY) =0.0;sq0.B(IZ) =0.0; sq0.v(IX) =0.0;sq0.v(IY) =0.0;sq0.v(IZ) =0.0; sq0.X=0.0;

        sq0.d -= q.d*div_u + dot(q.v,grad_d);
        sq0.p -= m_gamma*q.p*div_u + dot(q.v,grad_p);
        //
        sq0.v(IX) -= dot(q.v, grad_vx) + r*(dq[IX].p + dot(dq[IX].B,q.B) - dot(q.B,grad_Bx));
        sq0.v(IY) -= dot(q.v, grad_vy) + r*(dq[IY].p + dot(dq[IY].B,q.B) - dot(q.B,grad_By));
        sq0.v(IZ) -= dot(q.v, grad_vz) + r*(dq[IZ].p + dot(dq[IZ].B,q.B) - dot(q.B,grad_Bz));
        //
        sq0.B(IX) -= + q.B(IX)*div_u + dot(q.v,grad_Bx) - dot(q.B, grad_vx) ;
        sq0.B(IY) -= + q.B(IY)*div_u + dot(q.v,grad_By) - dot(q.B, grad_vy) ;
        sq0.B(IZ) -= + q.B(IZ)*div_u + dot(q.v,grad_Bz) - dot(q.B, grad_vz) ;

        sq0.X -= dot(q.v,grad_X);

        return sq0;
    }

    KOKKOS_INLINE_FUNCTION
    Real slope_unsplit_hydro_1d_scalar(Real const v_l, Real const v_c, Real const v_r) const noexcept
    {
        using namespace constants;

        const Real dlft = m_slope_type*(v_c - v_l);
        const Real drgt = m_slope_type*(v_r - v_c);
        const Real dcen = half * (v_r - v_l);
        const Real dsgn = (dcen > zero) ? one : -one;
        const Real slop = std::fmin(std::fabs(dlft), std::fabs(drgt));
        const Real dlim = (dlft*drgt) < zero ? zero : slop;

        return dsgn * std::fmin(dlim, std::fabs(dcen));
    }

    KOKKOS_INLINE_FUNCTION
    PrimState slope_unsplit_hydro_1d(const PrimState& q,
                                     const PrimState& qPlus,
                                     const PrimState& qMinus) const noexcept
    {
        PrimState dq;

        dq.d = slope_unsplit_hydro_1d_scalar(qMinus.d, q.d, qPlus.d);
        dq.p = slope_unsplit_hydro_1d_scalar(qMinus.p, q.p, qPlus.p);
        dq.v(IX) = slope_unsplit_hydro_1d_scalar(qMinus.v(IX), q.v(IX), qPlus.v(IX));
        dq.v(IY) = slope_unsplit_hydro_1d_scalar(qMinus.v(IY), q.v(IY), qPlus.v(IY));
        dq.v(IZ) = slope_unsplit_hydro_1d_scalar(qMinus.v(IZ), q.v(IZ), qPlus.v(IZ));
        dq.B(IX) = slope_unsplit_hydro_1d_scalar(qMinus.B(IX), q.B(IX), qPlus.B(IX));
        dq.B(IY) = slope_unsplit_hydro_1d_scalar(qMinus.B(IY), q.B(IY), qPlus.B(IY));
        dq.B(IZ) = slope_unsplit_hydro_1d_scalar(qMinus.B(IZ), q.B(IZ), qPlus.B(IZ));
        dq.X = slope_unsplit_hydro_1d_scalar(qMinus.X, q.X, qPlus.X);

        return dq;
    }

    Real m_gamma;
    Kokkos::Array<Array3d, 2*three_d> m_qr;
    ConstArray3d m_q;
    UniformGrid m_grid;
    RealVector m_dtdl;
    Real m_slope_type;
};

}
