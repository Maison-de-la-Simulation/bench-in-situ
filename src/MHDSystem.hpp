#pragma once

#include "HydroConstants.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "linalg/Vector.hpp"
#include "thermodynamics/PerfectGas.hpp"
#include "Utils.hpp"

#include <Kokkos_Core.hpp>
#include <vector>

namespace hydro
{

struct MHDVarPrimMD;

struct MHDVarPrimMD
{
  enum : int {ID=0, IP=1, IUx=2, IUy=3, IUz=4, IBx=5, IBy=6, IBz=7, I_X=8};
};

struct MHDVarConsMD;

struct MHDVarConsMD
{
  enum : int {ID=0, IE=1, IMx=2, IMy=3, IMz=4, IBx=5, IBy=6, IBz=7, I_dX=8};
};

struct MHDSystem
{
  static constexpr int dim {3};
  static constexpr int nbvar {2+2*3+1};

    using EquationOfState = thermodynamics::PerfectGas;

    using VarPrim = MHDVarPrimMD;
    using VarCons = MHDVarConsMD;

    struct ConsState
    {
        Real d;
        Real e; // change to E ?
        Vector<three_d, Real> m;
        Vector<three_d, Real> B;
        Real dX;

    };

    struct PrimState
    {
        Real d;
        Real p;
        Vector<three_d, Real> v;
        Vector<three_d, Real> B;
        Real X;

    };

    MHDSystem() = delete;
    MHDSystem(const MHDSystem& x) = delete;
    MHDSystem(MHDSystem&& x) = delete;
    ~MHDSystem() = delete;
    MHDSystem& operator=(const MHDSystem& x) = delete;
    MHDSystem& operator=(MHDSystem&& x) = delete;

    static std::vector<std::string> cons_names()
    {
        std::vector<std::string> names ( nbvar );

        names.at(VarCons::ID) = "d";
        names.at(VarCons::IE) = "E";
        names.at(VarCons::IMx+0) = "mx";
        names.at(VarCons::IBx+0) = "Bx";
        names.at(VarCons::IMx+1) = "my";
        names.at(VarCons::IBx+1) = "By";
        names.at(VarCons::IMx+2) = "mz";
        names.at(VarCons::IBx+2) = "Bz";
        names.at(VarCons::I_dX ) =  "dX";


        return names;
    }

    static std::vector<std::string> prim_names()
    {
        std::vector<std::string> names ( nbvar );
        names.at(VarPrim::ID) = "d";
        names.at(VarPrim::IP) = "p";
        names.at(VarPrim::IUx+0) = "ux";
        names.at(VarPrim::IBx+0) = "Bx";
        names.at(VarPrim::IUx+1) = "uy";
        names.at(VarPrim::IBx+1) = "By";
        names.at(VarPrim::IUx+2) = "uz";
        names.at(VarPrim::IBx+2) = "Bz";
        names.at(VarPrim::I_X  ) = "X";



        return names;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeMeanMolecularWeight(const PrimState& q, const EquationOfState& eos) noexcept
    {
        return eos.computeMeanMolecularWeight(q.X);
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeLogPotentialTemperature(const PrimState& q, const EquationOfState& eos) noexcept
    {
        Real T=computeTemperature(q,eos);
        return eos.computeLogPotentialTemperature(T,q.p);
    }
    KOKKOS_INLINE_FUNCTION
    static Real computeTemperature(const PrimState& q, const EquationOfState& eos) noexcept
    {
        return eos.computeTemperature(q.d, q.p, q.X);
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeInternalEnergy(const PrimState& q, const EquationOfState& eos) noexcept
    {
        return eos.computeInternalEnergy(q.d, q.p);
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeSpeedOfSound(const PrimState& q, const EquationOfState& eos) noexcept
    {
        return eos.computeSpeedOfSound(q.d, q.p);
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeAlfvenSpeed(const PrimState& q, const int& dir) noexcept
    {
        const Real ca = q.B(dir)/std::sqrt(q.d);
        return ca;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeSlowMagnetoAcousticSpeed(const PrimState& q, const EquationOfState& eos, const int& dir) noexcept
    {
      const Real B_norm22 = dot(q.B, q.B);
      const Real c    = eos.computeSpeedOfSound(q.d, q.p);
      const Real c02  = c*c;
      const Real ca2  = B_norm22/q.d;
      const Real cap2 = q.B(dir)*q.B(dir)/q.d;

      const Real cs = std::sqrt(0.5*(c02+ca2)-0.5*std::sqrt((c02+ca2)*(c02+ca2)-4.*c02*cap2));

      return cs;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeFastMagnetoAcousticSpeed(const PrimState& q, const EquationOfState& eos, const int& dir) noexcept
    {
      const Real B_norm22 = dot(q.B, q.B);
      const Real c    = eos.computeSpeedOfSound(q.d, q.p);
      const Real c02  = c*c;
      const Real ca2  = B_norm22/q.d;
      const Real cap2 = q.B(dir)*q.B(dir)/q.d;

      const Real cf = std::sqrt(0.5*(c02+ca2)+0.5*std::sqrt((c02+ca2)*(c02+ca2)-4.*c02*cap2));

      return cf;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeKineticEnergy(const PrimState& q) noexcept
    {
        const Real v_norm22 = dot(q.v, q.v);
        const Real ekin = constants::half * q.d * v_norm22;
        return ekin;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeKineticEnergy(const ConsState& u) noexcept
    {
        const Real m_norm22 = dot(u.m, u.m);
        const Real ekin = constants::half * m_norm22 / u.d;
        return ekin;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeMagneticEnergy(const PrimState& q) noexcept  // New for MHD
    {
        const Real B_norm22 = dot(q.B, q.B);
        const Real emag = constants::half * B_norm22;
        return emag;
    }

    KOKKOS_INLINE_FUNCTION
    static Real computeMagneticEnergy(const ConsState& u) noexcept  // New for MHD
    {
        const Real B_norm22 = dot(u.B, u.B);
        const Real emag = constants::half * B_norm22;
        return emag;
    }

    KOKKOS_INLINE_FUNCTION
    static PrimState conservativeToPrimitive(const ConsState& u, const EquationOfState& eos) noexcept
    {
        const Real ekin = computeKineticEnergy (u);
        const Real emag = computeMagneticEnergy(u);

        PrimState q;
        q.d = u.d;
        q.p = eos.computePressure(u.d, u.e - ekin - emag);
        q.v = u.m;
        q.v *= constants::one / u.d;
        q.B = u.B;
        q.X = u.dX/u.d;

        return q;
    }

    KOKKOS_INLINE_FUNCTION
    static ConsState primitiveToConservative(const PrimState& q, const EquationOfState& eos) noexcept
    {
        const Real eint = eos.computeInternalEnergy(q.d, q.p);
        const Real ekin = computeKineticEnergy (q);
        const Real emag = computeMagneticEnergy(q);

        ConsState u;
        u.d = q.d;
        u.e = eint + ekin + emag;
        u.m  = q.v;
        u.m *= q.d;
        u.B = q.B;
        u.dX = q.d*q.X;

        return u;
    }

    template <int dir>
    KOKKOS_INLINE_FUNCTION
    static ConsState Flux(const PrimState& q, const EquationOfState& eos) noexcept
    {
        const Real eint = eos.computeInternalEnergy(q.d, q.p);
        const Real ekin = computeKineticEnergy(q);
        const Real emag = computeMagneticEnergy(q);

        const Real vn = q.v(dir);
        const Real Bn = q.B(dir);

        const Real TotalEnergy = eint+ekin+emag;

        const Real p_total = q.p + emag;

        Vector<three_d,Real> tmp, tmp1;

        ConsState flux;
        flux.d = q.d * vn;
        flux.e = (TotalEnergy + p_total) * vn - Bn*dot(q.B,q.v);

        tmp = q.v; tmp*=q.d*vn;
        flux.m = tmp;

        tmp = q.B;  tmp*= Bn;
        flux.m -= tmp;

        flux.m(dir) += p_total;

        tmp = q.B; tmp*=vn;
        tmp1 = q.v; tmp1*=Bn;

        flux.B = tmp-tmp1;
        flux.dX = q.X*q.d*vn;

        return flux;
    }

    template <int dir>
    KOKKOS_INLINE_FUNCTION
    static ConsState Flux(const ConsState& u, const EquationOfState& eos) noexcept
    {
        const Real ekin = computeKineticEnergy (u);
        const Real emag = computeMagneticEnergy(u);
        const Real eint = u.e - ekin - emag;
        const Real p = eos.computePressure(u.d, eint);
        const Real vn = u.m(dir) / u.d;
        const Real Bn = u.B(dir);
        const Real p_total = p + emag;

        ConsState flux;

        Vector<three_d, Real> tmp, tmp1;
        flux.d = u.d * vn;
        flux.e = (u.e + p_total) * vn - Bn*dot(u.B,u.m)/u.d;

        tmp = u.m; tmp*=vn;
        flux.m = tmp;

        tmp = u.B; tmp*=Bn;
        flux.m -= tmp;

        flux.m(dir) += p_total;

        tmp=u.B; tmp*=vn;
        tmp1=u.m; tmp1*=Bn/u.d;

        flux.B = tmp-tmp1;
        flux.dX = u.dX * vn;

        return flux;
    }
};

}
