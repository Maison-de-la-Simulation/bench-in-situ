#pragma once

#include <Kokkos_Core.hpp>
#include "HydroConstants.hpp"
#include "HydroTypes.hpp"
#include "ThermoParams.hpp"

#include <cmath>

namespace hydro { namespace thermodynamics
{

//!
//! @class PerfectGas
//! @brief This class represents a perfect i.e. following:
//! @f$p = (\gamma - 1) \rho (e - e_0)@f$. Domain of validity
//! is @f$\mathcal{D} = \{\rho > 0, e > e_0\}@f$.
//!
class PerfectGas
{
public:
    //!
    //! @fn PerfectGas (Real gamma, Real, mmw, Real p0, Real e0) noexcept
    //! @brief Constructs a perfect gas.
    //! @param[in] gamma Adiabatic index of the perfect gas.
    //! @param[in] mmw Mean molecular weight of the perfect gas.
    //! @param[in] e0 Energy shift of the perfect gas.
    //!
    PerfectGas(Real gamma, Real mmw, Real mmw1, Real mmw2, Real e0, Real kB) noexcept;

    //!
    //! @fn PerfectGas(const ThermoParams& thermoParams)
    //!
    PerfectGas(const ThermoParams& thermoParams);

    //!
    //! @fn ~PerfectGas()
    //!
    ~PerfectGas() = default;

    //!
    //! @fn PerfectGas(const PerfectGas& x)
    //!
    PerfectGas(const PerfectGas& x) = default;

    //!
    //! @fn PerfectGas(PerfectGas&& x)
    //!
    PerfectGas(PerfectGas&& x) = default;

    //!
    //! @fn PerfectGas& operator=(const PerfectGas& x)
    //!
    PerfectGas& operator=(const PerfectGas& x) = default;

    //!
    //! @fn PerfectGas& operator=(PerfectGas&& x)
    //!
    PerfectGas& operator=(PerfectGas&& x) = default;

    //!
    //! @fn Real computeAdiabaticIndex() const noexcept
    //! @return Adiabatic index.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeAdiabaticIndex() const noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    Real computeInvGammam1() const noexcept;

    //!
    //! @fn Real computeInternalEnergy(Real d, Real p) const noexcept
    //! @param[in] d Density.
    //! @param[in] p Pressure.
    //! @return Internal energy.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeInternalEnergy(Real d, Real p) const noexcept;

    //!
    //! @fn Real computeMeanMolecularWeight() const noexcept
    //! @return Mean molecular weight.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeMeanMolecularWeight(Real X) const noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    Real computeMeanMolecularWeight() const noexcept;

    //!
    //! @fn Real computePressure(Real d, Real e) const noexcept
    //! @param[in] d Density.
    //! @param[in] e Internal energy.
    //! @return Pressure.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computePressure(Real d, Real e) const noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    Real computePressure(Real d, Real T, Real X) const noexcept;

    //!
    //! @fn Real computeSpeedOfSound(Real d, Real p) const noexcept
    //! @brief @f$ c = \sqrt{\gamma \frac{p}{\rho}}@f$.
    //! @param[in] d Density.
    //! @param[in] p Pressure.
    //! @return Speed of sound.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeSpeedOfSound(Real d, Real p) const noexcept;

    //!
    //! @fn Real computeTemperature(Real d, Real p) const noexcept
    //! @brief @f$ p = \rho R_* T@f$, for @f$R_* = \mu R@f$.
    //! @param[in] d Density.
    //! @param[in] p Pressure.
    //! @return Temperature.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeTemperature(Real d, Real p, Real X) const noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    Real computeLogPotentialTemperature(Real T, Real p) const noexcept;

private:
    Real m_gamma;
    Real m_gamma_m1;
    Real m_inv_gamma_m1;
    Real m_e0;
    Real m_mmw;
    Real m_mmw1;
    Real m_mmw2;
    Real m_kB;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation details
////////////////////////////////////////////////////////////////////////////////

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeAdiabaticIndex() const noexcept
{
    return m_gamma;
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeInvGammam1() const noexcept
{
    return m_inv_gamma_m1;
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeInternalEnergy(Real d, Real p) const noexcept
{
    return m_inv_gamma_m1 * p + d * m_e0;
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeMeanMolecularWeight(Real X) const noexcept
{
    return 1.0/(X/m_mmw1 + (1.0-X)/m_mmw2);
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeMeanMolecularWeight() const noexcept
{
    return m_mmw;
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computePressure(Real d, Real T, Real X) const noexcept
{
    return d*m_kB*T/computeMeanMolecularWeight(X);
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computePressure(Real d, Real e) const noexcept
{
    return m_gamma_m1 * (e - d * m_e0);
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeSpeedOfSound(Real d, Real p) const noexcept
{
    return std::sqrt(m_gamma * p / d);
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeTemperature(Real d, Real p, Real X) const noexcept
{
    return computeMeanMolecularWeight(X)*p/(m_kB*d);
}

KOKKOS_FORCEINLINE_FUNCTION
Real PerfectGas::computeLogPotentialTemperature(Real T, Real p) const noexcept
{
    return std::log(T) -((m_gamma-1.0)/m_gamma)*std::log(p);
}


}   // namespace thermodynamics

}   // namespace hydro
