#pragma once

#include <Kokkos_Core.hpp>
#include "HydroConstants.hpp"
#include "HydroTypes.hpp"
#include "ThermoParams.hpp"

#include <cmath>

namespace hydro { namespace thermodynamics
{

//!
//! @class StiffenedGas
//! @brief This class represents a stiffened gas i.e. following:
//! @f$p = (\gamma - 1) \rho (e - e_0) - \gamma \Pi@f$.
//!
class StiffenedGas
{
public:
    //!
    //! @fn StiffenedGas (Real gamma, Real, mmw, Real p0, Real e0) noexcept
    //! @brief  Constructs a stiffened gas.
    //! @param[in] gamma Adiabatic index of the stiffened gas.
    //! @param[in] mmw Mean molecular weight of the stiffened gas.
    //! @param[in] p0 Pressure shift of the stiffened gas.
    //! @param[in] e0 Energy shift of the stiffened gas.
    //!
    StiffenedGas(Real gamma, Real mmw, Real p0, Real e0) noexcept;

    //!
    //! @fn StiffenedGas (const ThermoParams& thermoParams)
    //!
    StiffenedGas(const ThermoParams& thermoParams);

    //!
    //! @fn ~StiffenedGas ()
    //!
    ~StiffenedGas() = default;

    //!
    //! @fn StiffenedGas (const StiffenedGas& x)
    //!
    StiffenedGas(const StiffenedGas& x) = default;

    //!
    //! @fn StiffenedGas (StiffenedGas&& x)
    //!
    StiffenedGas(StiffenedGas&& x) = default;

    //!
    //! @fn StiffenedGas& operator=(const StiffenedGas& x)
    //!
    StiffenedGas& operator=(const StiffenedGas& x) = default;

    //!
    //! @fn StiffenedGas& operator=(StiffenedGas&& x)
    //!
    StiffenedGas& operator=(StiffenedGas&& x) = default;

    //!
    //! @fn Real computeAdiabaticIndex() const noexcept
    //! @return Adiabatic index.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeAdiabaticIndex() const noexcept;

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
    Real computeMeanMolecularWeight() const noexcept;

    //!
    //! @fn Real computePressure(Real d, Real e) const noexcept
    //! @param[in] d Density.
    //! @param[in] e Internal energy.
    //! @return Pressure.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computePressure(Real d, Real e) const noexcept;

    //!
    //! @fn Real computeSpeedOfSound(Real d, Real p) const noexcept
    //! @param[in] d Density.
    //! @param[in] p Pressure.
    //! @return Speed of sound.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeSpeedOfSound(Real d, Real p) const noexcept;

    //!
    //! @fn Real computeTemperature(Real d, Real p) const noexcept
    //! @param[in] d Density.
    //! @param[in] p Pressure.
    //! @return Temperature.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeTemperature(Real d, Real p) const noexcept;

    //!
    //! @fn Real computeP0() const noexcept
    //! @return Constant p0.
    //!
    KOKKOS_FORCEINLINE_FUNCTION
    Real computeP0() const noexcept;

private:
    Real m_gamma;
    Real m_gamma_m1;
    Real m_inv_gamma_m1;
    Real m_mmw;
    Real m_Rstar;
    Real m_p0;
    Real m_e0;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation details
////////////////////////////////////////////////////////////////////////////////

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computeAdiabaticIndex() const noexcept
{
    return m_gamma;
}

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computeInternalEnergy(Real d, Real p) const noexcept
{
    return m_inv_gamma_m1 * (p + m_gamma * m_p0) + d * m_e0;
}

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computeMeanMolecularWeight() const noexcept
{
    return m_mmw;
}

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computePressure(Real d, Real e) const noexcept
{
    return m_gamma_m1 * (e - d * m_e0) - m_gamma * m_p0;
}

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computeSpeedOfSound(Real d, Real p) const noexcept
{
    return std::sqrt(m_gamma * (p + m_p0) / d);
}

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computeP0() const noexcept
{
    return m_p0;
}

KOKKOS_FORCEINLINE_FUNCTION
Real StiffenedGas::computeTemperature(Real d, Real p) const noexcept
{
    return (p + m_p0) / (d * m_Rstar);
}

}  // thermodynamics

}  // hydro
