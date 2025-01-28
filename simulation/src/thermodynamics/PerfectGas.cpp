#include "PerfectGas.hpp"
#include "HydroConstants.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "ThermoParams.hpp"

namespace hydro { namespace thermodynamics
{

PerfectGas::PerfectGas(Real gamma, Real mmw, Real mmw1, Real mmw2, Real e0, Real kB) noexcept
    : m_gamma {gamma}
    , m_gamma_m1 {gamma-constants::one}
    , m_inv_gamma_m1 {constants::one/(gamma-constants::one)}
    , m_e0 {e0}
    , m_mmw{mmw}
    , m_mmw1{mmw1}
    , m_mmw2{mmw2}
    , m_kB {kB}

{
}

PerfectGas::PerfectGas(const ThermoParams& thermoParams)
    : PerfectGas(thermoParams.gamma, thermoParams.mmw, thermoParams.mmw1, thermoParams.mmw2, thermoParams.e0, thermoParams.kB)
{
}

}  // namespace thermodynamics

}  // hydro
