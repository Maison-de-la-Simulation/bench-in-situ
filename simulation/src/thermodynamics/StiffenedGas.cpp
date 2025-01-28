#include "StiffenedGas.hpp"
#include "HydroConstants.hpp"
#include "HydroTypes.hpp"
#include "HydroUnits.hpp"
#include "ThermoParams.hpp"

namespace hydro { namespace thermodynamics
{

StiffenedGas::StiffenedGas(Real gamma, Real mmw, Real p0, Real e0) noexcept
    : m_gamma {gamma}
    , m_gamma_m1 {gamma-constants::one}
    , m_inv_gamma_m1 {constants::one/(gamma-constants::one)}
    , m_mmw {mmw}
    , m_Rstar {code_units::constants::Rstar_h / mmw}
    , m_p0 {p0}
    , m_e0 {e0}
{
}

StiffenedGas::StiffenedGas(const ThermoParams& thermoParams)
    : StiffenedGas {thermoParams.gamma, thermoParams.mmw,
                    thermoParams.p0, thermoParams.e0}
{
}

}  // namespace thermodynamics

}  // namespace hydro
