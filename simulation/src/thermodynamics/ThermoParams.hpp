#pragma once

#include <string>
#include "inih/INIReader.hpp"
#include "HydroConstants.hpp"
#include "HydroTypes.hpp"

namespace hydro { namespace thermodynamics
{

class ThermoParams
{
public:
    ThermoParams() = default;
    ThermoParams(const INIReader& reader);
    ThermoParams(const ThermoParams& x) = default;
    ThermoParams(ThermoParams&& x) = default;
    ~ThermoParams() = default;
    ThermoParams& operator=(const ThermoParams& x) = default;
    ThermoParams& operator=(ThermoParams&& x) = default;

    std::string type = {"perfect gas"};
    Real gamma = {constants::five_third};
    Real mmw = {constants::one};
    Real mmw1 = {1.0};
    Real mmw2 = {1.0};
    Real kB = {1.0};
    Real p0 = {constants::zero};
    Real e0 = {constants::zero};
    std::string section = {"thermo"};
};

}  // namespace thermodynamics

}  // hydro
