#pragma once

#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "inih/INIReader.hpp"

namespace hydro { namespace problems
{

struct RayleighTaylorParams
{
    RayleighTaylorParams(const INIReader& reader)
    {
        density_up      = reader.GetReal("problem", "density_up", density_up);
        density_bottom  = reader.GetReal("problem", "density_bottom", density_bottom);
        pressure_bottom = reader.GetReal("problem", "pressure_bottom", pressure_bottom);
        perturbation    = reader.GetReal("problem", "perturbation", perturbation);
    }

    Real density_up = constants::two;
    Real density_bottom = constants::one;
    Real pressure_bottom = constants::five*constants::half;
    Real perturbation = constants::tenth*constants::tenth;
};

}}
