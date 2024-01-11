#pragma once

#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include "inih/INIReader.hpp"

namespace hydro { namespace problems
{

struct GreshoParams
{
    GreshoParams(const INIReader& reader)
    {
        using namespace constants;

        OC[IX] = reader.GetReal("problem", "xcenter", OC[IX]);
        OC[IY] = reader.GetReal("problem", "ycenter", OC[IY]);
        rho    = reader.GetReal("problem", "density", rho);
        mach   = reader.GetReal("problem", "mach", mach);
        p0     = rho / (mach * mach);
    }

    RealVector3d OC = {{constants::half, constants::half}};
    Real mach = constants::tenth;
    Real rho = constants::one;
    Real p0 = constants::zero;
};

}}
