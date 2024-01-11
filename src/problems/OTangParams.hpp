#pragma once

#include "HydroParams.hpp"
#include "HydroTypes.hpp"

#include "inih/INIReader.hpp"

namespace hydro { namespace problems
{

struct OTangParams
{
    OTangParams(const INIReader& reader)
    {
        using namespace constants;

        OC[IX] = reader.GetReal("problem", "xcenter", OC[IX]);
        OC[IY] = reader.GetReal("problem", "ycenter", OC[IY]);
    }

    RealVector3d OC = {{constants::half, constants::half}};
};

}}
