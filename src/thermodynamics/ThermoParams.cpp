#include "ThermoParams.hpp"
#include <string>
#include "inih/INIReader.hpp"
#include "HydroConstants.hpp"
#include "HydroTypes.hpp"

namespace hydro { namespace thermodynamics
{

ThermoParams::ThermoParams(const INIReader& reader)
{
    type = reader.Get(section, "type", type);
    gamma = reader.GetReal(section, "gamma", gamma);
    mmw = reader.GetReal(section, "mmw", mmw);
    mmw1 = reader.GetReal(section, "mmw1", mmw1);
    mmw2 = reader.GetReal(section, "mmw2", mmw2);
    kB = reader.GetReal(section, "kB", kB);
    e0 = reader.GetReal(section, "e0", e0);
    p0 = reader.GetReal(section, "p0", p0);
}

}  // namespace thermodynamics

}  // hydro
