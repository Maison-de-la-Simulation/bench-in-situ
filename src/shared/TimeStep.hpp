#pragma once

#include "HydroTypes.hpp"

#include <cmath>
#include <string>

namespace hydro
{

class TimeStep
{
public:
    TimeStep() noexcept = default;
    explicit TimeStep(Real value) noexcept
        : hydro ( value )
    {
    }
    Real hydro;
    Real min() const
    {
        Real dt {hydro};
        return dt;
    }
    std::string cfl() const
    {
        std::string name {"hydro"};
        return name;
    }
};

}  // hydro
