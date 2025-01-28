#pragma once

#include "HydroTypes.hpp"

#include <cmath>
#include <string>

namespace hydro
{

class slice_average
{
public:
    slice_average() = default;

    explicit slice_average(const Real value) noexcept
        : m_Bx2(value)
        , m_By2(value)
        , m_Bz2(value)
        , m_u(value)
        , m_v(value)
        , m_logmu(value)
        , m_logtheta(value)
    {
    }

    Real m_Bx2;
    Real m_By2;
    Real m_Bz2;
    Real m_u;
    Real m_v;
    Real m_logmu;
    Real m_logtheta;
};

}  // hydro
