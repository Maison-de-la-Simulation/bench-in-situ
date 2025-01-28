#pragma once

#include "HydroTypes.hpp"

#include <cmath>
#include <string>

namespace hydro
{

class global_mean
{
public:
    global_mean() noexcept = default;

    explicit global_mean(Real value) noexcept
        : m_emag(value)
        , m_ekin(value)
    {
    }

    Real get_emag() const
    {
        return m_emag;
    }

    Real get_ekin() const
    {
        return m_ekin;
    }

    Real m_emag;
    Real m_ekin;
};

}  // hydro
