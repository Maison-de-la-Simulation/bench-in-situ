#pragma once

#include <string>

enum class writer_t : short
{
    unknown,
    vtk,
    pdi,
    hdf5
};

inline
writer_t s2writer(const std::string& name)
{
    writer_t type;
    if (name == "vtk")
    {
        type = writer_t::vtk;
    }
    else if (name == "pdi")
    {
        type = writer_t::pdi;
    }
    else if (name == "hdf5")
    {
        type = writer_t::hdf5;
    }
    else
    {
        type = writer_t::unknown;
    }
    return type;
}
