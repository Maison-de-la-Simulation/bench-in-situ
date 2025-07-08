#pragma once

#include <string>

enum class writer_t : short
{
    unknown,
    vtk,
    pdi,
    gpu_pdi
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
    else if (name == "gpu_pdi")
    {
        type = writer_t::gpu_pdi;
    }
    else
    {
        type = writer_t::unknown;
    }
    return type;
}
