#include "Writer.hpp"

#include "DistributedMemorySession.hpp"
#include "HydroUniformGrid.hpp"
#include "HydroParams.hpp"
#include "HydroTypes.hpp"
#include "WriterBase.hpp"
#if defined(Euler_ENABLE_PDI)
#include "WriterPDI.hpp"
#endif
#include "WriterTypes.hpp"
#include "WriterVTK.hpp"

#include <chrono>
#include <list>
#include <memory>
#include <stdexcept>
#include <vector>

namespace hydro { namespace io
{

std::shared_ptr<WriterBase> WriterFactory::New(const UniformGrid& grid, const Params& params,
                                                         const std::string& type, const std::string& prefix,
                                                         const std::vector<std::pair<int, std::string>>& variables)
{
    std::shared_ptr<WriterBase> ptr ( nullptr );
    if (s2writer(type) == writer_t::vtk)
    {
        ptr = std::make_shared<WriterVTK>(grid, params, prefix, variables);
    }
    else if (s2writer(type) == writer_t::pdi)
    {
#if defined(Euler_ENABLE_PDI)
        ptr = std::make_shared<WriterPDI>(grid, params, prefix, variables);
#else
        throw std::runtime_error("ARK has not been compiled with PDI support\n");
#endif
    }
    else
    {
        throw std::runtime_error("Unknown writer\n");
    }
    return ptr;
}

// template <dim_t three_d>
// void Writer<three_d>::write(HostConstArrayDyn u, const UniformGrid& grid,
//                         Int iStep, Real time, Real gamma, Real mmw) const
// {
//     Print() << " * I/O\n";
//     for (const auto& io_writer : io_writers)
//     {
//         io_writer->write(u, grid, iStep, time, gamma, mmw);
//         Session::synchronize();
//         double bandwidth = 1.0e-6 * (static_cast<double>(Session::getNProc()) *
//                                      static_cast<double>(grid.nbCells()) *
//                                      static_cast<double>(u.extent(1)) * // nbvar
//                                      static_cast<double>(sizeof(Real)) /
//                                      timer.template count<std::chrono::duration<double>>());
//         Print() << " ** Duration : " << timer.template count<std::chrono::duration<double>>() << " s\n";
//         Print() << " ** Bandwidth: " << bandwidth << " MB/s\n";
//     }
// }

struct WriterFactory;

}}
