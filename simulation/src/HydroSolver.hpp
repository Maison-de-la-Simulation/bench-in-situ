#pragma once

#include "HydroTypes.hpp"
#include "HydroConstants.hpp"

namespace hydro
{

class Solver
{
public:

    Solver() = default;
    Solver(const Solver& x) = default;
    Solver(Solver&& x) = default;
    virtual ~Solver() = default;
    Solver& operator=(const Solver& x) = default;
    Solver& operator=(Solver&& x) = default;
    virtual Real computeTimeStep() = 0;
    virtual void nextIteration(Real dt) = 0;
    virtual void prepareNextOutput(Real& dt) = 0;
    virtual void pdiExposeData() = 0;
    virtual bool finished() const = 0;
    virtual void printMonitoring(double t_tot) const = 0;
    virtual bool shouldPrintInformation() const = 0;
    virtual void printInformation(Real dt) const = 0;
    virtual double memoryUsage() const = 0;
    virtual Real time() const;
    virtual Int iteration() const;
    virtual void set_should_save() = 0;
    virtual void set_time_limit_reached() = 0;
    virtual void accumulate_compute_duration(const std::chrono::steady_clock::duration& duration) = 0;

    Real m_t = constants::zero;
    Int m_iteration = 0;
};

inline
Real Solver::time() const
{
    return m_t;
}

inline
Int Solver::iteration() const
{
    return m_iteration;
}

}
