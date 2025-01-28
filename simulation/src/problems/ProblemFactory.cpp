#include "ProblemFactory.hpp"

#include "HydroTypes.hpp"
#include "HydroParams.hpp"
#include "HydroProblem.hpp"
#include "RiemannProblem.hpp"
#include "Riemann2dProblem.hpp"
#include "GreshoProblem.hpp"
#include "OTangProblem.hpp"
#include "RayleighTaylorProblem.hpp"
#include "ImplodeProblem.hpp"
#include "BlastProblem.hpp"
#include "Blast_low_PB_Problem.hpp"
#include "RotorProblem.hpp"
#include "SkewedShockProblem.hpp"
#include "KevinHelmoltzMHDProblem.hpp"
#include "ConvectionProblem.hpp"




#include <exception>
#include <memory>

namespace hydro { namespace problems
{


std::shared_ptr<Problem> ProblemFactory::New(const std::string& problem_name,
                                                               const std::shared_ptr<Params>& params)
{
    std::shared_ptr<Problem> problem = nullptr;
     if (problem_name == "riemann")
    {
        problem = std::make_shared<RiemannProblem>(params);
    }
    else if (problem_name == "implode")
    {
        problem = std::make_shared<ImplodeProblem>(params);
    }
    else if (problem_name == "blast")
    {
        problem = std::make_shared<BlastProblem>(params);
    }
    else if (problem_name == "blast_low_pb")
    {
        problem = std::make_shared<Blast_low_PB_Problem>(params);
    }
    else if (problem_name == "rotor")
    {
        problem = std::make_shared<RotorProblem>(params);
    }
    else if (problem_name == "gresho")
    {
        problem = std::make_shared<GreshoProblem>(params);
    }
    else if (problem_name == "OTang")
    {
        problem = std::make_shared<OTangProblem>(params);
    }
    else if (problem_name == "rayleigh_taylor")
    {
        problem = std::make_shared<RayleighTaylorProblem>(params);
    }
    else if (problem_name == "riemann2d")
    {
        problem = std::make_shared<Riemann2dProblem>(params);
    }
    else if (problem_name == "skewed_shock")
    {
        problem = std::make_shared<SkewedShockProblem>(params);
    }
    else if (problem_name == "kevin_helmoltz_mhd")
    {
        problem = std::make_shared<KevinHelmoltzMHDProblem>(params);
    }
    else if (problem_name == "convection")
    {
        problem = std::make_shared<ConvectionProblem>(params);
    }
    else

    {
        throw std::runtime_error("Unknown problem");
    }
    return problem;
}

}}
