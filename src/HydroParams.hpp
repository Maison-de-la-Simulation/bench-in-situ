#pragma once

#include "HydroConstants.hpp"
#include "HydroTypes.hpp"
#include "Utils.hpp"
#include "inih/INIReader.hpp"
#include "thermodynamics/ThermoParams.hpp"

#include <iostream>
#include <string>

namespace hydro
{

struct RunParams
{
    RunParams() = default;

    RunParams(const INIReader& reader);

    Int nStepmax = 0;
    Int info = 100;
    Real tEnd = constants::infinity;
    Real cfl = constants::half;
    std::string solver = "unknown";
    std::string riemann = "unknown";
    bool useRangePolicy = true;
    bool restart = false;
    Real slope_type = 0.0;
    bool muscl_enabled = false;
    bool all_regime_correction = false;
    int all_regime_correction_type = 0;
    std::string restart_filename = "";
    std::string section = "run";
    Real beta_threshold = -1.0;
    bool powell_st_when_low_plasma_beta = false;
    bool entropic_correction=false;
    bool random_perturbation=false;
    bool three_d_perturbation=true;

};

struct HydroParams
{
    using RealVector = RealVector3d;

    HydroParams() = default;

    HydroParams(const INIReader& reader);

    bool hydro_enabled = true;
    bool convection_source_term_enabled = false;
    RealVector g = {utils::make_array<Real, three_d>(constants::zero)};
    Real K = constants::one+constants::tenth;


    Real resistivity_coefficient=0.0;
    bool magnetic_resistivity_enabled = false;

    bool splitted_diffusion_step = false;
    bool conservative_diffusion_update=false;

    std::string section = "hydro";
};

struct MeshParams
{
    using IntVector  = IntVectorNd<three_d>;
    using RealVector = RealVector3d;

    MeshParams() = default;

    MeshParams(const INIReader& reader);

    IntVector nbCells = {utils::make_array<Int, three_d>(2)};
    IntVector dom = {utils::make_array<Int, three_d>(1)};
    RealVector low = {utils::make_array<Real, three_d>(constants::zero)};
    RealVector up = {utils::make_array<Real, three_d>(constants::one)};
    Kokkos::Array<int, 2 * three_d> boundaryTypes = {utils::make_array<int, 2 * three_d>(boundary_t::PERIODIC)};
    std::string section = "mesh";
};

struct OutputParams
{
    OutputParams() = default;

    OutputParams(const INIReader& reader);

    std::string type = "vtk";
    std::string directory = ".";
    std::string format = "appended";
    std::string prefix = "output";
    Real dt_io = 0;
    Real dt_mean = 0;
    Real dt_profile = 0;
    Real dt_slice = 0;
    Int nOutput = 0;
    Int n_mean = 0;
    Int n_profile = 0;
    Int n_slice = 0;
    std::string section = "output";
};

struct Params
{
    using RealVector = RealVector3d;

    Params() = default;

    Params(const std::string& filename);

    void print(std::ostream& os) const;

    INIReader reader = INIReader("");
    RunParams run = RunParams();
    HydroParams hydro = HydroParams();
    MeshParams mesh = MeshParams();
    OutputParams output = OutputParams();
    thermodynamics::ThermoParams thermo = thermodynamics::ThermoParams();
    std::string problem = "unknown";
};

}   // namespace hydro
