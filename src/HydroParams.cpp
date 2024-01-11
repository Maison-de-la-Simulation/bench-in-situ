#include "HydroParams.hpp"

#include "HydroTypes.hpp"
#include "inih/INIReader.hpp"
#include "Print.hpp"
#include "utils/Stringify.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace hydro
{

const std::vector<std::string> X{"x", "y", "z"};

RunParams::RunParams(const INIReader& reader)
{
    nStepmax = reader.GetInteger(section, "nStepmax", nStepmax);
    info     = reader.GetInteger(section, "info", info);
    tEnd     = reader.GetReal(section, "tEnd", tEnd);
    slope_type = reader.GetReal(section, "slope_type", slope_type);
    muscl_enabled = reader.GetBoolean(section, "muscl_enabled", muscl_enabled);
    all_regime_correction = reader.GetBoolean(section, "all_regime_correction", all_regime_correction);
    all_regime_correction_type = reader.GetInteger(section, "all_regime_correction_type", all_regime_correction_type);
    cfl      = reader.GetReal(section, "cfl", cfl);
    solver   = reader.Get(section, "solver", solver);
    riemann  = reader.Get(section, "riemann", riemann);
    useRangePolicy = reader.GetBoolean(section, "useRangePolicy", useRangePolicy);

    restart = reader.GetBoolean(section, "restart", restart);
    restart_filename = reader.Get(section, "restart_filename", restart_filename);

    beta_threshold =  reader.GetReal(section, "beta_threshold", beta_threshold);
    powell_st_when_low_plasma_beta = reader.GetBoolean(section, "powell_st_when_low_plasma_beta", powell_st_when_low_plasma_beta);
    entropic_correction =  reader.GetBoolean(section, "entropic_correction", entropic_correction);
    three_d_perturbation =  reader.GetBoolean(section, "three_d_perturbation", three_d_perturbation);
    random_perturbation =  reader.GetBoolean(section, "random_perturbation", random_perturbation);

}

HydroParams::HydroParams(const INIReader& reader)
{
    hydro_enabled = reader.GetBoolean(section, "hydro_enabled", hydro_enabled);
    convection_source_term_enabled = reader.GetBoolean(section, "convection_source_term_enabled", convection_source_term_enabled);

    for (int idim = 0; idim < three_d; ++idim)
    {
        g[idim] = reader.GetReal(section, "g_" + X[idim], g[idim]);
    }
    K      = reader.GetReal(section, "K", K);

    magnetic_resistivity_enabled = reader.GetBoolean(section, "magnetic_resistivity_enabled", magnetic_resistivity_enabled);
    resistivity_coefficient  = reader.GetReal(section, "resistivity_coefficient", resistivity_coefficient);

    splitted_diffusion_step = reader.GetBoolean(section, "splitted_diffusion_step", splitted_diffusion_step);
    conservative_diffusion_update = reader.GetBoolean(section, "conservative_diffusion_update", conservative_diffusion_update);


}

MeshParams::MeshParams(const INIReader& reader)
{
    for (int idim = 0; idim < three_d; ++idim)
    {
        nbCells[idim] = reader.GetInteger(section, 'n' + X[idim], nbCells[idim]);
        dom[idim] = reader.GetInteger(section, 'm' + X[idim], dom[idim]);
        low[idim] = reader.GetReal(section, X[idim] + "min", low[idim]);
        up[idim]  = reader.GetReal(section, X[idim] + "max", up[idim]);
        boundaryTypes[IL+2*idim] = reader.GetInteger(
            section, "boundary_type_" + X[idim] + "min", boundaryTypes[IL+2*idim]);
        boundaryTypes[IR+2*idim] = reader.GetInteger(
            section, "boundary_type_" + X[idim] + "max", boundaryTypes[IR+2*idim]);
    }
}

OutputParams::OutputParams(const INIReader& reader)
{
    type      = reader.Get(section, "type", type);
    directory = reader.Get(section, "directory", directory);
    format    = reader.Get(section, "format", format);
    prefix    = reader.Get(section, "prefix", prefix);
    dt_io     = reader.GetReal(section, "dt_io", dt_io);
    dt_mean     = reader.GetReal(section, "dt_mean", dt_mean);
    dt_profile     = reader.GetReal(section, "dt_profile", dt_profile);
    dt_slice     = reader.GetReal(section, "dt_slice", dt_slice);
    nOutput   = reader.GetInteger(section, "nOutput", nOutput);
    n_mean   = reader.GetInteger(section, "n_mean", n_mean);
    n_profile   = reader.GetInteger(section, "n_profile", n_profile);
    n_slice   = reader.GetInteger(section, "n_slice", n_slice);
}

Params::Params(const std::string& filename)
    : reader{INIReader{filename}}
    , run{reader}
    , hydro{reader}
    , mesh{reader}
    , output{reader}
    , thermo{reader}
{
    problem = reader.Get("problem", "name", problem);
}

void Params::print(std::ostream& os) const
{
    const std::vector<std::string> v2{"i", "j", "k"};
    const int len1{40};
    const int len2{40};
    const int len3{10};

    os << std::scientific;
    os << std::setprecision(std::numeric_limits<Real>::digits10);

    os << std::string(len1+len2, '#') << std::endl;
    os << std::string(len1-len3, '#');
    os << std::left << std::setw(len2+len3) << std::setfill('#') << " Run parameters " << '\n';
    os << utils::stringify("solver", run.solver) << "\n";
    os << utils::stringify("riemann", run.riemann) << "\n";
    os << utils::stringify("cfl", run.cfl) << "\n";
    os << utils::stringify("tEnd", run.tEnd) << "\n";
    os << utils::stringify("slope_type", run.slope_type) << "\n";
    os << utils::stringify("muscl_enabled", run.muscl_enabled) << "\n";
    os << utils::stringify("all_regime_correction", run.all_regime_correction) << "\n";
    os << utils::stringify("all_regime_correction_type", run.all_regime_correction_type) << "\n";
    os << utils::stringify("nStepmax", run.nStepmax) << "\n";
    os << utils::stringify("info", run.info) << "\n";
  //  os << utils::stringify("useRangePolicy", run.useRangePolicy) << "\n";
    os << utils::stringify("restart", run.restart) << "\n";
    os << utils::stringify("restart_filename", run.restart_filename) << "\n";
    os << utils::stringify("powell_st_when_low_plasma_beta :", run.powell_st_when_low_plasma_beta) << "\n";
    os << utils::stringify("beta_threshold :", run.beta_threshold) << "\n";
    os << utils::stringify("random_perturbation :", run.random_perturbation) << "\n";
    os << utils::stringify("three_d_perturbation :", run.three_d_perturbation) << "\n";

    os << std::string(len1-len3, '#');
    os << std::left << std::setw(len2+len3) << std::setfill('#') << " Hydro parameters " << '\n';
  //  os << utils::stringify("magnetic_resistivity_enabled", hydro.magnetic_resistivity_enabled) << "\n";
  //  os << utils::stringify("resistivity_coefficient", hydro.resistivity_coefficient) << "\n";
  //  os << utils::stringify("splitted_diffusion_step", hydro.splitted_diffusion_step) << "\n";
  //  os << utils::stringify("conservative_diffusion_update", hydro.conservative_diffusion_update) << "\n";
    os << utils::stringify("convection_source_term_enabled", hydro.convection_source_term_enabled) << "\n";
    os << utils::stringify("H_source_term", reader.GetBoolean("hydro", "H_source_term", false)) << "\n";
    os << utils::stringify("R_source_term", reader.GetBoolean("hydro", "R_source_term", false)) << "\n";
    os << utils::stringify("Q_source_term", reader.GetBoolean("hydro", "Q_source_term", false)) << "\n";
    os << utils::stringify("HT", reader.GetReal("problem", "HT", 0.0)) << "\n";
    os << utils::stringify("RX", reader.GetReal("problem", "RX", 0.0)) << "\n";
    os << utils::stringify("QA", reader.GetReal("problem", "QA", 0.0)) << "\n";

    //os << utils::stringify("K", hydro.K) << "\n";
    for (int idim = 0; idim < three_d; ++idim)
    {
        os << utils::stringify('g' + X[idim], hydro.g[idim]) << "\n";
    }

    os << std::string(len1-len3, '#');
    os << std::left << std::setw(len2+len3) << std::setfill('#') << " Thermo parameters " << '\n';
    os << utils::stringify("Eos type", thermo.type) << "\n";
    os << utils::stringify("mmw", thermo.mmw) << "\n";
    os << utils::stringify("gamma", thermo.gamma) << "\n";
  //  os << utils::stringify("p0", thermo.p0) << "\n";

    os << std::string(len1-len3, '#');
    os << std::left << std::setw(len2+len3) << std::setfill('#') << " Mesh parameters " << '\n';
    for (int idim = 0; idim < three_d; ++idim)
    {
        os << utils::stringify('n' + X[idim], mesh.nbCells[idim]) << "\n";
    }
    for (int idim = 0; idim < three_d; ++idim)
    {
        os << utils::stringify('m' + X[idim], mesh.dom[idim]) << "\n";
    }
    for (int idim = 0; idim < three_d; ++idim)
    {
        os << utils::stringify(X[idim] + "min", mesh.low[idim]) << "\n";
        os << utils::stringify(X[idim] + "max", mesh.up[idim]) << "\n";
    }
    for (int idim = 0; idim < three_d; ++idim)
    {
        os << utils::stringify("boundary type " + X[idim] + "min", mesh.boundaryTypes[0+2*idim]) << "\n";
        os << utils::stringify("boundary type " + X[idim] + "max", mesh.boundaryTypes[1+2*idim]) << "\n";
    }

    os << std::string(len1-len3, '#');
    os << std::left << std::setw(len2+len3) << std::setfill('#') << " Output parameters " << '\n';
    os << utils::stringify("type", output.type) << "\n";
  //  os << utils::stringify("directory", output.directory) << "\n";
    os << utils::stringify("format", output.format) << "\n";
    os << utils::stringify("prefix", output.prefix) << "\n";
    os << utils::stringify("nOutput", output.nOutput) << "\n";
    os << utils::stringify("n_mean", output.n_mean) << "\n";
    os << utils::stringify("n_profile", output.n_profile) << "\n";
    os << utils::stringify("n_slice", output.n_slice) << "\n";
    os << utils::stringify("dt_io", output.dt_io) << "\n";
    os << utils::stringify("dt_mean :", output.dt_mean) << "\n";
    os << utils::stringify("dt_profile :", output.dt_profile) << "\n";
    os << utils::stringify("dt_slice :", output.dt_slice) << "\n";
    os << std::string(len1 + len2, '#') << std::endl;
}

}  // hydro
