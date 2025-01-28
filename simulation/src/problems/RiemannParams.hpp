#pragma once

#include "HydroTypes.hpp"
#include "HydroConstants.hpp"
#include "HydroParams.hpp"
#include "inih/INIReader.hpp"

namespace hydro { namespace problems
{

struct RiemannParams
{
    RiemannParams(const INIReader& reader)
    {
        density_left   = reader.GetReal("problem", "density_left", density_left);
        density_right  = reader.GetReal("problem", "density_right", density_right);

        pressure_left  = reader.GetReal("problem", "pressure_left", pressure_left);
        pressure_right = reader.GetReal("problem", "pressure_right", pressure_right);

        x_velocity_left  = reader.GetReal("problem", "x_velocity_left", x_velocity_left);
        x_velocity_right = reader.GetReal("problem", "x_velocity_right", x_velocity_right);

        y_velocity_left  = reader.GetReal("problem", "y_velocity_left", y_velocity_left);
        y_velocity_right = reader.GetReal("problem", "y_velocity_right", y_velocity_right);

        z_velocity_left  = reader.GetReal("problem", "z_velocity_left", z_velocity_left);
        z_velocity_right = reader.GetReal("problem", "z_velocity_right", z_velocity_right);

        x_Mag_Field_left  = reader.GetReal("problem", "x_Mag_Field_left", x_Mag_Field_left);
        x_Mag_Field_right = reader.GetReal("problem", "x_Mag_Field_right", x_Mag_Field_right);

        y_Mag_Field_left  = reader.GetReal("problem", "y_Mag_Field_left", y_Mag_Field_left);
        y_Mag_Field_right = reader.GetReal("problem", "y_Mag_Field_right", y_Mag_Field_right);

        z_Mag_Field_left  = reader.GetReal("problem", "z_Mag_Field_left", z_Mag_Field_left);
        z_Mag_Field_right = reader.GetReal("problem", "z_Mag_Field_right", z_Mag_Field_right);

        direction        = reader.GetInteger("problem", "direction", IX);

        center         = reader.GetReal("problem", "center", center);
    }

    Real density_left = constants::zero;
    Real pressure_left = constants::zero;
    Real density_right = constants::zero;
    Real pressure_right = constants::zero;
    Real x_velocity_left = constants::zero;
    Real x_velocity_right = constants::zero;
    Real y_velocity_left = constants::zero;
    Real y_velocity_right = constants::zero;
    Real z_velocity_left = constants::zero;
    Real z_velocity_right = constants::zero;
    Real x_Mag_Field_left = constants::zero;
    Real x_Mag_Field_right = constants::zero;
    Real y_Mag_Field_left = constants::zero;
    Real y_Mag_Field_right = constants::zero;
    Real z_Mag_Field_left = constants::zero;
    Real z_Mag_Field_right = constants::zero;
    Real center = constants::half;
    Int  direction = IX;
};

}}
