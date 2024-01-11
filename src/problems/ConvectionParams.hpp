#pragma once

#include "HydroConstants.hpp"
#include "HydroParams.hpp"
#include "inih/INIReader.hpp"

namespace hydro { namespace problems
{

struct ConvectionParams
{
    ConvectionParams(const INIReader& reader)
    {
        density_bottom = reader.GetReal("problem", "density_bottom", density_bottom);
        T_bottom       = reader.GetReal("problem", "T_bottom",             T_bottom);
        u0_bottom      = reader.GetReal("problem", "u0_bottom",           u0_bottom);
        Bx_bottom      = reader.GetReal("problem", "Bx_bottom",            Bx_bottom);
        By_bottom      = reader.GetReal("problem", "By_bottom",            By_bottom);
        Bz_bottom      = reader.GetReal("problem", "Bz_bottom",            Bz_bottom);
        X_bottom       = reader.GetReal("problem", "X_bottom",              X_bottom);
        grad_T         = reader.GetReal("problem", "grad_T",         grad_T);
        grad_u0        = reader.GetReal("problem", "grad_u0",        grad_u0);
        grad_X         = reader.GetReal("problem", "grad_X",         grad_X);
        amplitude_seed = reader.GetReal("problem", "amplitude_seed", amplitude_seed);
        phi_target = reader.GetReal("problem", "phi_target", phi_target);
        HT       = reader.GetReal("problem", "HT",              HT);
        QA       = reader.GetReal("problem", "QA",              QA);
        RX       = reader.GetReal("problem", "RX",              RX);
        g_z        = reader.GetReal("hydro", "g_z",              g_z);

        kx         = reader.GetReal("problem", "kx",         kx);
        ky         = reader.GetReal("problem", "ky",         ky);
        kz         = reader.GetReal("problem", "kz",         kz);


      //  s_ux = reader.GetInteger("problem", "s_ux", s_ux);
      //  s_Bx = reader.GetInteger("problem", "s_Bx", s_Bx);
      //  s_Bz = reader.GetInteger("problem", "s_Bz", s_Bz);

      //  extrapolation_Bx_bord = reader.GetBoolean("problem", "extrapolation_Bx_bord", extrapolation_Bx_bord);
      //  extrapolation_ux_bord = reader.GetBoolean("problem", "extrapolation_ux_bord", extrapolation_ux_bord);

      //  perio_Bx_bord = reader.GetBoolean("problem", "perio_Bx_bord", perio_Bx_bord);
    }

    Real density_bottom = 0.0;
    Real T_bottom       = 0.0;
    Real u0_bottom      = 0.0;
    Real Bx_bottom      = 0.0;
    Real By_bottom      = 0.0;
    Real Bz_bottom      = 0.0;
    Real X_bottom       = 0.0;

    Real grad_T         = 0.0;
    Real grad_u0        = 0.0;
    Real grad_X         = 0.0;
    Real amplitude_seed = 0.0;
    Real phi_target = 0.0;

    Real HT=0.0;
    Real QA=0.0;
    Real RX=0.0;
    Real g_z =0.0;

    Real kx = 1;
    Real ky = 0;
    Real kz = 1;


  //  Int s_ux = 1.0;
  //  Int s_Bx = 1.0;
  //  Int s_Bz = 1.0;

  //  bool extrapolation_ux_bord = false;
  //  bool extrapolation_Bx_bord = false;
  //  bool perio_Bx_bord = false;


};

}}
