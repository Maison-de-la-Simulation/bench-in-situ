#pragma once

#include "HydroConstants.hpp"
#include "HydroTypes.hpp"

// more constants at https://physics.nist.gov/cuu/index.html

namespace hydro { namespace code_units
{

#if defined(FLOAT64_EXTENDED)
class SI_Constants
{
public:
    // Avogadro number in SI (mol^-1)
    static constexpr Real Na      {6.02214085700000000000000000000000000000000000000000000000000E+23L};
    // Boltzmann constant in SI (J.K^-1)
    static constexpr Real k_b     {1.38064852000000000000000000000000000000000000000000000000000E-23L};
    // Universal gas constant in SI (J.mol^-1.K^-1)
    static constexpr Real R       {8.31445980000000000000000000000000000000000000000000000000000E+00L};
    // Mean molecular weight of hydrogen (no unit)
    static constexpr Real mmw_h   {constants::one};
    // Mass of hydrogen in SI (kg)
    static constexpr Real m_h     {1.67262189800000000000000000000000000000000000000000000000000E-27L};
    // Molar mass of hydrogen in SI (kg.mol^-1)
    static constexpr Real M_h     {m_h * Na};
    // Specific gas constant of hydrogen in SI (J.kg^-1.K^-1)
    static constexpr Real Rstar_h {R / M_h};
};
#elif defined(FLOAT64)
class SI_Constants
{
public:
    // Avogadro number in SI (mol^-1)
    static constexpr Real Na      {6.02214085700000000000000000000000000000000000000000000000000E+23};
    // Boltzmann constant in SI (J.K^-1)
    static constexpr Real k_b     {1.38064852000000000000000000000000000000000000000000000000000E-23};
    // Universal gas constant in SI (J.mol^-1.K^-1)
    static constexpr Real R       {8.31445980000000000000000000000000000000000000000000000000000E+00};
    // Mean molecular weight of hydrogen (no unit)
    static constexpr Real mmw_h   {constants::one};
    // Mass of hydrogen in SI (kg)
    static constexpr Real m_h     {1.67262189800000000000000000000000000000000000000000000000000E-27};
    // Molar mass of hydrogen in SI (kg.mol^-1)
    static constexpr Real M_h     {m_h * Na};
    // Specific gas constant of hydrogen in SI (J.kg^-1.K^-1)
    static constexpr Real Rstar_h {R / M_h};
};
#elif defined(FLOAT32)
class SI_Constants
{
public:
    // Avogadro number in SI (mol^-1)
    static constexpr Real Na      {6.02214085700000000000000000000000000000000000000000000000000E+23F};
    // Boltzmann constant in SI (J.K^-1)
    static constexpr Real k_b     {1.38064852000000000000000000000000000000000000000000000000000E-23F};
    // Universal gas constant in SI (J.mol^-1.K^-1)
    static constexpr Real R       {8.31445980000000000000000000000000000000000000000000000000000E+00F};
    // Mean molecular weight of hydrogen (no unit)
    static constexpr Real mmw_h   {constants::one};
    // Mass of hydrogen in SI (kg)
    static constexpr Real m_h     {1.67262189800000000000000000000000000000000000000000000000000E-27F};
    // Molar mass of hydrogen in SI (kg.mol^-1)
    static constexpr Real M_h     {m_h * Na};
    // Specific gas constant of hydrogen in SI (J.kg^-1.K^-1)
    static constexpr Real Rstar_h {R / M_h};
};
#else
    static_assert(false, "You must define a floating point type")
#endif

#if defined(CGS_UNITS)
static constexpr Real length          {constants::hundredth}; // in SI (m)
static constexpr Real mass            {constants::thousandth}; // in SI (kg
static constexpr Real time            {constants::one}; // in SI (s)
static constexpr Real temperature     {constants::one}; // in SI (K)
static constexpr Real matter_quantity {constants::one}; // in SI (mol)
#elif defined(SI_UNITS)
static constexpr Real length          {constants::one}; // in SI (m)
static constexpr Real mass            {constants::one}; // in SI (kg)
static constexpr Real time            {constants::one}; // in SI (s)
static constexpr Real temperature     {constants::one}; // in SI (K)
static constexpr Real matter_quantity {constants::one}; // in SI (mol)
#else
static_assert(false, "You must define a system unit !");
#endif

static constexpr Real volume          {length * length * length};
static constexpr Real density         {mass / volume};
static constexpr Real velocity        {length / time};
static constexpr Real energy          {mass * velocity * velocity};
static constexpr Real pressure        {energy / volume};

namespace constants
{
static constexpr Real Na              {SI_Constants::Na * matter_quantity};
static constexpr Real k_b             {SI_Constants::k_b / (energy / temperature)};
static constexpr Real R               {SI_Constants::R   / (energy / (matter_quantity * temperature))};
static constexpr Real mmw_h           {SI_Constants::mmw_h};
static constexpr Real m_h             {SI_Constants::m_h / mass};
static constexpr Real M_h             {m_h * Na};
static constexpr Real Rstar_h         {R / M_h};
}

}}
