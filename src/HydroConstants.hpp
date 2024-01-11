#pragma once

#include "HydroTypes.hpp"

namespace hydro { namespace constants
{

// 30 significant digits

#if defined(FLOAT64_EXTENDED)
// Integers
constexpr Real zero       {0.00000000000000000000000000000E+0L};
constexpr Real one        {1.00000000000000000000000000000E+0L};
constexpr Real two        {2.00000000000000000000000000000E+0L};
constexpr Real three      {3.00000000000000000000000000000E+0L};
constexpr Real four       {4.00000000000000000000000000000E+0L};
constexpr Real five       {5.00000000000000000000000000000E+0L};
constexpr Real six        {6.00000000000000000000000000000E+0L};
constexpr Real seven      {7.00000000000000000000000000000E+0L};
constexpr Real eight      {8.00000000000000000000000000000E+0L};
constexpr Real nine       {9.00000000000000000000000000000E+0L};
constexpr Real ten        {1.00000000000000000000000000000E+1L};
constexpr Real hundred    {1.00000000000000000000000000000E+2L};
constexpr Real thousand   {1.00000000000000000000000000000E+3L};

// Rationals
constexpr Real thousandth {1.00000000000000000000000000000E-3L};
constexpr Real hundredth  {1.00000000000000000000000000000E-2L};
constexpr Real tenth      {1.00000000000000000000000000000E-1L};
constexpr Real eigth      {1.25000000000000000000000000000E-1L};
constexpr Real quarter    {2.50000000000000000000000000000E-1L};
constexpr Real third      {3.33333333333333333333333333333E-1L};
constexpr Real half       {5.00000000000000000000000000000E-1L};
constexpr Real two_third  {6.66666666666666666666666666667E-1L};
constexpr Real five_third {1.66666666666666666666666666667E-1L};

// Irrationals
constexpr Real sqrt2      {1.41421356237309504880168872421E+0L};
constexpr Real phi        {1.61803398874989484820458683437E+0L};
constexpr Real e          {2.71828182845904523536028747135E+0L};
constexpr Real pi         {3.14159265358979323846264338328E+0L};
constexpr Real pi_2       {6.28318530717958647692528676656E+0L};
#elif defined(FLOAT64)
// Integers
constexpr Real zero       {0.00000000000000000000000000000E+0};
constexpr Real one        {1.00000000000000000000000000000E+0};
constexpr Real two        {2.00000000000000000000000000000E+0};
constexpr Real three      {3.00000000000000000000000000000E+0};
constexpr Real four       {4.00000000000000000000000000000E+0};
constexpr Real five       {5.00000000000000000000000000000E+0};
constexpr Real six        {6.00000000000000000000000000000E+0};
constexpr Real seven      {7.00000000000000000000000000000E+0};
constexpr Real eight      {8.00000000000000000000000000000E+0};
constexpr Real nine       {9.00000000000000000000000000000E+0};
constexpr Real ten        {1.00000000000000000000000000000E+1};
constexpr Real hundred    {1.00000000000000000000000000000E+2};
constexpr Real thousand   {1.00000000000000000000000000000E+3};

// Rationals
constexpr Real thousandth {1.00000000000000000000000000000E-3};
constexpr Real hundredth  {1.00000000000000000000000000000E-2};
constexpr Real tenth      {1.00000000000000000000000000000E-1};
constexpr Real eigth      {1.25000000000000000000000000000E-1};
constexpr Real quarter    {2.50000000000000000000000000000E-1};
constexpr Real third      {3.33333333333333333333333333333E-1};
constexpr Real half       {5.00000000000000000000000000000E-1};
constexpr Real two_third  {6.66666666666666666666666666667E-1};
constexpr Real five_third {1.66666666666666666666666666667E-1};

// Irrationals
constexpr Real sqrt2      {1.41421356237309504880168872421E+0};
constexpr Real phi        {1.61803398874989484820458683437E+0};
constexpr Real e          {2.71828182845904523536028747135E+0};
constexpr Real pi         {3.14159265358979323846264338328E+0};
constexpr Real pi_2       {6.28318530717958647692528676656E+0};
#elif defined(FLOAT32)
// Integers
constexpr Real zero       {0.00000000000000000000000000000E+0F};
constexpr Real one        {1.00000000000000000000000000000E+0F};
constexpr Real two        {2.00000000000000000000000000000E+0F};
constexpr Real three      {3.00000000000000000000000000000E+0F};
constexpr Real four       {4.00000000000000000000000000000E+0F};
constexpr Real five       {5.00000000000000000000000000000E+0F};
constexpr Real six        {6.00000000000000000000000000000E+0F};
constexpr Real seven      {7.00000000000000000000000000000E+0F};
constexpr Real eight      {8.00000000000000000000000000000E+0F};
constexpr Real nine       {9.00000000000000000000000000000E+0F};
constexpr Real ten        {1.00000000000000000000000000000E+1F};
constexpr Real hundred    {1.00000000000000000000000000000E+2F};
constexpr Real thousand   {1.00000000000000000000000000000E+3F};

// Rationals
constexpr Real thousandth {1.00000000000000000000000000000E-3F};
constexpr Real hundredth  {1.00000000000000000000000000000E-2F};
constexpr Real tenth      {1.00000000000000000000000000000E-1F};
constexpr Real eigth      {1.25000000000000000000000000000E-1F};
constexpr Real quarter    {2.50000000000000000000000000000E-1F};
constexpr Real third      {3.33333333333333333333333333333E-1F};
constexpr Real half       {5.00000000000000000000000000000E-1F};
constexpr Real two_third  {6.66666666666666666666666666667E-1F};
constexpr Real five_third {1.66666666666666666666666666667E-1F};

// Irrationals
constexpr Real sqrt2      {1.41421356237309504880168872421E+0F};
constexpr Real phi        {1.61803398874989484820458683437E+0F};
constexpr Real e          {2.71828182845904523536028747135E+0F};
constexpr Real pi         {3.14159265358979323846264338328E+0F};
constexpr Real pi_2       {6.28318530717958647692528676656E+0F};
#else
static_assert(false, "You must define a realing point type")
#endif

constexpr Real infinity   {std::numeric_limits<Real>::infinity()};

}}
