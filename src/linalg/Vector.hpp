#pragma once

namespace hydro
{

template <int dim, class T>
class Vector;

using Vector1f = Vector<1, float>;
using Vector2f = Vector<2, float>;
using Vector3f = Vector<3, float>;

using Vector1d = Vector<1, double>;
using Vector2d = Vector<2, double>;
using Vector3d = Vector<3, double>;

}  // hydro

#include "Vector1.hpp"
#include "Vector2.hpp"
#include "Vector3.hpp"
