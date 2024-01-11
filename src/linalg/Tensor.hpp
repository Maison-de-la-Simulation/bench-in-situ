#pragma once

namespace hydro
{

template <int dim, class T>
class Tensor;

using Tensor_1f = Tensor<1, float>;
using Tensor_2f = Tensor<2, float>;
using Tensor_3f = Tensor<3, float>;

using Tensor_1d = Tensor<1, double>;
using Tensor_2d = Tensor<2, double>;
using Tensor_3d = Tensor<3, double>;

}  // hydro

#include "Tensor_1.hpp"
#include "Tensor_2.hpp"
#include "Tensor_3.hpp"
