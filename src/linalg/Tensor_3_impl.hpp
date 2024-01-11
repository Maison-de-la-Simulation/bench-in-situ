#pragma once

#include "Tensor.hpp"
#include "Vector.hpp"

namespace hydro
{

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>::Tensor(const T& x00, const T& x01, const T& x02,
                     const T& x10, const T& x11, const T& x12,
                     const T& x20, const T& x21, const T& x22) noexcept
{
    this->m_data[0] = x00;
    this->m_data[1] = x01;
    this->m_data[2] = x02;
    this->m_data[3] = x10;
    this->m_data[4] = x11;
    this->m_data[5] = x12;
    this->m_data[6] = x20;
    this->m_data[7] = x21;
    this->m_data[8] = x22;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>::Tensor(const Tensor<3, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
    this->m_data[4] = x.m_data[4];
    this->m_data[5] = x.m_data[5];
    this->m_data[6] = x.m_data[6];
    this->m_data[7] = x.m_data[7];
    this->m_data[8] = x.m_data[8];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>::Tensor(Tensor<3, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
    this->m_data[4] = x.m_data[4];
    this->m_data[5] = x.m_data[5];
    this->m_data[6] = x.m_data[6];
    this->m_data[7] = x.m_data[7];
    this->m_data[8] = x.m_data[8];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>& Tensor<3, T>::operator=(const Tensor<3, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
    this->m_data[4] = x.m_data[4];
    this->m_data[5] = x.m_data[5];
    this->m_data[6] = x.m_data[6];
    this->m_data[7] = x.m_data[7];
    this->m_data[8] = x.m_data[8];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>& Tensor<3, T>::operator=(Tensor<3, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
    this->m_data[4] = x.m_data[4];
    this->m_data[5] = x.m_data[5];
    this->m_data[6] = x.m_data[6];
    this->m_data[7] = x.m_data[7];
    this->m_data[8] = x.m_data[8];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Tensor<3, T>::size() noexcept
{
    return 3 * 3;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T& Tensor<3, T>::operator()(int i, int j) noexcept
{
    return m_data[j + 3*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T Tensor<3, T>::operator()(int i, int j) const noexcept
{
    return m_data[j + 3*i];
}

template <class T>
template <int i, int j>
KOKKOS_FORCEINLINE_FUNCTION
T& Tensor<3, T>::operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) noexcept
{
    return m_data[j + 3*i];
}

template <class T>
template <int i, int j>
KOKKOS_FORCEINLINE_FUNCTION
T Tensor<3, T>::operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) const noexcept
{
    return m_data[j + 3*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>& Tensor<3, T>::operator+=(const Tensor<3, T>& rhs) noexcept
{
    this->m_data[0] += rhs.m_data[0];
    this->m_data[1] += rhs.m_data[1];
    this->m_data[2] += rhs.m_data[2];
    this->m_data[3] += rhs.m_data[3];
    this->m_data[4] += rhs.m_data[4];
    this->m_data[5] += rhs.m_data[5];
    this->m_data[6] += rhs.m_data[6];
    this->m_data[7] += rhs.m_data[7];
    this->m_data[8] += rhs.m_data[8];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>& Tensor<3, T>::operator-=(const Tensor<3, T>& rhs) noexcept
{
    this->m_data[0] -= rhs.m_data[0];
    this->m_data[1] -= rhs.m_data[1];
    this->m_data[2] -= rhs.m_data[2];
    this->m_data[3] -= rhs.m_data[3];
    this->m_data[4] -= rhs.m_data[4];
    this->m_data[5] -= rhs.m_data[5];
    this->m_data[6] -= rhs.m_data[6];
    this->m_data[7] -= rhs.m_data[7];
    this->m_data[8] -= rhs.m_data[8];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T>& Tensor<3, T>::operator*=(T val) noexcept
{
    this->m_data[0] *= val;
    this->m_data[1] *= val;
    this->m_data[2] *= val;
    this->m_data[3] *= val;
    this->m_data[4] *= val;
    this->m_data[5] *= val;
    this->m_data[6] *= val;
    this->m_data[7] *= val;
    this->m_data[8] *= val;
    return *this;
}

template <int i, int j, class T>
KOKKOS_FORCEINLINE_FUNCTION
T& get(Tensor<3, T>& x) noexcept
{
    return x.m_data[j + 3*i];
}

template <int i, int j, class T>
KOKKOS_FORCEINLINE_FUNCTION
T get(const Tensor<3, T>& x) noexcept
{
    return x.m_data[j + 3*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T> operator+(const Tensor<3, T>& lhs, const Tensor<3, T>& rhs) noexcept
{
    Tensor<3, T> result;
    result.m_data[0] = lhs.m_data[0] + rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] + rhs.m_data[1];
    result.m_data[2] = lhs.m_data[2] + rhs.m_data[2];
    result.m_data[3] = lhs.m_data[3] + rhs.m_data[3];
    result.m_data[4] = lhs.m_data[4] + rhs.m_data[4];
    result.m_data[5] = lhs.m_data[5] + rhs.m_data[5];
    result.m_data[6] = lhs.m_data[6] + rhs.m_data[6];
    result.m_data[7] = lhs.m_data[7] + rhs.m_data[7];
    result.m_data[8] = lhs.m_data[8] + rhs.m_data[8];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T> operator-(const Tensor<3, T>& lhs, const Tensor<3, T>& rhs) noexcept
{
    Tensor<3, T> result;
    result.m_data[0] = lhs.m_data[0] - rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] - rhs.m_data[1];
    result.m_data[2] = lhs.m_data[2] - rhs.m_data[2];
    result.m_data[3] = lhs.m_data[3] - rhs.m_data[3];
    result.m_data[4] = lhs.m_data[4] - rhs.m_data[4];
    result.m_data[5] = lhs.m_data[5] - rhs.m_data[5];
    result.m_data[6] = lhs.m_data[6] - rhs.m_data[6];
    result.m_data[7] = lhs.m_data[7] - rhs.m_data[7];
    result.m_data[8] = lhs.m_data[8] - rhs.m_data[8];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<3, T> operator*(const T& val, const Tensor<3, T>& rhs) noexcept
{
    Tensor<3, T> result;
    result.m_data[0] = val * rhs.m_data[0];
    result.m_data[1] = val * rhs.m_data[1];
    result.m_data[2] = val * rhs.m_data[2];
    result.m_data[3] = val * rhs.m_data[3];
    result.m_data[4] = val * rhs.m_data[4];
    result.m_data[5] = val * rhs.m_data[5];
    result.m_data[6] = val * rhs.m_data[6];
    result.m_data[7] = val * rhs.m_data[7];
    result.m_data[8] = val * rhs.m_data[8];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T trace(const Tensor<3, T>& x) noexcept
{
    return x.m_data[0] + x.m_data[4] + x.m_data[8];
}

template <int i, class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T> sym(std::integral_constant<int, i>, const Tensor<3, T>& x) noexcept
{
    Vector<3, T> sym;
    sym.m_data[0] = 0.5 * (x.m_data[i + 3*0] + x.m_data[0 + 3*i]);
    sym.m_data[1] = 0.5 * (x.m_data[i + 3*1] + x.m_data[1 + 3*i]);
    sym.m_data[2] = 0.5 * (x.m_data[i + 3*2] + x.m_data[2 + 3*i]);
    return sym;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Tensor<3, T>& rhs)
{
    os << "[";
    os << "[" << get<0, 0>(rhs) << "," << get<0, 1>(rhs) << "," << get<0, 2>(rhs) << "]";
    os << ",";
    os << "[" << get<1, 0>(rhs) << "," << get<1, 1>(rhs) << "," << get<1, 2>(rhs) << "]";
    os << ",";
    os << "[" << get<2, 0>(rhs) << "," << get<2, 1>(rhs) << "," << get<2, 2>(rhs) << "]";
    os << "]";
    return os;
}

}  // hydro
