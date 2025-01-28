#pragma once

#include "Tensor.hpp"
#include "Vector.hpp"

namespace hydro
{

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>::Tensor(const T& x00, const T& x01,
                     const T& x10, const T& x11) noexcept
{
    this->m_data[0] = x00;
    this->m_data[1] = x01;
    this->m_data[2] = x10;
    this->m_data[3] = x11;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>::Tensor(const Tensor<2, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>::Tensor(Tensor<2, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>& Tensor<2, T>::operator=(const Tensor<2, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>& Tensor<2, T>::operator=(Tensor<2, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    this->m_data[3] = x.m_data[3];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Tensor<2, T>::size() noexcept
{
    return 2 * 2;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T& Tensor<2, T>::operator()(int i, int j) noexcept
{
    return m_data[j + 2*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T Tensor<2, T>::operator()(int i, int j) const noexcept
{
    return m_data[j + 2*i];
}

template <class T>
template <int i, int j>
KOKKOS_FORCEINLINE_FUNCTION
T& Tensor<2, T>::operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) noexcept
{
    return m_data[j + 2*i];
}

template <class T>
template <int i, int j>
KOKKOS_FORCEINLINE_FUNCTION
T Tensor<2, T>::operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) const noexcept
{
    return m_data[j + 2*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>& Tensor<2, T>::operator+=(const Tensor<2, T>& rhs) noexcept
{
    this->m_data[0] += rhs.m_data[0];
    this->m_data[1] += rhs.m_data[1];
    this->m_data[2] += rhs.m_data[2];
    this->m_data[3] += rhs.m_data[3];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>& Tensor<2, T>::operator-=(const Tensor<2, T>& rhs) noexcept
{
    this->m_data[0] -= rhs.m_data[0];
    this->m_data[1] -= rhs.m_data[1];
    this->m_data[2] -= rhs.m_data[2];
    this->m_data[3] -= rhs.m_data[3];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T>& Tensor<2, T>::operator*=(T val) noexcept
{
    this->m_data[0] *= val;
    this->m_data[1] *= val;
    this->m_data[2] *= val;
    this->m_data[3] *= val;
    return *this;
}

template <int i, int j, class T>
KOKKOS_FORCEINLINE_FUNCTION
T& get(Tensor<2, T>& x) noexcept
{
    return x.m_data[j + 2*i];
}

template <int i, int j, class T>
KOKKOS_FORCEINLINE_FUNCTION
T get(const Tensor<2, T>& x) noexcept
{
    return x.m_data[j + 2*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T> operator+(const Tensor<2, T>& lhs, const Tensor<2, T>& rhs) noexcept
{
    Tensor<2, T> result;
    result.m_data[0] = lhs.m_data[0] + rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] + rhs.m_data[1];
    result.m_data[2] = lhs.m_data[2] + rhs.m_data[2];
    result.m_data[3] = lhs.m_data[3] + rhs.m_data[3];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T> operator-(const Tensor<2, T>& lhs, const Tensor<2, T>& rhs) noexcept
{
    Tensor<2, T> result;
    result.m_data[0] = lhs.m_data[0] - rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] - rhs.m_data[1];
    result.m_data[2] = lhs.m_data[2] - rhs.m_data[2];
    result.m_data[3] = lhs.m_data[3] - rhs.m_data[3];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<2, T> operator*(const T& val, const Tensor<2, T>& rhs) noexcept
{
    Tensor<2, T> result;
    result.m_data[0] = val * rhs.m_data[0];
    result.m_data[1] = val * rhs.m_data[1];
    result.m_data[2] = val * rhs.m_data[2];
    result.m_data[3] = val * rhs.m_data[3];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T trace(const Tensor<2, T>& x) noexcept
{
    return x.m_data[0] + x.m_data[3];
}

template <int i, class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T> sym(std::integral_constant<int, i>, const Tensor<2, T>& x) noexcept
{
    Vector<2, T> sym;
    sym.m_data[0] = 0.5 * (x.m_data[i + 2*0] + x.m_data[0 + 2*i]);
    sym.m_data[1] = 0.5 * (x.m_data[i + 2*1] + x.m_data[1 + 2*i]);
    return sym;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Tensor<2, T>& rhs)
{
    os << "[";
    os << "[" << get<0, 0>(rhs) << "," << get<0, 1>(rhs) << "]";
    os << ",";
    os << "[" << get<1, 0>(rhs) << "," << get<1, 1>(rhs) << "]";
    os << "]";
    return os;
}

}  // hydro
