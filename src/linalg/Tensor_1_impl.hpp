#pragma once

#include "Tensor.hpp"
#include "Vector.hpp"

namespace hydro
{

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>::Tensor(const T& x00) noexcept
{
    this->m_data[0] = x00;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>::Tensor(const Tensor<1, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>::Tensor(Tensor<1, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>& Tensor<1, T>::operator=(const Tensor<1, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>& Tensor<1, T>::operator=(Tensor<1, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Tensor<1, T>::size() noexcept
{
    return 1 * 1;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T& Tensor<1, T>::operator()(int i, int j) noexcept
{
    return m_data[j + 0*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T Tensor<1, T>::operator()(int i, int j) const noexcept
{
    return m_data[j + 0*i];
}

template <class T>
template <int i, int j>
KOKKOS_FORCEINLINE_FUNCTION
T& Tensor<1, T>::operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) noexcept
{
    return m_data[j + 0*i];
}

template <class T>
template <int i, int j>
KOKKOS_FORCEINLINE_FUNCTION
T Tensor<1, T>::operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) const noexcept
{
    return m_data[j + 0*i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>& Tensor<1, T>::operator+=(const Tensor<1, T>& rhs) noexcept
{
    this->m_data[0] += rhs.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>& Tensor<1, T>::operator-=(const Tensor<1, T>& rhs) noexcept
{
    this->m_data[0] -= rhs.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T>& Tensor<1, T>::operator*=(T val) noexcept
{
    this->m_data[0] *= val;
    return *this;
}

template <int i, int j, class T>
KOKKOS_FORCEINLINE_FUNCTION
T& get(Tensor<1, T>& x) noexcept
{
    return x.m_data[0];
}

template <int i, int j, class T>
KOKKOS_FORCEINLINE_FUNCTION
T get(const Tensor<1, T>& x) noexcept
{
    return x.m_data[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T> operator+(const Tensor<1, T>& lhs, const Tensor<1, T>& rhs) noexcept
{
    Tensor<1, T> result;
    result.m_data[0] = lhs.m_data[0] + rhs.m_data[0];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T> operator-(const Tensor<1, T>& lhs, const Tensor<1, T>& rhs) noexcept
{
    Tensor<1, T> result;
    result.m_data[0] = lhs.m_data[0] - rhs.m_data[0];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Tensor<1, T> operator*(const T& val, const Tensor<1, T>& rhs) noexcept
{
    Tensor<1, T> result;
    result.m_data[0] = val * rhs.m_data[0];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T trace(const Tensor<1, T>& x) noexcept
{
    return x.m_data[0];
}

template <int i, class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T> sym(std::integral_constant<int, i>, const Tensor<1, T>& x) noexcept
{
    Vector<1, T> sym;
    sym.m_data[0] = 0.5 * (x.m_data[i] + x.m_data[0]);
    return sym;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Tensor<1, T>& rhs)
{
    os << "[" << "[" << get<0, 0>(rhs) << "]" << "]";
    return os;
}

}  // hydro
