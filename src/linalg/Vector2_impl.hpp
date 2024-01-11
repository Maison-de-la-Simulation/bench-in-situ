#pragma once

#include "Vector.hpp"

namespace hydro
{

// Member functions
// ----------------
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<2, T>::dim() noexcept
{
    return 2;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<2, T>::rank() noexcept
{
    return 1;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<2, T>::size() noexcept
{
    return 2;
}

// Constructors / destructor
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>::Vector(const T& x1, const T& x2) noexcept
{
    this->m_data[0] = x1;
    this->m_data[1] = x2;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>::Vector(const Vector<2, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>::Vector(Vector<2, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>& Vector<2, T>::operator=(const Vector<2, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>& Vector<2, T>::operator=(Vector<2, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    return *this;
}

// Accessing operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T& Vector<2, T>::operator()(int i) noexcept
{
    return m_data[i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T Vector<2, T>::operator()(int i) const noexcept
{
    return m_data[i];
}

template <class T>
template <int i>
KOKKOS_FORCEINLINE_FUNCTION
T& Vector<2, T>::operator()(std::integral_constant<int, i>) noexcept
{
    return m_data[i];
}

template <class T>
template <int i>
KOKKOS_FORCEINLINE_FUNCTION
T Vector<2, T>::operator()(std::integral_constant<int, i>) const noexcept
{
    return m_data[i];
}

// Linear algebra operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>& Vector<2, T>::operator+=(const Vector<2, T>& rhs) noexcept
{
    this->m_data[0] += rhs.m_data[0];
    this->m_data[1] += rhs.m_data[1];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>& Vector<2, T>::operator-=(const Vector<2, T>& rhs) noexcept
{
    this->m_data[0] -= rhs.m_data[0];
    this->m_data[1] -= rhs.m_data[1];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T>& Vector<2, T>::operator*=(T val) noexcept
{
    this->m_data[0] *= val;
    this->m_data[1] *= val;
    return *this;
}

// Friend functions
// ----------------
// Linear algebra operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T> operator+(const Vector<2, T>& lhs, const Vector<2, T>& rhs) noexcept
{
    Vector<2, T> result;
    result.m_data[0] = lhs.m_data[0] + rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] + rhs.m_data[1];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T> operator-(const Vector<2, T>& lhs, const Vector<2, T>& rhs) noexcept
{
    Vector<2, T> result;
    result.m_data[0] = lhs.m_data[0] - rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] - rhs.m_data[1];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, T> operator*(const T& val, const Vector<2, T>& rhs) noexcept
{
    Vector<2, T> result;
    result.m_data[0] = val * rhs.m_data[0];
    result.m_data[1] = val * rhs.m_data[1];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T dot(const Vector<2, T>& x, const Vector<2, T>& y) noexcept
{
    return x.m_data[0] * y.m_data[0] + x.m_data[1] * y.m_data[1];
}

template <class U>
KOKKOS_FORCEINLINE_FUNCTION
void gradient(const Vector<2, U>& x, const Vector<2, U>& y, const U& invdl, Vector<2, U>& grad) noexcept
{
    grad.m_data[0] = (y.m_data[0] - x.m_data[0]) * invdl;
    grad.m_data[1] = (y.m_data[1] - x.m_data[1]) * invdl;
}

template <class U>
KOKKOS_FORCEINLINE_FUNCTION
Vector<2, U> mean(const Vector<2, U>& x, const Vector<2, U>& y) noexcept
{
    Vector<2, U> mean;
    mean.m_data[0] = 0.5 * (x.m_data[0] + y.m_data[0]);
    mean.m_data[1] = 0.5 * (x.m_data[1] + y.m_data[1]);
    return mean;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<2, T>& rhs)
{
    os << "[" << rhs.m_data[0] << "," << rhs.m_data[1] << "]";
    return os;
}

}  // hydro
