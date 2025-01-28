#pragma once

#include "Vector.hpp"

namespace hydro
{

// Member functions
// ----------------
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<1, T>::dim() noexcept
{
    return 1;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<1, T>::rank() noexcept
{
    return 1;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<1, T>::size() noexcept
{
    return 1;
}

// Constructors / destructor
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>::Vector(const T& x1) noexcept
{
    this->m_data[0] = x1;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>::Vector(const Vector<1, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>::Vector(Vector<1, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>& Vector<1, T>::operator=(const Vector<1, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>& Vector<1, T>::operator=(Vector<1, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    return *this;
}

// Accessing operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T& Vector<1, T>::operator()(int i) noexcept
{
    return m_data[i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T Vector<1, T>::operator()(int i) const noexcept
{
    return m_data[i];
}

template <class T>
template <int i>
KOKKOS_FORCEINLINE_FUNCTION
T& Vector<1, T>::operator()(std::integral_constant<int, i>) noexcept
{
    return m_data[i];
}

template <class T>
template <int i>
KOKKOS_FORCEINLINE_FUNCTION
T Vector<1, T>::operator()(std::integral_constant<int, i>) const noexcept
{
    return m_data[i];
}

// Linear algebra operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>& Vector<1, T>::operator+=(const Vector<1, T>& rhs) noexcept
{
    this->m_data[0] += rhs.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>& Vector<1, T>::operator-=(const Vector<1, T>& rhs) noexcept
{
    this->m_data[0] -= rhs.m_data[0];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T>& Vector<1, T>::operator*=(T val) noexcept
{
    this->m_data[0] *= val;
    return *this;
}

// Friend functions
// ----------------
// Linear algebra operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T> operator+(const Vector<1, T>& lhs, const Vector<1, T>& rhs) noexcept
{
    Vector<1, T> result;
    result.m_data[0] = lhs.m_data[0] + rhs.m_data[0];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T> operator-(const Vector<1, T>& lhs, const Vector<1, T>& rhs) noexcept
{
    Vector<1, T> result;
    result.m_data[0] = lhs.m_data[0] - rhs.m_data[0];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, T> operator*(const T& val, const Vector<1, T>& rhs) noexcept
{
    Vector<1, T> result;
    result.m_data[0] = val * rhs.m_data[0];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T dot(const Vector<1, T>& x, const Vector<1, T>& y) noexcept
{
    return x.m_data[0] * y.m_data[0];
}

template <class U>
KOKKOS_FORCEINLINE_FUNCTION
void gradient(const Vector<1, U>& x, const Vector<1, U>& y, const U invdl, Vector<1, U>& grad) noexcept
{
    grad.m_data[0] = (y.m_data[0] - x.m_data[0]) * invdl;
}

template <class U>
KOKKOS_FORCEINLINE_FUNCTION
Vector<1, U> mean(const Vector<1, U>& x, const Vector<1, U>& y) noexcept
{
    Vector<1, U> mean;
    mean.m_data[0] = 0.5 * (x.m_data[0] + y.m_data[0]);
    return mean;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<1, T>& rhs)
{
    os << "[" << rhs.m_data[0] << "]";
    return os;
}

}  // hydro
