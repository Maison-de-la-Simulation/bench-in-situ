#pragma once

#include "Vector.hpp"

namespace hydro
{

// Member functions
// ----------------
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<3, T>::dim() noexcept
{
    return 3;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<3, T>::rank() noexcept
{
    return 1;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
constexpr int Vector<3, T>::size() noexcept
{
    return 3;
}

// Constructors / destructor
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>::Vector(const T& x1, const T& x2, const T& x3) noexcept
{
    this->m_data[0] = x1;
    this->m_data[1] = x2;
    this->m_data[2] = x3;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>::Vector(const Vector<3, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>::Vector(Vector<3, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>& Vector<3, T>::operator=(const Vector<3, T>& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>& Vector<3, T>::operator=(Vector<3, T>&& x) noexcept
{
    this->m_data[0] = x.m_data[0];
    this->m_data[1] = x.m_data[1];
    this->m_data[2] = x.m_data[2];
    return *this;
}

// Accessing operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T& Vector<3, T>::operator()(int i) noexcept
{
    return m_data[i];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T Vector<3, T>::operator()(int i) const noexcept
{
    return m_data[i];
}

template <class T>
template <int i>
KOKKOS_FORCEINLINE_FUNCTION
T& Vector<3, T>::operator()(std::integral_constant<int, i>) noexcept
{
    return m_data[i];
}

template <class T>
template <int i>
KOKKOS_FORCEINLINE_FUNCTION
T Vector<3, T>::operator()(std::integral_constant<int, i>) const noexcept
{
    return m_data[i];
}

// Linear algebra operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>& Vector<3, T>::operator+=(const Vector<3, T>& rhs) noexcept
{
    this->m_data[0] += rhs.m_data[0];
    this->m_data[1] += rhs.m_data[1];
    this->m_data[2] += rhs.m_data[2];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>& Vector<3, T>::operator-=(const Vector<3, T>& rhs) noexcept
{
    this->m_data[0] -= rhs.m_data[0];
    this->m_data[1] -= rhs.m_data[1];
    this->m_data[2] -= rhs.m_data[2];
    return *this;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T>& Vector<3, T>::operator*=(T val) noexcept
{
    this->m_data[0] *= val;
    this->m_data[1] *= val;
    this->m_data[2] *= val;
    return *this;
}

// Friend functions
// ----------------
// Linear algebra operators
template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T> operator+(const Vector<3, T>& lhs, const Vector<3, T>& rhs) noexcept
{
    Vector<3, T> result;
    result.m_data[0] = lhs.m_data[0] + rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] + rhs.m_data[1];
    result.m_data[2] = lhs.m_data[2] + rhs.m_data[2];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T> operator-(const Vector<3, T>& lhs, const Vector<3, T>& rhs) noexcept
{
    Vector<3, T> result;
    result.m_data[0] = lhs.m_data[0] - rhs.m_data[0];
    result.m_data[1] = lhs.m_data[1] - rhs.m_data[1];
    result.m_data[2] = lhs.m_data[2] - rhs.m_data[2];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, T> operator*(const T& val, const Vector<3, T>& rhs) noexcept
{
    Vector<3, T> result;
    result.m_data[0] = val * rhs.m_data[0];
    result.m_data[1] = val * rhs.m_data[1];
    result.m_data[2] = val * rhs.m_data[2];
    return result;
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION
T dot(const Vector<3, T>& x, const Vector<3, T>& y) noexcept
{
    return x.m_data[0] * y.m_data[0] + x.m_data[1] * y.m_data[1] + x.m_data[2] * y.m_data[2];
}

template <class U>
KOKKOS_FORCEINLINE_FUNCTION
void gradient(const Vector<3, U>& x, const Vector<3, U>& y, const U& invdl, Vector<3, U>& grad) noexcept
{
    grad.m_data[0] = (y.m_data[0] - x.m_data[0]) * invdl;
    grad.m_data[1] = (y.m_data[1] - x.m_data[1]) * invdl;
    grad.m_data[2] = (y.m_data[2] - x.m_data[2]) * invdl;
}

template <class U>
KOKKOS_FORCEINLINE_FUNCTION
Vector<3, U> mean(const Vector<3, U>& x, const Vector<3, U>& y) noexcept
{
    Vector<3, U> mean;
    mean.m_data[0] = 0.5 * (x.m_data[0] + y.m_data[0]);
    mean.m_data[1] = 0.5 * (x.m_data[1] + y.m_data[1]);
    mean.m_data[2] = 0.5 * (x.m_data[2] + y.m_data[2]);
    return mean;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<3, T>& rhs)
{
    os << "[" << rhs.m_data[0] << "," << rhs.m_data[1] << "," << rhs.m_data[2] << "]";
    return os;
}

}  // hydro
