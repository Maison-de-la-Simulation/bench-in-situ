#pragma once

#include "Vector.hpp"

#include <iostream>
#include <Kokkos_Core.hpp>
#include <type_traits>

namespace hydro
{

template <class T>
class Vector<3, T>
{
public:
    KOKKOS_FORCEINLINE_FUNCTION
    static constexpr int dim() noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    static constexpr int rank() noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    static constexpr int size() noexcept;

    // Member functions
    // ----------------
public:
    // Constructors / destructor
    KOKKOS_DEFAULTED_FUNCTION
    Vector<3, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>(const T& x1, const T& x2, const T& x3) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>(const Vector<3, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>(Vector<3, T>&& x) noexcept;
    KOKKOS_DEFAULTED_FUNCTION
    ~Vector<3, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>& operator=(const Vector<3, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>& operator=(Vector<3, T>&& x) noexcept;

    // Accessing operators
    KOKKOS_FORCEINLINE_FUNCTION
    T& operator()(int i) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    T operator()(int i) const noexcept;
    template <int i>
    KOKKOS_FORCEINLINE_FUNCTION
    T& operator()(std::integral_constant<int, i>) noexcept;
    template <int i>
    KOKKOS_FORCEINLINE_FUNCTION
    T operator()(std::integral_constant<int, i>) const noexcept;

    // Linear algebra operators
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>& operator+=(const Vector<3, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>& operator-=(const Vector<3, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<3, T>& operator*=(T val) noexcept;

    // Friend functions
    // ----------------
public:
    // Linear algebra operators
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<3, U> operator+(const Vector<3, U>& lhs, const Vector<3, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<3, U> operator-(const Vector<3, U>& lhs, const Vector<3, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<3, U> operator*(const U& val, const Vector<3, U>& rhs) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U dot(const Vector<3, U>& x, const Vector<3, U>& y) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend void gradient(const Vector<3, U>& x, const Vector<3, U>& y, const U& invdl, Vector<3, U>& grad) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<3, U> mean(const Vector<3, U>& x, const Vector<3, U>& y) noexcept;

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const Vector<3, U>& rhs);

public:
    T m_data[3];
};

}  // hydro

#include "Vector3_impl.hpp"
