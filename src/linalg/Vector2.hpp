#pragma once

#include "Vector.hpp"

#include <iostream>
#include <Kokkos_Core.hpp>
#include <type_traits>

namespace hydro
{

template <class T>
class Vector<2, T>
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
    Vector<2, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>(const T& x1, const T& x2) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>(const Vector<2, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>(Vector<2, T>&& x) noexcept;
    KOKKOS_DEFAULTED_FUNCTION
    ~Vector<2, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>& operator=(const Vector<2, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>& operator=(Vector<2, T>&& x) noexcept;

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
    Vector<2, T>& operator+=(const Vector<2, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>& operator-=(const Vector<2, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<2, T>& operator*=(T val) noexcept;

    // Friend functions
    // ----------------
public:
    // Linear algebra operators
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<2, U> operator+(const Vector<2, U>& lhs, const Vector<2, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<2, U> operator-(const Vector<2, U>& lhs, const Vector<2, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<2, U> operator*(const U& val, const Vector<2, U>& rhs) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U dot(const Vector<2, U>& x, const Vector<2, U>& y) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend void gradient(const Vector<2, U>& x, const Vector<2, U>& y, const U& invdl, Vector<2, U>& grad) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<2, U> mean(const Vector<2, U>& x, const Vector<2, U>& y) noexcept;

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const Vector<2, U>& rhs);

public:
    T m_data[2];
};

}  // hydro

#include "Vector2_impl.hpp"
