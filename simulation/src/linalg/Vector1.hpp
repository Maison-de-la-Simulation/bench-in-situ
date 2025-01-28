#pragma once

#include "Vector.hpp"

#include <iostream>
#include <Kokkos_Core.hpp>
#include <type_traits>

namespace hydro
{

template <class T>
class Vector<1, T>
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
    Vector<1, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>(const T& x1) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>(const Vector<1, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>(Vector<1, T>&& x) noexcept;
    KOKKOS_DEFAULTED_FUNCTION
    ~Vector<1, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>& operator=(const Vector<1, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>& operator=(Vector<1, T>&& x) noexcept;

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
    Vector<1, T>& operator+=(const Vector<1, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>& operator-=(const Vector<1, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Vector<1, T>& operator*=(T val) noexcept;

    // Friend functions
    // ----------------
public:
    // Linear algebra operators
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<1, U> operator+(const Vector<1, U>& lhs, const Vector<1, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<1, U> operator-(const Vector<1, U>& lhs, const Vector<1, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<1, U> operator*(const U& val, const Vector<1, U>& rhs) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U dot(const Vector<1, U>& x, const Vector<1, U>& y) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend void gradient(const Vector<1, U>& x, const Vector<1, U>& y, const U& invdl, Vector<1, U>& grad) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<1, U> mean(const Vector<1, U>& x, const Vector<1, U>& y) noexcept;

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const Vector<1, U>& rhs);

public:
    T m_data[1];
};

}  // hydro

#include "Vector1_impl.hpp"
