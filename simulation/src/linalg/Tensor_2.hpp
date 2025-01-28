#pragma once

#include "Tensor.hpp"
#include "Vector.hpp"

#include <iostream>
#include <Kokkos_Core.hpp>

namespace hydro
{

template <class T>
class Tensor<2, T>
{
public:
    KOKKOS_DEFAULTED_FUNCTION
    Tensor<2, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>(const T& x00, const T& x01,
                 const T& x10, const T& x11) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>(const Tensor<2, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>(Tensor<2, T>&& x) noexcept;
    KOKKOS_DEFAULTED_FUNCTION
    ~Tensor<2, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>& operator=(const Tensor<2, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>& operator=(Tensor<2, T>&& x) noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    static constexpr int size() noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    T& operator()(int i, int j) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    T operator()(int i, int j) const noexcept;
    template <int i, int j>
    KOKKOS_FORCEINLINE_FUNCTION
    T& operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) noexcept;
    template <int i, int j>
    KOKKOS_FORCEINLINE_FUNCTION
    T operator()(std::integral_constant<int, i>, std::integral_constant<int, j>) const noexcept;

    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>& operator+=(const Tensor<2, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>& operator-=(const Tensor<2, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<2, T>& operator*=(T val) noexcept;

    template <int i, int j, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U& get(Tensor<2, U>& x) noexcept;
    template <int i, int j, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U get(const Tensor<2, U>& x) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<2, U> operator+(const Tensor<2, U>& lhs, const Tensor<2, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<2, U> operator-(const Tensor<2, U>& lhs, const Tensor<2, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<2, U> operator*(const U& val, const Tensor<2, U>& rhs) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U trace(const Tensor<2, U>& x) noexcept;
    template <int i, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<2, U> sym(std::integral_constant<int, i>, const Tensor<2, U>& x) noexcept;

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const Tensor<2, U>& rhs);

public:
    T m_data[4];
};

}  // hydro

#include "Tensor_2_impl.hpp"
