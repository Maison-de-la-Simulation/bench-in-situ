#pragma once

#include "Tensor.hpp"
#include "Vector.hpp"

#include <iostream>
#include <Kokkos_Core.hpp>

namespace hydro
{

template <class T>
class Tensor<1, T>
{
public:
    KOKKOS_DEFAULTED_FUNCTION
    Tensor<1, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>(const T& x00) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>(const Tensor<1, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>(Tensor<1, T>&& x) noexcept;
    KOKKOS_DEFAULTED_FUNCTION
    ~Tensor<1, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>& operator=(const Tensor<1, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>& operator=(Tensor<1, T>&& x) noexcept;

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
    Tensor<1, T>& operator+=(const Tensor<1, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>& operator-=(const Tensor<1, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<1, T>& operator*=(T val) noexcept;

    template <int i, int j, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U& get(Tensor<1, U>& x) noexcept;
    template <int i, int j, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U get(const Tensor<1, U>& x) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<1, U> operator+(const Tensor<1, U>& lhs, const Tensor<1, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<1, U> operator-(const Tensor<1, U>& lhs, const Tensor<1, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<1, U> operator*(const U& val, const Tensor<1, U>& rhs) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U trace(const Tensor<1, U>& x) noexcept;
    template <int i, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<1, U> sym(std::integral_constant<int, i>, const Tensor<1, U>& x) noexcept;

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const Tensor<1, U>& rhs);

public:
    T m_data[1];
};

}  // hydro

#include "Tensor_1_impl.hpp"
