#pragma once

#include "Tensor.hpp"
#include "Vector.hpp"

#include <iostream>
#include <Kokkos_Core.hpp>

namespace hydro
{

template <class T>
class Tensor<3, T>
{
public:
    KOKKOS_DEFAULTED_FUNCTION
    Tensor<3, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>(const T& x00, const T& x01, const T& x02,
                 const T& x10, const T& x11, const T& x12,
                 const T& x20, const T& x21, const T& x22) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>(const Tensor<3, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>(Tensor<3, T>&& x) noexcept;
    KOKKOS_DEFAULTED_FUNCTION
    ~Tensor<3, T>() noexcept = default;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>& operator=(const Tensor<3, T>& x) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>& operator=(Tensor<3, T>&& x) noexcept;

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
    Tensor<3, T>& operator+=(const Tensor<3, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>& operator-=(const Tensor<3, T>& rhs) noexcept;
    KOKKOS_FORCEINLINE_FUNCTION
    Tensor<3, T>& operator*=(T val) noexcept;

    template <int i, int j, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U& get(Tensor<3, U>& x) noexcept;
    template <int i, int j, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U get(const Tensor<3, U>& x) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<3, U> operator+(const Tensor<3, U>& lhs, const Tensor<3, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<3, U> operator-(const Tensor<3, U>& lhs, const Tensor<3, U>& rhs) noexcept;
    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Tensor<3, U> operator*(const U& val, const Tensor<3, U>& rhs) noexcept;

    template <class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend U trace(const Tensor<3, U>& x) noexcept;
    template <int i, class U>
    KOKKOS_FORCEINLINE_FUNCTION
    friend Vector<3, U> sym(std::integral_constant<int, i>, const Tensor<3, U>& x) noexcept;

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const Tensor<3, U>& rhs);

public:
    T m_data[9];
};

}  // hydro

#include "Tensor_3_impl.hpp"
