#pragma once

#include <Kokkos_Core.hpp>
#include <cmath>
#include <limits>

namespace utils
{

inline
bool isLittleEndian()
{
    const unsigned int i {1};
    const char* c {reinterpret_cast<const char*>(&i)};
    return static_cast<bool>(*c);
}

template <class T>
T adjust(T v1, T v2)
{
    T delta {v2-v1};
    while (delta + v1 < v2)
    {
        delta = std::nextafter(delta, std::numeric_limits<T>::infinity());
    }
    return delta;
}

template <class T, class Functor>
T adjust(T v1, T v2, const Functor& f)
{
    T delta {v2-v1};
    while (f(v1, v2, delta))
    {
        delta = std::nextafter(delta, std::numeric_limits<T>::infinity());
    }
    return delta;
}

template <class T,
          typename Kokkos::Array<T, 0>::size_type N>
Kokkos::Array<T, N> make_array(const T& value)
{
    Kokkos::Array<T, N> array{};
    for (typename Kokkos::Array<T, 0>::size_type i = 0; i < N; ++i)
    {
        array[i] = value;
    }
    return array;
}

namespace
{

template <class T,
          typename Kokkos::Array<T, 0>::size_type I,
          typename Kokkos::Array<T, 0>::size_type... J>
struct _MultiArray
{
    using Nested = typename _MultiArray<T, J...>::type;
    using type = Kokkos::Array<Nested, I>;
};

template <class T,
          typename Kokkos::Array<T, 0>::size_type I>
struct _MultiArray<T, I>
{
    using type = Kokkos::Array<T, I>;
};

}

template <class T,
          typename Kokkos::Array<T, 0>::size_type ...N>
using MultiArray = typename _MultiArray<T, N...>::type;

}  // utils
