#pragma once

#include "HydroTypes.hpp"

namespace hydro
{

class StridedFillingCurve
{
    using  IntVector =  IntVectorNd<three_d>;

public:
    StridedFillingCurve();
    StridedFillingCurve(IntVector nbCells, const IntVector& ghostWidths);
    StridedFillingCurve(IntVector nbCells, Int ghostWidth);
    StridedFillingCurve(const StridedFillingCurve& x) = default;
    StridedFillingCurve(StridedFillingCurve&& x) = default;
    ~StridedFillingCurve() = default;
    StridedFillingCurve& operator=(const StridedFillingCurve& x) = default;
    StridedFillingCurve& operator=(StridedFillingCurve&& x) = default;

    KOKKOS_INLINE_FUNCTION
    IntVector indexToCoord(Int index) const;

    KOKKOS_INLINE_FUNCTION
    Int coordToIndex(const IntVector& coord) const;

    KOKKOS_INLINE_FUNCTION
    Int faceToDim(Int face) const;

    KOKKOS_INLINE_FUNCTION
    Int faceToSide(Int face) const;

    KOKKOS_INLINE_FUNCTION
    constexpr Int faces(Int) const;

    KOKKOS_INLINE_FUNCTION
    Int getNeighbour(Int index, Int face) const;

    KOKKOS_INLINE_FUNCTION
    constexpr Int nbNeighbours() const;

    KOKKOS_INLINE_FUNCTION
    bool belongsToInnerDomain(const IntVector& coords, Int width=0) const;

    KOKKOS_INLINE_FUNCTION
    bool belongsToInnerDomain(Int index, Int width=0) const;

    KOKKOS_INLINE_FUNCTION
    Int nbCells() const;

    Int m_size;
    IntVector m_ghostWidths;
    IntVector m_nbCells;
    IntVectorNd<three_d+1> m_strides;
    IntVectorNd<2*three_d> m_neighbour_strides;
};
inline
StridedFillingCurve::StridedFillingCurve()
    : m_size        {1}
    , m_ghostWidths {}
    , m_nbCells     {}
    , m_strides     {}
    , m_neighbour_strides {}
{
}
inline
StridedFillingCurve::StridedFillingCurve(IntVector nbCells, const IntVector& ghostWidths)
    : m_size       {1}
    , m_ghostWidths {ghostWidths}
    , m_nbCells    {nbCells}
    , m_strides    {}
    , m_neighbour_strides {}
{
    for (Int idim=0; idim<three_d; ++idim)
    {
        m_size *= nbCells[idim] + 2*m_ghostWidths[idim];
    }
    for (Int idim=0; idim<=three_d; ++idim)
    {
        m_strides[idim] = (idim==0) ? 1 : m_strides[idim-1]*(nbCells[idim-1] + 2*m_ghostWidths[idim-1]);
    }
    for (Int face=0; face<2*three_d; ++face)
    {
        m_neighbour_strides[face] = faceToSide(face)*m_strides[faceToDim(face)];
    }
}
inline
StridedFillingCurve::StridedFillingCurve(IntVector nbCells, Int ghostWidth)
    : m_size       {1}
    , m_ghostWidths {}
    , m_nbCells    {nbCells}
    , m_strides    {}
    , m_neighbour_strides {}
{
    for (Int idim = 0; idim < three_d; ++idim)
    {
        m_ghostWidths[idim] = ghostWidth;
    }
    for (Int idim=0; idim<three_d; ++idim)
    {
        m_size *= nbCells[idim] + 2*m_ghostWidths[idim];
    }
    for (Int idim=0; idim<=three_d; ++idim)
    {
        m_strides[idim] = (idim==0) ? 1 : m_strides[idim-1]*(nbCells[idim-1] + 2*m_ghostWidths[idim-1]);
    }
    for (Int face=0; face<2*three_d; ++face)
    {
        m_neighbour_strides[face] = faceToSide(face)*m_strides[faceToDim(face)];
    }
}

KOKKOS_INLINE_FUNCTION
IntVectorNd<three_d> StridedFillingCurve::indexToCoord(Int index) const
{
    IntVectorNd<three_d> coord;
    coord[IZ] = index / m_strides[IZ];
    coord[IY] = (index - coord[IZ] * m_strides[IZ]) / m_strides[IY];
    coord[IX] = index - coord[IY] * m_strides[IY] - coord[IZ] * m_strides[IZ];
    return coord;
}

KOKKOS_INLINE_FUNCTION
Int StridedFillingCurve::coordToIndex(const IntVector& coord) const
{
    return coord[IX] * m_strides[IX] + coord[IY] * m_strides[IY] + coord[IZ] * m_strides[IZ];
}

KOKKOS_INLINE_FUNCTION
Int StridedFillingCurve::faceToDim(Int face) const
{
    return face / 2;
}

KOKKOS_INLINE_FUNCTION
Int StridedFillingCurve::faceToSide(Int face) const
{
    Int side = face - 2*faceToDim(face);
    return (side==0) ? -1 : 1;
}

KOKKOS_INLINE_FUNCTION
constexpr Int StridedFillingCurve::faces(Int) const
{
    return 2*three_d;
}

KOKKOS_INLINE_FUNCTION
Int StridedFillingCurve::getNeighbour(Int index, Int face) const
{
    return index + m_neighbour_strides[face];
}

KOKKOS_INLINE_FUNCTION
constexpr Int StridedFillingCurve::nbNeighbours() const
{
    return 2*three_d;
}

KOKKOS_FORCEINLINE_FUNCTION
Int StridedFillingCurve::nbCells() const
{
    return m_size;
}

KOKKOS_INLINE_FUNCTION
bool StridedFillingCurve::belongsToInnerDomain(const IntVector& coords, Int width) const
{
    bool belongsToInnerDomain_b {true};
    for (Int idim=0; idim<three_d; ++idim)
    {
        belongsToInnerDomain_b = belongsToInnerDomain_b && ((coords[idim] >= m_ghostWidths[idim]-width) &&
                                                            (coords[idim] < m_ghostWidths[idim]+m_nbCells[idim]+width));
    }
    return belongsToInnerDomain_b;
}

KOKKOS_INLINE_FUNCTION
bool StridedFillingCurve::belongsToInnerDomain(Int index, Int width) const
{
    return belongsToInnerDomain(indexToCoord(index), width);
}

}
