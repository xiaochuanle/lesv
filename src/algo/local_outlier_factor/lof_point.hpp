#ifndef __LOF_POINT_HPP
#define __LOF_POINT_HPP

#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <exception>
#include <iostream>

namespace hbn_lof {

using std::cerr;
using std::endl;
using std::exception;

template <int DIM>
class LofPoint
{
public:
    LofPoint() {}
    LofPoint(double x0, ...);
    LofPoint(const LofPoint<DIM>& pnt);
    LofPoint<DIM>& operator=(const LofPoint<DIM>& pnt);
    ~LofPoint() {}
    double operator[](int i) const;
    double& operator[](int i);
    bool operator==(const LofPoint<DIM>& y) const;
    bool operator!=(const LofPoint<DIM>& y) const;

private:
    std::array<double, DIM> m_coord;
};

template <int DIM>
bool LofPoint<DIM>::operator==(const LofPoint<DIM>& y) const
{
    for (int i = 0; i < DIM; ++i) {
        if (fabs(m_coord[i]-y[i]) > 1.0e-6) return false;
    }
    return true;
}

template <int DIM>
bool LofPoint<DIM>::operator!=(const LofPoint<DIM>& y) const
{
    return !((*this) == y);
}

template <int DIM>
double LofPoint<DIM>::operator[](int i) const
{
    return m_coord.at(i);
}

template <int DIM>
double& LofPoint<DIM>::operator[](int i)
{
    return m_coord.at(i);
}

template <int DIM>
LofPoint<DIM>::LofPoint(double x0, ...)
{
    m_coord.at(0) = x0;
    va_list ap;
    va_start(ap, x0);
    for (int i = 1; i < DIM; ++i) {
        double x = va_arg(ap, double);
        m_coord.at(i) = x;
    }
    va_end(ap);
}

template <int DIM>
LofPoint<DIM>::LofPoint(const LofPoint<DIM>& pnt)
{
    for (int i = 0; i < DIM; ++i) m_coord.at(i) = pnt[i];
}

template <int DIM>
LofPoint<DIM>& LofPoint<DIM>::operator=(const LofPoint<DIM>& pnt)
{
    for (int i = 0; i < DIM; ++i) m_coord[i] = pnt[i];
    return *this;
}

template <int DIM>
std::ostream& operator<<(std::ostream& os, const LofPoint<DIM>& pnt)
{
    os << "(";
    for (int i = 0; i < DIM - 1; ++i) os << pnt[i] << ", ";
    os << pnt[DIM-1] <<  ")";
    return os;
}

} // hbn_lof

#endif // __LOF_POINT_HPP