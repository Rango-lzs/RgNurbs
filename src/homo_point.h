#ifndef HPoint_HH
#define HPoint_HH

#include "nurbs_export.h"
#include "point.h"
#include <iostream>

class RG_API HPoint3d
{
public:
    typedef Point3d Point;

    HPoint3d();
    HPoint3d(double x, double y, double z, double w);

    double x() const;
    double y() const;
    double z() const;

    double wx() const;
    double wy() const;
    double wz() const;
    double wi() const;

    double DistanceTo(const HPoint3d& other) const;
    Point ToPoint() const;

    HPoint3d operator +(const HPoint3d& pt) const;
    HPoint3d operator -(const HPoint3d& pt) const;
    HPoint3d operator *(double d) const;
    HPoint3d operator /(double d) const;
    HPoint3d& operator +=(const HPoint3d& pt);
    HPoint3d& operator -=(const HPoint3d& pt);
private:
    double m_wx, m_wy, m_wz, m_w;
};

RG_API HPoint3d operator *(double d, const HPoint3d& pt);
RG_API std::ostream& operator << (std::ostream& out, const HPoint3d& pt);


class RG_API HPoint2d
{
public:
    typedef Point2d Point;

    HPoint2d();
    HPoint2d(double x, double y, double w);

    double x() const;
    double y() const;

    double wx() const;
    double wy() const;
    double wi() const;

    double DistanceTo(const HPoint2d& other) const;
    Point ToPoint() const;

    HPoint2d operator +(const HPoint2d& pt) const;
    HPoint2d operator -(const HPoint2d& pt) const;
    HPoint2d operator *(double d) const;
    HPoint2d operator /(double d) const;
    HPoint2d& operator +=(const HPoint2d& pt);
    HPoint2d& operator -=(const HPoint2d& pt);
private:
    double m_wx, m_wy, m_w;
};

RG_API HPoint2d operator *(double d, const HPoint2d& pt);
RG_API std::ostream& operator << (std::ostream& out, const HPoint2d& pt);

#endif
