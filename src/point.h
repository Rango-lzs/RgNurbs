#ifndef Point_HH
#define Point_HH

#include "nurbs_export.h"
#include <iostream>

class RG_API Point3d
{
public:
    Point3d();
    Point3d(double x, double y, double z);

    double x() const;
    double y() const;
    double z() const;

    double DistanceTo(const Point3d& other) const;
    double DotProduct(const Point3d& other) const;
    double Length() const;

    Point3d operator +(const Point3d& pt) const;
    Point3d operator -(const Point3d& pt) const;
    Point3d operator *(double d) const;
    Point3d operator /(double d) const;
    Point3d& operator +=(const Point3d& pt);
    Point3d& operator -=(const Point3d& pt);
    
private:
    double mx, my, mz;
};

RG_API Point3d operator *(double d, const Point3d& pt);
RG_API std::ostream& operator << (std::ostream& out, const Point3d& pt);


class RG_API Point2d
{
public:
    Point2d();
    Point2d(double x, double y);

    double x() const;
    double y() const;

    double DistanceTo(const Point2d& other) const;
    double DotProduct(const Point2d& other) const;
    double Length() const;

    Point2d operator +(const Point2d& pt) const;
    Point2d operator -(const Point2d& pt) const;
    Point2d operator *(double d) const;
    Point2d operator /(double d) const;
    Point2d& operator +=(const Point2d& pt);
    Point2d& operator -=(const Point2d& pt);
private:
    double mx, my;
};

RG_API Point2d operator *(double d, const Point2d& pt);
RG_API std::ostream& operator << (std::ostream& out, const Point2d& pt);

#endif
