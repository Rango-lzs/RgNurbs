#ifndef Point_HH
#define Point_HH

#include "nurbs_export.h"

class RG_API Point3d
{
public:
    Point3d();
    Point3d(double x, double y, double z);

    double x() const;
    double y() const;
    double z() const;

    Point3d operator +(const Point3d& pt) const;
private:
    double mx, my, mz;
};

RG_API Point3d operator *(const Point3d& pt, double d);
RG_API Point3d operator *(double d, const Point3d& pt);
RG_API Point3d operator /(const Point3d& pt, double d);

class RG_API Point2d
{
public:
    Point2d();
    Point2d(double x, double y);

    double x() const;
    double y() const;

    Point2d operator +(const Point2d& pt) const;
private:
    double mx, my;
};

RG_API Point2d operator *(const Point2d& pt, double d);
RG_API Point2d operator *(double d, const Point2d& pt);
RG_API Point2d operator /(const Point2d& pt, double d);

#endif
