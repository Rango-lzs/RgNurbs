#ifndef Point_HH
#define Point_HH

class Point3d
{
public:
    Point3d(double x, double y, double z);

    double x() const;
    double y() const;
    double z() const;

    Point3d operator +(const Point3d& pt) const;
private:
    double mx, my, mz;
};

Point3d operator *(const Point3d& pt, double d);
Point3d operator *(double d, const Point3d& pt);
Point3d operator /(const Point3d& pt, double d);


class Point2d
{
public:
    Point2d(double x, double y);

    double x() const;
    double y() const;

    Point2d operator +(const Point2d& pt) const;
private:
    double mx, my;
};

Point2d operator *(const Point2d& pt, double d);
Point2d operator *(double d, const Point2d& pt);
Point2d operator /(const Point2d& pt, double d);


#endif
