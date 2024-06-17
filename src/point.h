#ifndef POINT_HH
#define POINT_HH


class Point
{
public:
    Point(double x, double y, double z)
        :mx(x), my(y), mz(z)
    {

    }

    double x() const
    {
        return mx;
    }

    double y() const
    {
        return my;
    }

    double z() const
    {
        return mz;
    }

    Point operator +(const Point& pt) const
    {
        return Point(mx + pt.mx, my + pt.my, mz + pt.mz);
    }
private:
    double mx, my, mz;
};

Point operator *(const Point& pt, double d);
Point operator *(double d, const Point& pt);
Point operator /(const Point& pt, double d);

#endif
