#include "point.h"

Point3d::Point3d()
	:Point3d(0.0,0.0,0.0)
{

}

Point3d::Point3d(double x, double y, double z)
	:mx(x), my(y), mz(z)
{

}

double Point3d::x() const
{
	return mx;
}

double Point3d::y() const
{
	return my;
}

double Point3d::z() const
{
	return mz;
}

Point3d Point3d::operator +(const Point3d& pt) const
{
	return Point3d(mx + pt.mx, my + pt.my, mz + pt.mz);
}


Point3d operator *(const Point3d& pt, double d)
{
	return Point3d{ pt.x() * d, pt.y() * d, pt.z() * d };
}

Point3d operator *(double d, const Point3d& pt)
{
	return  pt * d;
}

Point3d operator /(const Point3d& pt, double d)
{
	return pt * (1 / d);
}

std::iostream& operator << (std::iostream out, const Point3d& pt)
{
	out << pt.x() << pt.y() << pt.z() << std::endl;
	return out;
}

Point2d::Point2d()
	:Point2d(0.0, 0.0)
{

}

Point2d::Point2d(double x, double y)
	:mx(x), my(y)
{

}

double Point2d::x() const
{
	return mx;
}

double Point2d::y() const
{
	return my;
}

Point2d Point2d::operator +(const Point2d& pt) const
{
	return Point2d(mx + pt.mx, my + pt.my);
}


Point2d operator *(const Point2d& pt, double d)
{
	return Point2d{ pt.x() * d, pt.y() * d};
}

Point2d operator *(double d, const Point2d& pt)
{
	return  pt * d;
}

Point2d operator /(const Point2d& pt, double d)
{
	return pt * (1 / d);
}