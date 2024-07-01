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

double Point3d::DistanceTo(const Point3d& p) const
{
	return sqrt((mx - p.mx) * (mx - p.mx) + (my - p.my) * (my - p.my) + (mz - p.mz) * (mz - p.mz));
}

double Point3d::DotProduct(const Point3d& p) const
{
	return mx * p.mx + my * p.my + mz * p.mz;
}

double Point3d::Length() const
{
	return std::sqrt(mx * mx + my * my + mz*mz);
}

Point3d Point3d::operator +(const Point3d& pt) const
{
	return Point3d(mx + pt.mx, my + pt.my, mz + pt.mz);
}

Point3d Point3d::operator -(const Point3d& pt) const
{
	return  *this + (-1.0) * pt;
}

Point3d& Point3d::operator +=(const Point3d& pt)
{
	*this = *this + pt;
	return *this;
}

Point3d& Point3d::operator -=(const Point3d& pt)
{
	*this = *this - pt;
	return *this;
}

Point3d Point3d::operator *(double d) const
{
	return Point3d{ mx * d, my * d, mz * d };
}

Point3d Point3d::operator /(double d) const
{
	return *this * (1 / d);
}

Point3d operator *(double d, const Point3d& pt)
{
	return  pt * d;
}

std::ostream& operator << (std::ostream& out, const Point3d& pt)
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

double Point2d::DistanceTo(const Point2d& p) const
{
	return sqrt((mx - p.mx) * (mx - p.mx) + (my - p.my) * (my - p.my));
}

double Point2d::DotProduct(const Point2d& p) const
{
	return mx * p.mx + my * p.my;
}

double Point2d::Length() const
{
	return std::sqrt(mx*mx+my*my);
}

Point2d Point2d::operator +(const Point2d& pt) const
{
	return Point2d(mx + pt.mx, my + pt.my);
}

Point2d Point2d::operator -(const Point2d& pt) const
{
	return Point2d(mx - pt.mx, my - pt.my);
}

Point2d& Point2d::operator +=(const Point2d& pt)
{
	*this = *this + pt;
	return *this;
}

Point2d Point2d::operator -(const Point2d& pt) const
{
	return  *this + (-1.0) * pt;
}

Point2d& Point2d::operator -=(const Point2d& pt)
{
	*this = *this - pt;
	return *this;
}

Point2d Point2d::operator *(double d) const
{
	return Point2d{ mx * d, my * d};
}

Point2d Point2d::operator /(double d) const
{
	return *this * (1 / d);
}

std::ostream& operator << (std::ostream& out, const Point2d& pt)
{
	out << pt.x() << pt.y() << std::endl;
	return out;
}