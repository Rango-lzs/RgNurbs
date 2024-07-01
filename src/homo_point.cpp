#include "homo_point.h"

HPoint3d::HPoint3d()
	:HPoint3d(0.0, 0.0, 0.0, 1)
{

}

HPoint3d::HPoint3d(double x, double y, double z, double w)
	: mx(x), my(y), mz(z), mw(w)
{

}

double HPoint3d::x() const
{
	return mx;
}

double HPoint3d::y() const
{
	return my;
}

double HPoint3d::z() const
{
	return mz;
}

double HPoint3d::wx() const
{
	return mx * mw;
}

double HPoint3d::wy() const
{
	return my * mw;
}

double HPoint3d::wz() const
{
	return mz * mw;
}

double HPoint3d::wi() const
{
	return mw;
}

double HPoint3d::DistanceTo(const HPoint3d& p) const
{
	return sqrt((mx - p.mx) * (mx - p.mx) + (my - p.my) * (my - p.my) +
		(mz - p.mz) * (mz - p.mz) + (mw - p.mw) * (mw - p.mw));
}

HPoint3d HPoint3d::operator +(const HPoint3d& pt) const
{
	return HPoint3d(mx + pt.mx, my + pt.my, mz + pt.mz, mw + pt.mw);
}

HPoint3d HPoint3d::operator -(const HPoint3d& pt) const
{
	return  *this + (-1.0) * pt;
}

HPoint3d& HPoint3d::operator +=(const HPoint3d& pt)
{
	*this = *this + pt;
	return *this;
}

HPoint3d& HPoint3d::operator -=(const HPoint3d& pt)
{
	*this = *this - pt;
	return *this;
}

HPoint3d HPoint3d::operator *(double d) const
{
	return HPoint3d{ mx * d, my * d, mz * d, mw * d };
}

HPoint3d HPoint3d::operator /(double d) const
{
	return *this * (1 / d);
}

HPoint3d operator *(double d, const HPoint3d& pt)
{
	return  pt * d;
}

std::ostream& operator << (std::ostream& out, const HPoint3d& pt)
{
	out << pt.x() << pt.y() << pt.z() << pt.wi() << std::endl;
	return out;
}

HPoint2d::HPoint2d()
	:HPoint2d(0.0, 0.0, 1)
{

}

HPoint2d::HPoint2d(double x, double y, double w)
	: mx(x), my(y), mw(w)
{

}

double HPoint2d::x() const
{
	return mx;
}

double HPoint2d::y() const
{
	return my;
}

double HPoint2d::wx() const
{
	return mx * mw;
}

double HPoint2d::wy() const
{
	return my * mw;
}

double HPoint2d::wi() const
{
	return mw;
}

double HPoint2d::DistanceTo(const HPoint2d& p) const
{
	return sqrt((mx - p.mx) * (mx - p.mx) + (my - p.my) * (my - p.my) +
		(mw - p.mw) * (mw - p.mw));
}

HPoint2d HPoint2d::operator +(const HPoint2d& pt) const
{
	return HPoint2d(mx + pt.mx, my + pt.my, mw + pt.mw);
}

HPoint2d HPoint2d::operator -(const HPoint2d& pt) const
{
	return  *this + (-1.0) * pt;
}

HPoint2d& HPoint2d::operator +=(const HPoint2d& pt)
{
	*this = *this + pt;
	return *this;
}

HPoint2d& HPoint2d::operator -=(const HPoint2d& pt)
{
	*this = *this - pt;
	return *this;
}

HPoint2d HPoint2d::operator *(double d) const
{
	return HPoint2d{ mx * d, my * d, mw * d };
}

HPoint2d HPoint2d::operator /(double d) const
{
	return *this * (1 / d);
}

HPoint2d operator *(double d, const HPoint2d& pt)
{
	return  pt * d;
}

std::ostream& operator << (std::ostream& out, const HPoint2d& pt)
{
	out << pt.x() << pt.y() << pt.wi() << std::endl;
	return out;
}