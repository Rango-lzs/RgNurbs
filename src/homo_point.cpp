#include "homo_point.h"

HPoint3d::HPoint3d()
	:HPoint3d(0.0, 0.0, 0.0, 1)
{

}

HPoint3d::HPoint3d(double x, double y, double z, double w)
	: m_wx(w*x), m_wy(w*y), m_wz(w*z), m_w(w)
{

}

double HPoint3d::x() const
{
	return m_wx/m_w;
}

double HPoint3d::y() const
{
	return m_wy/m_w;
}

double HPoint3d::z() const
{
	return m_wz/m_w;
}

double HPoint3d::wx() const
{
	return m_wx;
}

double HPoint3d::wy() const
{
	return m_wy;
}

double HPoint3d::wz() const
{
	return m_wz;
}

double HPoint3d::wi() const
{
	return m_w;
}

HPoint3d::Point HPoint3d::ToPoint() const
{
	return HPoint3d::Point(x(), y(), z());
}

double HPoint3d::DistanceTo(const HPoint3d& p) const
{
	return sqrt((m_wx - p.m_wx) * (m_wx - p.m_wx) + (m_wy - p.m_wy) * (m_wy - p.m_wy) +
		(m_wz - p.m_wz) * (m_wz - p.m_wz) + (m_w - p.m_w) * (m_w - p.m_w));
}

HPoint3d HPoint3d::operator +(const HPoint3d& pt) const
{
	return HPoint3d(m_wx + pt.m_wx, m_wy + pt.m_wy, m_wz + pt.m_wz, m_w + pt.m_w);
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
	return HPoint3d{ m_wx * d, m_wy * d, m_wz * d, m_w * d };
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
	: m_wx(w*x), m_wy(w*y), m_w(w)
{

}

double HPoint2d::x() const
{
	return m_wx/m_w;
}

double HPoint2d::y() const
{
	return m_wy / m_w;
}

double HPoint2d::wx() const
{
	return m_wx;
}

double HPoint2d::wy() const
{
	return m_wy;
}

double HPoint2d::wi() const
{
	return m_w;
}

double HPoint2d::DistanceTo(const HPoint2d& p) const
{
	return sqrt((m_wx - p.m_wx) * (m_wx - p.m_wx) + (m_wy - p.m_wy) * (m_wy - p.m_wy) +
		(m_w - p.m_w) * (m_w - p.m_w));
}

HPoint2d HPoint2d::operator +(const HPoint2d& pt) const
{
	return HPoint2d(m_wx + pt.m_wx, m_wy + pt.m_wy, m_w + pt.m_w);
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
	return HPoint2d{ m_wx * d, m_wy * d, m_w * d };
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