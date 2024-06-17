#ifndef BEIZIER_CURVE_HH
#define BEIZIER_CURVE_HH

#include <vector>

class Point
{
public:
    Point(double x, double y, double z)
        :mx(x),my(y),mz(z)
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


class BeizierCurve
{
public:
	//beizier curve was total defined by control points
	BeizierCurve(int degree, std::vector<Point> ctrlPts);

    //
    Point EvalPoint(double param);

    Point EvalPointDirect(double param);

    Point EvalPointByDeCasteljau(double param)
    {
        std::vector<Point> temp =m_ctrlPts;
        for (int k = 1; k <= m_degree; k++)
        {
            for (int i = 0; i <= m_degree - k; i++)
            {
                temp[i] = (1.0 - param) * temp[i] + param * temp[i + 1];
            }
        }
        return temp[0];
    }

private:
    int m_degree;
    std::vector<Point> m_ctrlPts;
};

#endif
