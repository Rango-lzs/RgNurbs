#include "beizier_curve.h"


namespace
{
	int comb(int n, int r)
	{
		if (r == 0 || r == n)
			return 1;
		else
			return comb(n - 1, r - 1) + comb(n - 1, r);
	}

	std::vector<double> calBernstein(int degree, double param)
	{
		std::vector<double> B(degree + 1);
		B[0] = 1.0;

		double t1 = 1.0 - param;
		for (int j = 1; j <= degree; j++)
		{
			double saved = 0.0;
			for (int k = 0; k < j; k++)
			{
				double temp = B[k];
				B[k] = saved + t1 * temp;
				saved = param * temp;
			}
			B[j] = saved;
		}
		return B;
	}
}

Point operator *(const Point& pt, double d)
{
	return Point{ pt.x() * d, pt.y() * d, pt.z() * d};
}

Point operator *(double d, const Point& pt)
{
	return  pt * d;
}

Point operator /(const Point& pt, double d)
{
	return pt * (1 / d);
}


BeizierCurve::BeizierCurve(int degree, std::vector<Point> ctrlPts)
	:m_degree(degree)
	,m_ctrlPts(ctrlPts)
{
	//std::assert(degree == ctrlPts.size() - 1);
}

//
Point BeizierCurve::EvalPoint(double param)
{
	std::vector<double> B = calBernstein(m_degree, param);
	Point temp(0,0,0);
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + B[k] * m_ctrlPts[k];
	}
	return temp;
}

Point BeizierCurve::EvalPointDirect(double t)
{
	Point temp(0, 0, 0);
	int n = m_degree;
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + comb(n, k) * pow(t,k) * (pow((1 - t) ,(n - k))) * m_ctrlPts[k];
	}
	return temp;
}
