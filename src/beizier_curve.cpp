#include "beizier_curve.h"
#include <assert.h>


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

BeizierCurve::BeizierCurve(int degree, std::vector<Point3d> ctrlPts)
	:m_degree(degree)
	,m_ctrlPts(ctrlPts)
{
	assert(degree == ctrlPts.size() - 1);
}

//
Point3d BeizierCurve::EvalPoint(double param)
{
	std::vector<double> B = calBernstein(m_degree, param);
	Point3d temp(0,0,0);
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + B[k] * m_ctrlPts[k];
	}
	return temp;
}

Point3d BeizierCurve::EvalPointDirect(double t)
{
	Point3d temp(0, 0, 0);
	int n = m_degree;
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + comb(n, k) * pow(t,k) * (pow((1 - t) ,(n - k))) * m_ctrlPts[k];
	}
	return temp;
}
