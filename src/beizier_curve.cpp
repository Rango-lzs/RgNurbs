#include "beizier_curve.h"
#include "nurbs_export.h"
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

template<class T>
BeizierCurve<T>::BeizierCurve(int degree, std::vector<T> ctrlPts)
	:m_degree(degree)
	,m_ctrlPts(ctrlPts)
{
	assert(degree == ctrlPts.size() - 1);
}

//
template<class T>
T BeizierCurve<T>::EvalPoint(double param)
{
	std::vector<double> B = calBernstein(m_degree, param);
	T temp;
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + B[k] * m_ctrlPts[k];
	}
	return temp;
}

template<class T>
T BeizierCurve<T>::EvalPointDirect(double t)
{
	T temp;
	int n = m_degree;
	for (int k = 0; k <= m_degree; k++)
	{
		temp = temp + comb(n, k) * pow(t,k) * (pow((1 - t) ,(n - k))) * m_ctrlPts[k];
	}
	return temp;
}

template class RG_API BeizierCurve<Point3d>;
template class RG_API BeizierCurve<Point2d>;
