//#include "bspline_curve.h"
//
//
//std::vector<double> basis_funcs(int span, int degree, const std::vector<double>& knots, double t)
//{
//	std::vector<double> N(degree + 1);
//	N[0] = 1.0;
//
//	std::vector<double> left(degree + 1);
//	std::vector<double> right(degree + 1);
//
//	for (int j = 1; j <= degree; j++)
//	{
//		left[j] = t - m_knots[span + 1 - j];
//		right[j] = m_knots[span + j] - t;
//
//		double saved = 0.0;
//
//		for (int r = 0; r < j; r++)
//		{
//			double temp = N[r] / (right[r + 1] + left[j - r]);
//			N[r] = saved + right[r + 1] * temp;
//			saved = left[j - r] * temp;
//		}
//		N[j] = saved;
//	}
//	return N;
//}

